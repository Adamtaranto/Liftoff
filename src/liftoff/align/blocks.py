"""Convert SAM alignments of gene sequences into gap-free aligned blocks.

Each SAM record aligns an extracted reference gene (the *query*) to the
target genome. The record's CIGAR is split at every insertion and deletion into
:class:`~liftoff.models.AlignedSegment` blocks; blocks that do not touch any
exon/CDS of the gene are discarded. The mapping stage later chains blocks into
the best-scoring gene model.
"""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

from liftoff.io.sam import CigarOp, SamRecord, read_sam
from liftoff.models import AlignedSegment, Feature, FeatureHierarchy, LiftoverType
from liftoff.utils import (
    convert_id_to_original,
    count_overlap,
    get_relative_child_coord,
    merge_children_intervals,
)

#: Aligned blocks keyed by copy-aware feature name (e.g. ``"gene1_0"``).
AlignedSegments = dict[str, list[AlignedSegment]]


def parse_all_sam_files(
    feature_hierarchy: FeatureHierarchy,
    unmapped_features: list[Feature],
    liftover_type: LiftoverType,
    sam_files: Iterable[str | Path],
) -> AlignedSegments:
    """Parse several SAM files into a single block dictionary.

    Parameters
    ----------
    feature_hierarchy : FeatureHierarchy
        Reference features.
    unmapped_features : list of Feature
        Receives top-level features without any usable alignment.
    liftover_type : LiftoverType
        Pipeline stage; ``COPIES`` numbers every alignment as a separate copy.
    sam_files : iterable of str or pathlib.Path
        SAM files in processing order.

    Returns
    -------
    dict of str to list of AlignedSegment
        Blocks for every aligned feature copy.
    """
    aligned_segments: AlignedSegments = {}
    for sam_file in sam_files:
        aligned_segments.update(
            parse_alignment(sam_file, feature_hierarchy, unmapped_features, liftover_type)
        )
    return aligned_segments


def parse_alignment(
    sam_file: str | Path,
    feature_hierarchy: FeatureHierarchy,
    unmapped_features: list[Feature],
    liftover_type: LiftoverType,
) -> AlignedSegments:
    """Parse one SAM file into aligned blocks.

    Parameters
    ----------
    sam_file : str or pathlib.Path
        SAM file of gene sequences aligned to the target.
    feature_hierarchy : FeatureHierarchy
        Reference features.
    unmapped_features : list of Feature
        Receives features whose alignments are unmapped or cover no children.
    liftover_type : LiftoverType
        Pipeline stage.

    Returns
    -------
    dict of str to list of AlignedSegment
        Blocks keyed by copy-aware feature name.
    """
    all_aligned_blocks: AlignedSegments = {}
    copy_counts: dict[str, int] = {}
    aln_id = 0
    for record in read_sam(sam_file):
        if record.is_unmapped:
            unmapped_features.append(feature_hierarchy.parents[record.query_name])
            continue
        record.query_name = _copy_aware_name(liftover_type, record.query_name, copy_counts)
        aln_id += 1
        blocks = get_aligned_blocks(record, aln_id, feature_hierarchy, liftover_type)
        all_aligned_blocks.setdefault(record.query_name, []).extend(blocks)
    _remove_alignments_without_children(all_aligned_blocks, unmapped_features, feature_hierarchy)
    return all_aligned_blocks


def _copy_aware_name(liftover_type: LiftoverType, name: str, copy_counts: dict[str, int]) -> str:
    """Suffix a query name with its copy number.

    Outside the copies stage every alignment is copy 0; when searching for
    copies each alignment of a feature receives the next number from 1.
    """
    if liftover_type != LiftoverType.COPIES:
        return name + '_0'
    copy_counts[name] = copy_counts.get(name, 0) + 1
    return f'{name}_{copy_counts[name]}'


def get_aligned_blocks(
    alignment: SamRecord,
    aln_id: int,
    feature_hierarchy: FeatureHierarchy,
    liftover_type: LiftoverType,
) -> list[AlignedSegment]:
    """Split one alignment into gap-free blocks overlapping child features.

    Parameters
    ----------
    alignment : SamRecord
        A mapped record whose ``query_name`` already carries a copy suffix.
    aln_id : int
        Identifier assigned to all blocks of this record.
    feature_hierarchy : FeatureHierarchy
        Reference features.
    liftover_type : LiftoverType
        Pipeline stage; in the copies stage only end-to-end alignments count.

    Returns
    -------
    list of AlignedSegment
        Blocks in query order.
    """
    cigar = alignment.cigartuples
    original_name = convert_id_to_original(alignment.query_name)
    parent = feature_hierarchy.parents[original_name]
    children = feature_hierarchy.children[original_name]

    query_start = alignment.query_alignment_start
    query_end = alignment.query_alignment_end
    if cigar[0][0] == CigarOp.HARD_CLIP:
        # Hard-clipped bases are absent from SEQ; shift into gene coordinates.
        query_start += cigar[0][1]
        query_end += cigar[0][1]

    end_to_end = parent.end - parent.start + 1 == query_end - query_start
    if liftover_type == LiftoverType.COPIES and not end_to_end:
        return []

    merged_children_coords = merge_children_intervals(children)
    builder = _BlockBuilder(
        alignment, aln_id, parent, merged_children_coords, query_start, alignment.reference_start
    )
    for operation, length in cigar:
        if operation in (CigarOp.EQUAL, CigarOp.DIFF):
            builder._advance_aligned(operation, length)
            if builder.query_pos == query_end:
                builder._close_block()
                break
        elif operation in (CigarOp.INSERTION, CigarOp.DELETION):
            builder._close_block()
            builder._skip_gap(operation, length)
        # Other operations (M, S, H, N, P) do not move the block cursors.
    return builder.blocks


class _BlockBuilder:
    """Cursor over a CIGAR string that emits blocks at every gap.

    Parameters
    ----------
    alignment : SamRecord
        Record being split.
    aln_id : int
        Identifier assigned to every emitted block.
    parent : Feature
        Reference top-level feature.
    merged_children_coords : list of list of int
        Merged child feature intervals of ``parent``.
    query_start : int
        Query position of the first aligned base.
    reference_start : int
        Target position of the first aligned base.
    """

    def __init__(
        self,
        alignment: SamRecord,
        aln_id: int,
        parent: Feature,
        merged_children_coords: list[list[int]],
        query_start: int,
        reference_start: int,
    ) -> None:
        self.alignment = alignment
        self.aln_id = aln_id
        self.parent = parent
        self.children_coords = merged_children_coords
        self.query_block_start = self.query_pos = query_start
        self.reference_block_start = self.reference_pos = reference_start
        self.mismatches: list[int] = []
        self.blocks: list[AlignedSegment] = []

    def _advance_aligned(self, operation: int, length: int) -> None:
        """Consume ``length`` aligned bases, recording mismatch positions."""
        if operation == CigarOp.DIFF:
            self.mismatches.extend(range(self.query_pos, self.query_pos + length))
        self.query_pos += length
        self.reference_pos += length

    def _skip_gap(self, operation: int, length: int) -> None:
        """Consume a gap and start a new block after it."""
        if operation == CigarOp.INSERTION:
            self.query_pos += length
        else:
            self.reference_pos += length
        self.mismatches = []
        self.query_block_start = self.query_pos
        self.reference_block_start = self.reference_pos

    def _close_block(self) -> None:
        """Emit the current block if it overlaps any child feature."""
        block = AlignedSegment(
            aln_id=self.aln_id,
            query_name=self.alignment.query_name,
            reference_name=self.alignment.reference_name or '',
            query_block_start=self.query_block_start,
            query_block_end=self.query_pos - 1,
            reference_block_start=self.reference_block_start,
            reference_block_end=self.reference_pos - 1,
            is_reverse=self.alignment.is_reverse,
            mismatches=list(self.mismatches),
        )
        if _overlaps_children(block, self.children_coords, self.parent):
            self.blocks.append(block)


def _overlaps_children(
    aln: AlignedSegment, children_coords: list[list[int]], parent: Feature
) -> bool:
    """Return whether a block overlaps any merged child interval."""
    for child_start_abs, child_end_abs in children_coords:
        relative_start = get_relative_child_coord(parent, child_start_abs, aln.is_reverse)
        relative_end = get_relative_child_coord(parent, child_end_abs, aln.is_reverse)
        child_start, child_end = (
            min(relative_start, relative_end),
            max(relative_start, relative_end),
        )
        if count_overlap(child_start, child_end, aln.query_block_start, aln.query_block_end) > 0:
            return True
    return False


def _remove_alignments_without_children(
    all_aligned_blocks: AlignedSegments,
    unmapped_features: list[Feature],
    feature_hierarchy: FeatureHierarchy,
) -> None:
    """Drop features with no child-overlapping blocks and mark them unmapped."""
    for name in [name for name, blocks in all_aligned_blocks.items() if not blocks]:
        unmapped_features.append(feature_hierarchy.parents[convert_id_to_original(name)])
        del all_aligned_blocks[name]
