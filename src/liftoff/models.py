"""Core data structures shared across the pipeline.

These replace the attribute-bag classes of earlier versions (``new_feature``,
``aligned_seg`` and ``feature_hierarchy``) with typed, slotted dataclasses.
Instances are intentionally mutable: the lift-over algorithm adjusts feature
coordinates and attributes in place as it refines mappings.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum


class LiftoverType(StrEnum):
    """Stages of the lift-over pipeline.

    The value of each member is used to name intermediate files, so it must
    remain stable.
    """

    CHROM_BY_CHROM = 'chrm_by_chrm'
    UNMAPPED = 'unmapped'
    UNPLACED = 'unplaced'
    COPIES = 'copies'


#: GFF attributes are multi-valued: each key maps to a list of strings.
Attributes = dict[str, list[str]]


@dataclass(slots=True, eq=False)
class Feature:
    """A single annotation feature (gene, transcript, exon, CDS, ...).

    Equality is identity-based (``eq=False``) because the algorithm tracks
    specific feature objects in lists and compares them with ``==``/``!=``.

    Attributes
    ----------
    id : str
        Unique feature identifier.
    featuretype : str
        Feature type from column 3 of the GFF (e.g. ``"gene"``).
    seqid : str
        Name of the sequence (chromosome/contig) the feature is on.
    source : str
        Annotation source from column 2.
    strand : str
        ``"+"``, ``"-"`` or ``"."``.
    start : int
        One-based inclusive start coordinate.
    end : int
        One-based inclusive end coordinate.
    frame : str
        CDS phase from column 8 (``"0"``, ``"1"``, ``"2"`` or ``"."``).
    attributes : dict of str to list of str
        Column 9 key/value attributes.
    score : float, default 0.0
        Mapping score; for lifted top-level features this is
        ``1 - sequence_identity`` (lower is better).
    original_bounds : tuple of int or None, default None
        Annotated ``(start, end)`` of a reference top-level feature whose
        ``start``/``end`` were widened by flanking sequence for alignment.
    """

    id: str
    featuretype: str
    seqid: str
    source: str
    strand: str
    start: int
    end: int
    frame: str
    attributes: Attributes
    score: float = 0.0
    original_bounds: tuple[int, int] | None = None

    @property
    def annotated_bounds(self) -> tuple[int, int]:
        """Annotated coordinates, ignoring any alignment flank.

        Returns
        -------
        tuple of int
            ``(start, end)`` as written in the reference annotation.
        """
        return self.original_bounds or (self.start, self.end)


@dataclass(slots=True, eq=False)
class AlignedSegment:
    """A gap-free block of an alignment between a reference gene and the target.

    Coordinates on the query (reference gene sequence) are zero-based and
    relative to the start of the extracted gene sequence; coordinates on the
    target are zero-based positions in the target sequence. Both ends are
    inclusive.

    Attributes
    ----------
    aln_id : int
        Identifier shared by all blocks that come from the same SAM record.
        Negative values are reserved for graph sentinel nodes.
    query_name : str
        Name of the aligned reference feature, including its copy suffix
        (e.g. ``"gene1_0"``).
    reference_name : str
        Name of the target sequence the block aligns to.
    query_block_start : int
        First aligned query position.
    query_block_end : int
        Last aligned query position.
    reference_block_start : int
        First aligned target position.
    reference_block_end : int
        Last aligned target position.
    is_reverse : bool
        Whether the query aligned to the reverse strand of the target.
    mismatches : list of int
        Sorted query positions of mismatched bases within the block.
    """

    aln_id: int
    query_name: str
    reference_name: str
    query_block_start: int
    query_block_end: int
    reference_block_start: int
    reference_block_end: int
    is_reverse: bool
    mismatches: list[int] = field(default_factory=list)


@dataclass(slots=True)
class FeatureHierarchy:
    """Reference features grouped by their position in the feature tree.

    Attributes
    ----------
    parents : dict of str to Feature
        Top-level features to lift over, keyed by ID.
    intermediates : dict of str to Feature
        Features that are both a child and a parent (e.g. transcripts).
    children : dict of str to list of Feature
        Lowest-level features (e.g. exons, CDS) keyed by the ID of their
        top-level parent.
    """

    parents: dict[str, Feature]
    intermediates: dict[str, Feature]
    children: dict[str, list[Feature]]
