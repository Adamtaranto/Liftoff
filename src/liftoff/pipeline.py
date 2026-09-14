"""The end-to-end Liftoff pipeline.

Stages, in order:

1. **Chromosome-by-chromosome** (or whole genome) lift-over of all genes.
2. **Unmapped**: genes that failed a per-chromosome lift-over are aligned to
   the whole target genome.
3. **Unplaced**: genes on unplaced reference sequences are lifted.
4. **Copies** (optional): additional gene copies are searched for.
5. **CDS check / polishing**: ORFs are validated, and optionally broken genes
   are re-aligned exon by exon.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass
from pathlib import Path

import gffutils

from liftoff import cds, polish
from liftoff.align.base import Aligner
from liftoff.align.blocks import parse_all_sam_files
from liftoff.align.external import PrecomputedSamAligner
from liftoff.align.minimap2 import Minimap2Aligner
from liftoff.align.runner import align_features_to_target
from liftoff.config import STDOUT, LiftoffConfig
from liftoff.io.fasta import open_fasta, write_gene_sequences
from liftoff.io.gff_db import annotation_format, build_database, separate_parents_and_children
from liftoff.io.gff_writer import write_new_gff
from liftoff.log import get_logger
from liftoff.mapping.lift import lift_all_features
from liftoff.mapping.overlaps import fix_incorrectly_overlapping_features
from liftoff.models import Feature, FeatureHierarchy, LiftoverType
from liftoff.utils import LiftedFeatures, ParentOrder, clear_scores

logger = get_logger(__name__)

#: Relaxed thresholds used where partial mappings are reported rather than excluded.
PARTIAL_MAPPING_THRESHOLD = 0.05

#: Distance factor used when re-lifting polished genes (effectively unbounded).
POLISH_DISTANCE_FACTOR = 100_000_000


@dataclass(slots=True)
class LiftoffResult:
    """Outcome of a Liftoff run.

    Attributes
    ----------
    lifted_features : dict of str to list of Feature
        Lifted features keyed by copy ID, top-level feature first.
    unmapped_features : list of Feature
        Reference top-level features that could not be lifted.
    output_paths : list of str
        Annotation files written (``"stdout"`` when written to the terminal).
    """

    lifted_features: LiftedFeatures
    unmapped_features: list[Feature]
    output_paths: list[str]


@dataclass(slots=True)
class _Context:
    """State shared by all pipeline stages."""

    config: LiftoffConfig
    aligner: Aligner
    feature_db: gffutils.FeatureDB
    hierarchy: FeatureHierarchy
    parent_order: ParentOrder
    lifted: LiftedFeatures


def default_aligner(config: LiftoffConfig) -> Aligner:
    """Choose the alignment backend implied by the configuration.

    Parameters
    ----------
    config : LiftoffConfig
        Run configuration.

    Returns
    -------
    Aligner
        :class:`PrecomputedSamAligner` when ``config.alignments`` is set,
        otherwise :class:`Minimap2Aligner`.
    """
    if config.alignments is not None:
        return PrecomputedSamAligner(Path(config.alignments))
    return Minimap2Aligner(config.minimap2 or 'minimap2')


def run_liftoff(config: LiftoffConfig, aligner: Aligner | None = None) -> LiftoffResult:
    """Run the complete lift-over.

    Parameters
    ----------
    config : LiftoffConfig
        Run configuration.
    aligner : Aligner or None, default None
        Alignment backend; see :func:`default_aligner`.

    Returns
    -------
    LiftoffResult
        Lifted and unmapped features and the paths written.
    """
    aligner = aligner or default_aligner(config)
    if config.chroms is not None:
        ref_chroms, target_chroms = parse_chrom_file(config.chroms)
    else:
        ref_chroms, target_chroms = [config.reference], [config.target]
    parent_types = (
        None if config.all_feature_types else get_parent_features_to_lift(config.feature_types)
    )

    config.intermediate_path.mkdir(parents=True, exist_ok=True)
    logger.info('extracting features')
    feature_db = build_database(config.gff, config.db, config.infer_genes, config.infer_transcripts)
    try:
        hierarchy, parent_order = separate_parents_and_children(
            feature_db, parent_types, config.exclude_feature_types
        )
        ctx = _Context(config, aligner, feature_db, hierarchy, parent_order, {})
        return _run_stages(ctx, ref_chroms, target_chroms)
    finally:
        feature_db.conn.close()


def _run_stages(ctx: _Context, ref_chroms: list[str], target_chroms: list[str]) -> LiftoffResult:
    """Run every pipeline stage and write the outputs."""
    config = ctx.config
    unmapped: list[Feature] = []
    _lift_original_annotation(ctx, ref_chroms, target_chroms, unmapped)
    unmapped = _map_unmapped_features(ctx, target_chroms, unmapped)
    _map_features_from_unplaced_seq(ctx, unmapped)
    write_unmapped_features_file(config.unmapped, unmapped)
    if config.copies:
        _map_extra_copies(ctx)

    out_type = annotation_format(ctx.feature_db)
    outputs: list[str] = []
    if config.polish:
        logger.info('polishing annotations')
        cds.check_cds(ctx.lifted, ctx.hierarchy, config.reference, config.target)
        _write_output(ctx, config.output, out_type)
        outputs.append(config.output)
        _find_and_polish_broken_cds(ctx, unmapped)
        final_output = config.output if config.output == STDOUT else config.output + '_polished'
    else:
        if config.cds:
            cds.check_cds(ctx.lifted, ctx.hierarchy, config.reference, config.target)
        final_output = config.output
    _write_output(ctx, final_output, out_type)
    outputs.append(final_output)
    return LiftoffResult(ctx.lifted, unmapped, outputs)


def parse_chrom_file(chroms_file: str) -> tuple[list[str], list[str]]:
    """Read a ``reference,target`` chromosome correspondence file.

    Parameters
    ----------
    chroms_file : str
        File with one ``reference[,target]`` pair per line.

    Returns
    -------
    ref_chroms : list of str
        Reference chromosome names.
    target_chroms : list of str
        Target chromosome names (lines without a target are skipped).
    """
    ref_chroms: list[str] = []
    target_chroms: list[str] = []
    with Path(chroms_file).open() as handle:
        for line in handle:
            if not line.strip():
                continue
            names = line.rstrip().split(',')
            ref_chroms.append(names[0])
            if len(names) > 1:
                target_chroms.append(names[1])
    return ref_chroms, target_chroms


def get_parent_features_to_lift(feature_types_file: str | None) -> list[str]:
    """Return the top-level feature types to lift.

    Parameters
    ----------
    feature_types_file : str or None
        File listing extra types, one per line.

    Returns
    -------
    list of str
        ``"gene"`` followed by any listed types.
    """
    feature_types = ['gene']
    if feature_types_file is not None:
        with Path(feature_types_file).open() as handle:
            feature_types.extend(line.rstrip() for line in handle)
    return feature_types


def write_unmapped_features_file(path: str, unmapped_features: list[Feature]) -> None:
    """Write the IDs of unmapped features, one per line.

    Parameters
    ----------
    path : str
        Destination file.
    unmapped_features : list of Feature
        Features to list.
    """
    with Path(path).open('w') as handle:
        handle.writelines(feature.id + '\n' for feature in unmapped_features)


def _write_output(ctx: _Context, output: str, out_type: str) -> None:
    write_new_gff(
        ctx.lifted,
        output,
        out_type,
        ctx.config.min_coverage,
        ctx.config.min_identity,
        ctx.config.command_line,
    )


def _align_and_lift_features(
    ctx: _Context,
    ref_chroms: list[str],
    target_chroms: list[str],
    liftover_type: LiftoverType,
    unmapped_features: list[Feature],
    min_cov: float,
    min_seqid: float,
) -> None:
    """Align, lift and resolve overlaps for one stage."""
    config = ctx.config
    aligned_segments = align_features_to_target(
        ref_chroms,
        target_chroms,
        config.reference,
        config.target,
        config.intermediate_path,
        config.minimap2_arguments,
        config.threads,
        ctx.aligner,
        ctx.hierarchy,
        liftover_type,
        unmapped_features,
    )
    logger.info('lifting features')
    lift_all_features(
        aligned_segments,
        min_cov,
        ctx.feature_db,
        ctx.hierarchy,
        unmapped_features,
        ctx.lifted,
        min_seqid,
        None,
        config,
        ctx.parent_order,
    )
    fix_incorrectly_overlapping_features(
        ctx.lifted,
        ctx.lifted,
        aligned_segments,
        unmapped_features,
        min_cov,
        ctx.hierarchy,
        ctx.feature_db,
        ctx.parent_order,
        min_seqid,
        config,
        config.max_overlap,
    )


def _extract_sequences(
    ctx: _Context,
    parents: dict[str, Feature],
    ref_chroms: list[str],
    liftover_type: LiftoverType,
) -> None:
    write_gene_sequences(
        parents,
        ref_chroms,
        ctx.config.reference,
        ctx.config.intermediate_path,
        liftover_type,
        ctx.config.flank,
    )


def _lift_original_annotation(
    ctx: _Context, ref_chroms: list[str], target_chroms: list[str], unmapped: list[Feature]
) -> None:
    """Stage 1: lift all genes, per chromosome when a correspondence is given."""
    config = ctx.config
    if target_chroms[0] == config.target and not config.exclude_partial:
        min_cov = min_seqid = PARTIAL_MAPPING_THRESHOLD
    else:
        min_cov, min_seqid = config.min_coverage, config.min_identity
    _extract_sequences(ctx, ctx.hierarchy.parents, ref_chroms, LiftoverType.CHROM_BY_CHROM)
    _align_and_lift_features(
        ctx, ref_chroms, target_chroms, LiftoverType.CHROM_BY_CHROM, unmapped, min_cov, min_seqid
    )


def _stage_thresholds(config: LiftoffConfig) -> tuple[float, float]:
    """Thresholds for the unmapped and unplaced stages."""
    if config.exclude_partial:
        return config.min_coverage, config.min_identity
    return PARTIAL_MAPPING_THRESHOLD, PARTIAL_MAPPING_THRESHOLD


def _map_unmapped_features(
    ctx: _Context, target_chroms: list[str], unmapped: list[Feature]
) -> list[Feature]:
    """Stage 2: align genes that failed per-chromosome lift-over to the whole genome."""
    config = ctx.config
    if not unmapped or target_chroms[0] == config.target:
        return unmapped
    logger.info('mapping unaligned features to whole genome')
    clear_scores(ctx.lifted, ctx.hierarchy.parents)
    unmapped_dict = {feature.id: feature for feature in unmapped}
    min_cov, min_seqid = _stage_thresholds(config)
    ref_chroms, genome = [config.reference], [config.target]
    _extract_sequences(ctx, unmapped_dict, ref_chroms, LiftoverType.UNMAPPED)
    still_unmapped: list[Feature] = []
    _align_and_lift_features(
        ctx, ref_chroms, genome, LiftoverType.UNMAPPED, still_unmapped, min_cov, min_seqid
    )
    return still_unmapped


def _map_features_from_unplaced_seq(ctx: _Context, unmapped: list[Feature]) -> None:
    """Stage 3: lift genes located on unplaced reference sequences."""
    config = ctx.config
    if config.unplaced is None or config.chroms is None:
        return
    logger.info('mapping unplaced genes')
    ref_chroms, _ = parse_chrom_file(config.unplaced)
    clear_scores(ctx.lifted, ctx.hierarchy.parents)
    unplaced_dict = {
        feature.id: feature
        for feature in ctx.hierarchy.parents.values()
        if feature.seqid in ref_chroms
    }
    _extract_sequences(ctx, unplaced_dict, ref_chroms, LiftoverType.UNPLACED)
    min_cov, min_seqid = _stage_thresholds(config)
    _align_and_lift_features(
        ctx, ref_chroms, [config.target], LiftoverType.UNPLACED, unmapped, min_cov, min_seqid
    )


def _map_extra_copies(ctx: _Context) -> None:
    """Stage 4: search the whole target genome for additional gene copies."""
    config = ctx.config
    logger.info('mapping gene copies')
    clear_scores(ctx.lifted, ctx.hierarchy.parents)
    ref_chroms, genome = [config.reference], [config.target]
    _extract_sequences(ctx, ctx.hierarchy.parents, ref_chroms, LiftoverType.COPIES)
    _align_and_lift_features(
        ctx, ref_chroms, genome, LiftoverType.COPIES, [], 0, config.copy_identity
    )


def _find_and_polish_broken_cds(ctx: _Context, unmapped: list[Feature]) -> None:
    """Stage 5: re-align exons of genes with broken ORFs and keep improvements."""
    config = ctx.config
    polish_config = dataclasses.replace(config, distance_factor=POLISH_DISTANCE_FACTOR)
    intermediate_dir = config.intermediate_path
    polished: LiftedFeatures = {}
    ref_fasta, target_fasta = open_fasta(config.reference), open_fasta(config.target)
    matrix = polish.make_scoring_matrix(3)
    for target_feature in list(ctx.lifted):
        if not polish.polish_annotations(
            ctx.lifted,
            ref_fasta,
            target_fasta,
            intermediate_dir,
            ctx.hierarchy,
            target_feature,
            matrix,
        ):
            continue
        aligned_segments = parse_all_sam_files(
            ctx.hierarchy,
            unmapped,
            LiftoverType.CHROM_BY_CHROM,
            [intermediate_dir / polish.POLISH_SAM_NAME],
        )
        if not aligned_segments:
            continue
        # All records in polish.sam belong to this gene copy.
        blocks = next(iter(aligned_segments.values()))
        for block in blocks:
            block.query_name = target_feature
        lift_all_features(
            {target_feature: blocks},
            config.min_coverage,
            ctx.feature_db,
            ctx.hierarchy,
            unmapped,
            polished,
            config.min_identity,
            None,
            polish_config,
            ctx.parent_order,
        )

    cds.check_cds(polished, ctx.hierarchy, config.reference, config.target)
    for feature_id, polished_features in polished.items():
        original = ctx.lifted[feature_id][0].attributes
        candidate = polished_features[0].attributes
        if polish.polished_is_better(original, candidate):
            ctx.lifted[feature_id] = polished_features
