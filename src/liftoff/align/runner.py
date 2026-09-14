"""Plan and dispatch the alignments needed by a pipeline stage."""

from __future__ import annotations

from collections.abc import Sequence
import math
from pathlib import Path
import sys

from liftoff.align.base import Aligner, AlignmentJob
from liftoff.align.blocks import AlignedSegments, parse_all_sam_files
from liftoff.io.fasta import features_file_name, get_genome_size, split_target_sequence
from liftoff.log import get_logger
from liftoff.models import Feature, FeatureHierarchy, LiftoverType

logger = get_logger(__name__)

#: Targets larger than this are aligned without a single pre-built index.
MAX_SINGLE_INDEX_SIZE = 4_000_000_000


def can_use_multiprocessing() -> bool:
    """Return whether worker processes can be started on this platform.

    Returns
    -------
    bool
        ``False`` under Emscripten/WASI (Pyodide) or when the multiprocessing
        machinery is unavailable.
    """
    if sys.platform in ('emscripten', 'wasi'):
        return False
    try:
        import multiprocessing.synchronize  # noqa: F401  # needs working semaphores
    except ImportError:  # pragma: no cover - platform specific
        return False
    return True


def plan_alignment_jobs(
    ref_chroms: Sequence[str],
    target_chroms: Sequence[str],
    reference: str,
    target: str,
    intermediate_dir: Path,
    minimap2_options: Sequence[str],
    threads: int,
    liftover_type: LiftoverType,
) -> list[AlignmentJob]:
    """Describe the alignments for one pipeline stage.

    One job is created per target chromosome (or a single job for the whole
    genome). Per-chromosome target FASTA files are written as a side effect.

    Parameters
    ----------
    ref_chroms : sequence of str
        Reference chromosomes, parallel to ``target_chroms``.
    target_chroms : sequence of str
        Target chromosomes; ``[target]`` means the whole genome.
    reference : str
        Reference genome FASTA.
    target : str
        Target genome FASTA.
    intermediate_dir : pathlib.Path
        Directory for intermediate files.
    minimap2_options : sequence of str
        Arguments passed to minimap2.
    threads : int
        Total threads available.
    liftover_type : LiftoverType
        Pipeline stage.

    Returns
    -------
    list of AlignmentJob
        Jobs in chromosome order.
    """
    target_fasta = split_target_sequence(target_chroms, target, intermediate_dir)
    genome_size = get_genome_size(target_fasta)
    threads_per_alignment = max(1, math.floor(threads / len(ref_chroms)))
    whole_genome = liftover_type != LiftoverType.CHROM_BY_CHROM or target_chroms[0] == target
    jobs = []
    for index in range(len(target_chroms)):
        features_name = features_file_name(ref_chroms[index], reference, liftover_type)
        if whole_genome:
            target_file = Path(target)
            target_name = 'target_all'
        else:
            target_file = intermediate_dir / f'{target_chroms[index]}.fa'
            target_name = target_chroms[index]
        split_prefix = None
        if genome_size > MAX_SINGLE_INDEX_SIZE:
            split_prefix = intermediate_dir / f'{features_name}_to_{target_name}_split'
        jobs.append(
            AlignmentJob(
                features_fasta=intermediate_dir / f'{features_name}_genes.fa',
                target_fasta=target_file,
                output_sam=intermediate_dir / f'{features_name}_to_{target_name}.sam',
                minimap2_options=tuple(minimap2_options),
                threads=threads_per_alignment,
                liftover_type=liftover_type,
                split_prefix=split_prefix,
            )
        )
    return jobs


def run_alignment_jobs(
    aligner: Aligner, jobs: Sequence[AlignmentJob], processes: int
) -> list[Path]:
    """Execute alignment jobs, in parallel when possible.

    Parameters
    ----------
    aligner : Aligner
        Backend performing each alignment.
    jobs : sequence of AlignmentJob
        Jobs to run.
    processes : int
        Maximum number of worker processes.

    Returns
    -------
    list of pathlib.Path
        SAM files, in the same order as ``jobs`` (so results are
        deterministic regardless of completion order).
    """
    for job in jobs:
        logger.debug(
            'alignment job: %s -> %s (%s)', job.features_fasta, job.target_fasta, job.output_sam
        )
    if processes > 1 and len(jobs) > 1 and aligner.parallel_safe and can_use_multiprocessing():
        from multiprocessing import Pool

        with Pool(min(processes, len(jobs))) as pool:
            return list(pool.imap(aligner.align, jobs))
    return [aligner.align(job) for job in jobs]


def align_features_to_target(
    ref_chroms: Sequence[str],
    target_chroms: Sequence[str],
    reference: str,
    target: str,
    intermediate_dir: Path,
    minimap2_options: Sequence[str],
    threads: int,
    aligner: Aligner,
    feature_hierarchy: FeatureHierarchy,
    liftover_type: LiftoverType,
    unmapped_features: list[Feature],
) -> AlignedSegments:
    """Align extracted gene sequences to the target and parse the blocks.

    Parameters
    ----------
    ref_chroms : sequence of str
        Reference chromosomes, parallel to ``target_chroms``.
    target_chroms : sequence of str
        Target chromosomes; ``[target]`` means the whole genome.
    reference : str
        Reference genome FASTA.
    target : str
        Target genome FASTA.
    intermediate_dir : pathlib.Path
        Directory for intermediate files.
    minimap2_options : sequence of str
        Arguments passed to minimap2.
    threads : int
        Total threads / worker processes available.
    aligner : Aligner
        Alignment backend.
    feature_hierarchy : FeatureHierarchy
        Reference features.
    liftover_type : LiftoverType
        Pipeline stage.
    unmapped_features : list of Feature
        Receives features that failed to align.

    Returns
    -------
    dict of str to list of AlignedSegment
        Aligned blocks keyed by copy-aware feature name.
    """
    jobs = plan_alignment_jobs(
        ref_chroms,
        target_chroms,
        reference,
        target,
        intermediate_dir,
        minimap2_options,
        threads,
        liftover_type,
    )
    logger.info('aligning features')
    sam_files = run_alignment_jobs(aligner, jobs, threads)
    return parse_all_sam_files(feature_hierarchy, unmapped_features, liftover_type, sam_files)
