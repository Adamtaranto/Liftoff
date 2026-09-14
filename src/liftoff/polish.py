"""Polish lifted genes with broken coding sequences.

For every lifted gene where at least one transcript lost its valid ORF, each
reference exon (with its splice-site dinucleotides) is re-aligned to the
corresponding target region using a semi-global alignment that rewards
matches in coding sequence and splice sites more than elsewhere. The resulting
alignments are written as SAM records and lifted again; the polished gene
replaces the original lift-over when it has more valid ORFs or better identity.
"""

from __future__ import annotations

from collections.abc import Sequence
from itertools import groupby
from pathlib import Path
from typing import Any, TextIO

from Bio.Seq import reverse_complement

from liftoff.align.semiglobal import Alignment, ScoringMatrix, sg_dx_trace
from liftoff.cds import features_by_type, features_lookup
from liftoff.models import Feature, FeatureHierarchy
from liftoff.utils import LiftedFeatures, merge_children_intervals, merge_intervals

#: Gap penalties used when re-aligning exons.
GAP_OPEN = 10
GAP_EXTEND = 1

#: Name of the SAM file written for each polished gene.
POLISH_SAM_NAME = 'polish.sam'


def make_scoring_matrix(match_reward_coding: int) -> ScoringMatrix:
    """Build the case-aware nucleotide scoring matrix.

    Reference bases in CDS and splice sites are upper case, all other bases
    (including the whole target) are lower case. Matching an upper-case
    reference base to its lower-case target base scores
    ``match_reward_coding``; other matches score 2 and mismatches -4.

    Parameters
    ----------
    match_reward_coding : int
        Reward for matching a coding/splice-site base.

    Returns
    -------
    ScoringMatrix
        Matrix over the alphabet ``ACGTacgt*``.
    """
    matrix = ScoringMatrix.create('ACGTacgt*', 2, -4, case_sensitive=True)
    upper_lower_pairs = [(0, 4), (1, 5), (2, 6), (3, 7)]
    updates = [(a, b, match_reward_coding) for a, b in upper_lower_pairs]
    updates += [(b, a, match_reward_coding) for a, b in upper_lower_pairs]
    return matrix.with_scores(updates)


def polish_annotations(
    feature_list: LiftedFeatures,
    ref_fasta: Any,
    target_fasta: Any,
    intermediate_dir: Path,
    feature_hierarchy: FeatureHierarchy,
    target_feature: str,
    matrix: ScoringMatrix | None = None,
) -> bool:
    """Re-align the exons of one lifted gene if any transcript is broken.

    Parameters
    ----------
    feature_list : dict of str to list of Feature
        Lifted features with ORF annotations from :func:`liftoff.cds.check_cds`.
    ref_fasta : pyfaidx.Fasta
        Reference genome.
    target_fasta : pyfaidx.Fasta
        Target genome.
    intermediate_dir : pathlib.Path
        Directory receiving ``polish.sam``.
    feature_hierarchy : FeatureHierarchy
        Reference features (not modified).
    target_feature : str
        Copy ID of the gene to polish.
    matrix : ScoringMatrix or None, default None
        Scoring matrix; built with :func:`make_scoring_matrix` when omitted.

    Returns
    -------
    bool
        ``True`` if a SAM file was written and the gene should be re-lifted.
    """
    target_sub_features = feature_list[target_feature]
    gene_id = target_sub_features[0].id
    ref_sub_features = feature_hierarchy.children[gene_id]
    ref_gene = feature_hierarchy.parents[gene_id]
    transcripts = [f for f in target_sub_features if 'matches_ref_protein' in f.attributes]
    good_transcripts = [t for t in transcripts if t.attributes['valid_ORF'] == ['True']]
    if len(transcripts) == len(good_transcripts):
        return False
    with (intermediate_dir / POLISH_SAM_NAME).open('w') as output_sam:
        write_sam_header(target_fasta, output_sam)
        polish_annotation(
            ref_gene,
            target_sub_features[0],
            ref_sub_features,
            target_sub_features,
            ref_fasta,
            target_fasta,
            output_sam,
            matrix or make_scoring_matrix(3),
        )
    return True


def write_sam_header(target_fasta: Any, output_sam: TextIO) -> None:
    """Write ``@SQ`` header lines for every target sequence.

    Parameters
    ----------
    target_fasta : pyfaidx.Fasta
        Target genome.
    output_sam : TextIO
        Destination stream.
    """
    for name in target_fasta.keys():  # noqa: SIM118 - iterating a Fasta yields records
        output_sam.write(f'@SQ\tSN:{name}\tLN:{len(target_fasta[name])}\n')


def polish_annotation(
    ref_gene: Feature,
    target_gene: Feature,
    ref_children: Sequence[Feature],
    target_children: Sequence[Feature],
    ref_fasta: Any,
    target_fasta: Any,
    output_sam: TextIO,
    matrix: ScoringMatrix,
) -> None:
    """Align each merged reference exon to its target region and write SAM records.

    Parameters
    ----------
    ref_gene : Feature
        Reference top-level feature.
    target_gene : Feature
        Lifted top-level feature.
    ref_children : sequence of Feature
        Reference child features.
    target_children : sequence of Feature
        Lifted features of the gene.
    ref_fasta : pyfaidx.Fasta
        Reference genome.
    target_fasta : pyfaidx.Fasta
        Target genome.
    output_sam : TextIO
        Destination for SAM records.
    matrix : ScoringMatrix
        Scoring matrix for the alignments.
    """
    ref_exons = features_by_type(ref_children, 'exon')
    target_exons = features_by_type(target_children, 'exon')
    ref_cds = features_by_type(ref_children, 'CDS')
    if not ref_exons:
        # CDS-only annotations: treat CDS features as exons.
        ref_exons = ref_cds
        target_exons = features_by_type(target_children, 'CDS')
    ref_cds_intervals = merge_children_intervals(ref_cds)
    # Reference exons are not modified; extended copies of their intervals are
    # used for alignment.
    exon_spans, splice_sites = add_splice_sites(ref_exons, ref_gene)
    merged_ref_intervals = merge_intervals((start, end) for _, start, end in exon_spans)
    exon_groups = find_overlapping_exon_groups(merged_ref_intervals, exon_spans)
    target_exons_by_id = features_lookup(target_exons)
    is_reverse = ref_gene.strand != target_gene.strand
    for ref_interval, exon_group in zip(merged_ref_intervals, exon_groups, strict=True):
        target_interval = get_target_interval(
            exon_group, target_exons, target_exons_by_id, ref_interval
        )
        # Allow extra target sequence when the lifted exon is shorter than the
        # reference, plus a fixed margin.
        flank = (
            max(0, (ref_interval[1] - ref_interval[0]) - (target_interval[1] - target_interval[0]))
            + 10
        )
        if target_interval == [0, 0]:
            continue
        ref_seq = get_feature_sequence(ref_interval, ref_fasta, ref_gene.seqid, 0)
        capitalized_ref_seq = cds_and_splice_sites_to_upper(
            ref_interval, splice_sites, ref_cds_intervals, ref_seq, is_reverse
        )
        target_seq = get_feature_sequence(target_interval, target_fasta, target_gene.seqid, flank)
        alignment = sg_dx_trace(capitalized_ref_seq, target_seq, GAP_OPEN, GAP_EXTEND, matrix)
        write_sam_record(
            ref_interval,
            target_interval,
            ref_gene,
            target_gene,
            capitalized_ref_seq,
            alignment,
            output_sam,
            flank,
        )


#: An exon ID with its splice-site-extended ``start`` and ``end``.
ExonSpan = tuple[str, int, int]


def add_splice_sites(
    exons: Sequence[Feature], parent: Feature
) -> tuple[list[ExonSpan], list[list[int]]]:
    """Extend exon intervals by their two-base splice sites.

    The exon features themselves are not modified.

    Parameters
    ----------
    exons : sequence of Feature
        Reference exons.
    parent : Feature
        Reference gene bounding the extension; boundaries that would extend
        beyond the gene are left unchanged.

    Returns
    -------
    exon_spans : list of tuple
        ``(exon_id, start, end)`` for each exon including its splice sites.
    splice_sites : list of list of int
        ``[start, end]`` of each splice-site dinucleotide added.
    """
    exon_spans: list[ExonSpan] = []
    splice_sites: list[list[int]] = []
    for exon in exons:
        start, end = exon.start, exon.end
        if start - 2 >= parent.start:
            start -= 2
            splice_sites.append([start, start + 1])
        if end + 2 <= parent.end:
            end += 2
            splice_sites.append([end - 1, end])
        exon_spans.append((exon.id, start, end))
    return exon_spans, splice_sites


def find_overlapping_exon_groups(
    merged_ref_intervals: Sequence[Sequence[int]], exon_spans: Sequence[ExonSpan]
) -> list[list[str]]:
    """List the exon IDs starting inside each merged interval.

    Parameters
    ----------
    merged_ref_intervals : sequence of sequence of int
        Merged exon intervals.
    exon_spans : sequence of tuple
        ``(exon_id, start, end)`` extended exon intervals.

    Returns
    -------
    list of list of str
        Exon IDs for each interval, in exon order.
    """
    return [
        [exon_id for exon_id, exon_start, _ in exon_spans if start <= exon_start <= end]
        for start, end in merged_ref_intervals
    ]


def get_target_interval(
    exon_group: Sequence[str],
    target_exons: Sequence[Feature],
    target_exons_by_id: dict[str, Feature] | Any,
    ref_interval: Sequence[int],
) -> list[int]:
    """Return the target region to search for a reference exon group.

    Parameters
    ----------
    exon_group : sequence of str
        Reference exon IDs in the interval.
    target_exons : sequence of Feature
        Lifted exons of the gene.
    target_exons_by_id : mapping of str to Feature
        ``target_exons`` keyed by ID (first occurrence).
    ref_interval : sequence of int
        Reference ``[start, end]`` interval.

    Returns
    -------
    list of int
        Target ``[start, end]``. When none of the exons were lifted, the span
        of all lifted exons widened by the reference interval length.
    """
    target_group = [
        target_exons_by_id[exon_id] for exon_id in exon_group if exon_id in target_exons_by_id
    ]
    if target_group:
        return [min(exon.start for exon in target_group), max(exon.end for exon in target_group)]
    ref_length = ref_interval[1] - ref_interval[0]
    return [
        min(exon.start for exon in target_exons) - ref_length,
        max(exon.end for exon in target_exons) + ref_length,
    ]


def get_feature_sequence(interval: Sequence[int], fasta: Any, chrom: str, flank: float) -> str:
    """Extract a lower-case sequence with optional flanks.

    Parameters
    ----------
    interval : sequence of int
        One-based inclusive ``[start, end]``.
    fasta : pyfaidx.Fasta
        Genome.
    chrom : str
        Sequence name.
    flank : float
        Bases to add on each side (clipped to the sequence).

    Returns
    -------
    str
        Lower-case sequence.
    """
    chrom_seq = fasta[chrom]
    flank_start = get_flank_start(interval[0], flank)
    flank_end = get_flank_end(interval[1], flank, len(chrom_seq))
    return str(chrom_seq[flank_start - 1 : flank_end].seq).lower()


def get_flank_start(start: int, flank: float) -> int:
    """Return the flanked start coordinate, clipped at 1.

    Parameters
    ----------
    start : int
        One-based start.
    flank : float
        Flank length.

    Returns
    -------
    int
        ``max(1, start - flank)`` rounded to an integer.
    """
    return round(max(1, start - flank))


def get_flank_end(end: int, flank: float, chrom_length: int) -> int:
    """Return the flanked end coordinate, clipped at the sequence length.

    Parameters
    ----------
    end : int
        One-based end.
    flank : float
        Flank length.
    chrom_length : int
        Length of the sequence.

    Returns
    -------
    int
        ``min(end + flank, chrom_length)`` rounded to an integer.
    """
    return round(min(end + flank, chrom_length))


def cds_and_splice_sites_to_upper(
    ref_interval: Sequence[int],
    splice_sites: Sequence[Sequence[int]],
    cds_intervals: Sequence[Sequence[int]],
    ref_seq: str,
    is_reverse: bool,
) -> str:
    """Upper-case CDS and splice-site bases of a reference exon sequence.

    Parameters
    ----------
    ref_interval : sequence of int
        Genomic ``[start, end]`` of ``ref_seq``.
    splice_sites : sequence of sequence of int
        Splice-site intervals.
    cds_intervals : sequence of sequence of int
        Merged CDS intervals.
    ref_seq : str
        Lower-case reference sequence.
    is_reverse : bool
        Reverse-complement the result (gene changed strand in the target).

    Returns
    -------
    str
        Mixed-case sequence.
    """
    bases = list(ref_seq)
    for start, end in [*splice_sites, *cds_intervals]:
        # Only segments entirely within this exon interval are capitalised.
        if start >= ref_interval[0] and end <= ref_interval[1]:
            rel_start, rel_end = start - ref_interval[0], end - ref_interval[0]
            bases[rel_start : rel_end + 1] = [b.upper() for b in bases[rel_start : rel_end + 1]]
    result = ''.join(bases)
    if is_reverse:
        return str(reverse_complement(result))
    return result


def write_sam_record(
    ref_interval: Sequence[int],
    target_interval: Sequence[int],
    ref_gene: Feature,
    target_gene: Feature,
    ref_seq: str,
    alignment: Alignment,
    output_sam: TextIO,
    flank: float,
) -> None:
    """Write one exon alignment as a SAM record in gene coordinates.

    Parameters
    ----------
    ref_interval : sequence of int
        Reference exon interval.
    target_interval : sequence of int
        Target search interval (before flanking).
    ref_gene : Feature
        Reference gene; the query offset is expressed as a leading hard clip.
    target_gene : Feature
        Lifted gene.
    ref_seq : str
        Aligned reference sequence.
    alignment : Alignment
        Exon-to-target alignment.
    output_sam : TextIO
        Destination stream.
    flank : float
        Flank added to the target interval.
    """
    is_reverse = ref_gene.strand != target_gene.strand
    if is_reverse:
        hard_clip_start = ref_gene.end - ref_interval[1]
    else:
        hard_clip_start = ref_interval[0] - ref_gene.start
    cigar_string, target_start_offset = make_cigar(
        alignment.query_traceback, alignment.ref_traceback, hard_clip_start
    )
    fields = (
        target_gene.id,
        16 if is_reverse else 0,
        target_gene.seqid,
        get_flank_start(target_interval[0], flank) + target_start_offset,
        0,
        cigar_string,
        '*',
        0,
        0,
        ref_seq.upper(),
        '*',
    )
    output_sam.write('\t'.join(map(str, fields)) + '\n')


def make_cigar(
    reference_traceback: str, target_traceback: str, hard_clip_start: int
) -> tuple[str, int]:
    """Convert gapped alignment strings into an extended CIGAR.

    Parameters
    ----------
    reference_traceback : str
        Gapped reference (query) sequence.
    target_traceback : str
        Gapped target sequence.
    hard_clip_start : int
        Reference bases preceding this exon in the gene, emitted as ``H``.

    Returns
    -------
    cigar : str
        CIGAR using ``=``, ``X``, ``I``, ``D`` and ``H``.
    target_start_offset : int
        Number of leading alignment columns before the first reference base.
    """
    start = next(pos for pos, char in enumerate(reference_traceback) if char != '-')
    ops = []
    for ref_char, target_char in zip(
        reference_traceback[start:], target_traceback[start:], strict=True
    ):
        if ref_char.upper() == target_char.upper():
            ops.append('=')
        elif ref_char == '-':
            ops.append('D')
        elif target_char == '-':
            ops.append('I')
        else:
            ops.append('X')
    return condense_cigar_string(''.join(ops), hard_clip_start), start


def condense_cigar_string(expanded_cigar: str, hard_clip_start: int) -> str:
    """Run-length encode a per-column CIGAR.

    Parameters
    ----------
    expanded_cigar : str
        One operation character per alignment column.
    hard_clip_start : int
        Leading hard clip length (omitted when zero).

    Returns
    -------
    str
        Compact CIGAR string.

    Examples
    --------
    >>> condense_cigar_string('===XX=', 3)
    '3H3=2X1='
    """
    runs = [(sum(1 for _ in group), op) for op, group in groupby(expanded_cigar)]
    if runs and runs[0][1] == 'D':
        runs = runs[1:]
    cigar = ''.join(f'{length}{op}' for length, op in runs)
    if hard_clip_start > 0:
        cigar = f'{hard_clip_start}H{cigar}'
    return cigar
