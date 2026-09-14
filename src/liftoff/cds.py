"""Assess whether lifted coding sequences still encode valid ORFs."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
import re
from typing import Any
import warnings

from Bio.Seq import reverse_complement, translate

from liftoff.io.fasta import open_fasta
from liftoff.models import Feature, FeatureHierarchy
from liftoff.utils import LiftedFeatures

#: Minimum length (nt, including the stop codon) of an internal ORF that a
#: broken CDS is trimmed to.
MIN_ORF_LENGTH = 180

_START_CODON = re.compile(r'(?=ATG)')
_STOP_CODON = re.compile(r'(?=TAA|TAG|TGA)')


def check_cds(
    feature_list: LiftedFeatures,
    feature_hierarchy: FeatureHierarchy,
    reference: str,
    target: str,
) -> None:
    """Annotate ORF status on every lifted gene.

    Parameters
    ----------
    feature_list : dict of str to list of Feature
        Lifted features; transcript and gene attributes are updated in place.
    feature_hierarchy : FeatureHierarchy
        Reference features.
    reference : str
        Reference genome FASTA.
    target : str
        Target genome FASTA.
    """
    ref_fasta, target_fasta = open_fasta(reference), open_fasta(target)
    for target_sub_features in feature_list.values():
        ref_sub_features = feature_hierarchy.children[target_sub_features[0].id]
        find_and_check_cds(target_sub_features, ref_sub_features, ref_fasta, target_fasta)


def find_and_check_cds(
    target_sub_features: list[Feature],
    ref_sub_features: Sequence[Feature],
    ref_fasta: Any,
    target_fasta: Any,
) -> tuple[int, int]:
    """Check each transcript's CDS and record the number of valid ORFs.

    Parameters
    ----------
    target_sub_features : list of Feature
        All lifted features of one gene, top-level feature first.
    ref_sub_features : sequence of Feature
        Reference child features of the same gene.
    ref_fasta : pyfaidx.Fasta
        Reference genome.
    target_fasta : pyfaidx.Fasta
        Target genome.

    Returns
    -------
    total_cds : int
        Number of transcripts with CDS features.
    good_cds : int
        Number of those with a valid ORF.
    """
    cds_features = features_by_type(target_sub_features, 'CDS')
    if not cds_features:
        return 0, 0
    total_cds, num_good_cds = _count_good_cds(
        cds_features, ref_fasta, target_fasta, ref_sub_features, target_sub_features
    )
    target_sub_features[0].attributes['valid_ORFs'] = [str(num_good_cds)]
    return total_cds, num_good_cds


def _count_good_cds(
    cds_features: list[Feature],
    ref_fasta: Any,
    target_fasta: Any,
    ref_children: Sequence[Feature],
    target_sub_features: list[Feature],
) -> tuple[int, int]:
    """Translate each transcript's CDS and tag the transcript with its status.

    If a transcript's CDS is not a valid ORF but contains an ORF of at least
    :data:`MIN_ORF_LENGTH` bases, its CDS features are trimmed to that ORF and
    the ORF status is evaluated on the trimmed sequence.
    """
    grouped_cds = group_cds_by_transcript(cds_features)
    features_by_id = features_lookup(target_sub_features)
    good_cds_count = 0
    for cds_group in grouped_cds.values():
        cds_seq = get_seq(cds_group, target_fasta)
        transcript = features_by_id[cds_group[0].attributes['Parent'][0]]
        protein = _translate(cds_seq)
        matches = matches_ref(cds_group, ref_fasta, ref_children, protein)
        transcript.attributes['matches_ref_protein'] = [str(matches)]
        if not _is_valid_orf(protein):
            orf = find_longest_orf(cds_seq)
            if orf is not None:
                trim_cds_to_orf(cds_group, orf, target_sub_features)
                protein = _translate(cds_seq[orf[0] : orf[1]])
        transcript.attributes['valid_ORF'] = ['False']
        if len(protein) < 3:
            transcript.attributes['partial_ORF'] = ['True']
        elif missing_start(protein):
            transcript.attributes['missing_start_codon'] = ['True']
        elif missing_stop(protein):
            transcript.attributes['missing_stop_codon'] = ['True']
        elif inframe_stop(protein):
            transcript.attributes['inframe_stop_codon'] = ['True']
        else:
            transcript.attributes['valid_ORF'] = ['True']
            good_cds_count += 1
    return len(grouped_cds), good_cds_count


def _is_valid_orf(protein: str) -> bool:
    """Return whether a translation is a complete ORF without premature stops."""
    return (
        len(protein) >= 3
        and not missing_start(protein)
        and not missing_stop(protein)
        and not inframe_stop(protein)
    )


def find_longest_orf(cds_seq: str, min_length: int = MIN_ORF_LENGTH) -> tuple[int, int] | None:
    """Locate the longest ATG-initiated ORF within a coding sequence.

    All three reading frames are searched. An ORF runs from an ATG to the
    first in-frame stop codon (inclusive); ties are broken by the earliest
    start.

    Parameters
    ----------
    cds_seq : str
        Upper-case nucleotide sequence in transcript orientation.
    min_length : int, default MIN_ORF_LENGTH
        Minimum ORF length in bases, including the stop codon.

    Returns
    -------
    tuple of int or None
        Half-open ``(start, end)`` offsets of the longest ORF, or ``None`` if
        no ORF reaches ``min_length`` or the ORF spans the whole sequence
        (nothing to trim).

    Examples
    --------
    >>> find_longest_orf('CCATGAAATAGCC', min_length=6)
    (2, 11)
    """
    starts = [match.start() for match in _START_CODON.finditer(cds_seq)]
    stops = [match.start() for match in _STOP_CODON.finditer(cds_seq)]
    best: tuple[int, int] | None = None
    for frame in range(3):
        frame_stops = [pos for pos in stops if pos % 3 == frame]
        stop_index = 0
        covered_until = -1
        for start in (pos for pos in starts if pos % 3 == frame):
            if start < covered_until:
                # A later ATG before the same stop only gives a nested, shorter ORF.
                continue
            while stop_index < len(frame_stops) and frame_stops[stop_index] < start:
                stop_index += 1
            if stop_index == len(frame_stops):
                break
            orf = (start, frame_stops[stop_index] + 3)
            covered_until = orf[1]
            if best is None or (orf[1] - orf[0], -orf[0]) > (best[1] - best[0], -best[0]):
                best = orf
    if best is None or best[1] - best[0] < min_length or best[1] - best[0] == len(cds_seq):
        return None
    return best


def trim_cds_to_orf(
    cds_group: list[Feature], orf: tuple[int, int], gene_features: list[Feature]
) -> None:
    """Restrict a transcript's CDS features to an ORF, in place.

    CDS features entirely outside the ORF are removed from ``gene_features``;
    the CDS features containing the ORF boundaries are shortened, and the
    phase of every remaining CDS is recomputed.

    Parameters
    ----------
    cds_group : list of Feature
        CDS features of one transcript, sorted by genomic start.
    orf : tuple of int
        Half-open ORF offsets in transcript orientation (as returned by
        :func:`find_longest_orf`).
    gene_features : list of Feature
        All lifted features of the gene; removed CDS features are deleted from
        this list.
    """
    orf_first, orf_last = orf[0], orf[1] - 1
    reverse = cds_group[-1].strand == '-'
    # Walk CDS features in transcript order, tracking transcript offsets.
    ordered = list(reversed(cds_group)) if reverse else list(cds_group)
    offset = 0
    kept: list[Feature] = []
    for cds in ordered:
        first, last = offset, offset + cds.end - cds.start
        offset = last + 1
        if last < orf_first or first > orf_last:
            gene_features.remove(cds)
            cds_group.remove(cds)
            continue
        trim_head = max(first, orf_first) - first  # bases removed at the 5' side
        trim_tail = last - min(last, orf_last)  # bases removed at the 3' side
        if reverse:
            cds.end -= trim_head
            cds.start += trim_tail
        else:
            cds.start += trim_head
            cds.end -= trim_tail
        kept.append(cds)
    # GFF3 phase: bases to skip at the 5' end of a CDS to reach a codon start.
    coding_length = 0
    for cds in kept:
        cds.frame = str((3 - coding_length % 3) % 3)
        coding_length += cds.end - cds.start + 1


def group_cds_by_transcript(cds_features: Iterable[Feature]) -> dict[str, list[Feature]]:
    """Group CDS features by their parent transcript.

    Parameters
    ----------
    cds_features : iterable of Feature
        CDS features with ``Parent`` attributes.

    Returns
    -------
    dict of str to list of Feature
        CDS features keyed by transcript ID, in input order.
    """
    groups: dict[str, list[Feature]] = {}
    for cds in cds_features:
        groups.setdefault(cds.attributes['Parent'][0], []).append(cds)
    return groups


def get_seq(cds_group: list[Feature], fasta: Any) -> str:
    """Concatenate a transcript's CDS sequence in transcript orientation.

    Parameters
    ----------
    cds_group : list of Feature
        CDS features of one transcript; sorted by start in place.
    fasta : pyfaidx.Fasta
        Genome containing the features.

    Returns
    -------
    str
        Upper-case coding sequence (reverse-complemented on the minus strand).
    """
    cds_group.sort(key=lambda x: x.start)
    chrom_seq = fasta[cds_group[0].seqid]
    seq = ''.join(str(chrom_seq[cds.start - 1 : cds.end].seq) for cds in cds_group)
    if cds_group[-1].strand == '-':
        seq = str(reverse_complement(seq))
    return seq.upper()


def _translate(seq: str) -> str:
    """Translate, silencing Biopython's partial-codon warnings."""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return str(translate(seq))


def missing_start(protein: str) -> bool:
    """Return whether a protein does not begin with methionine.

    Parameters
    ----------
    protein : str
        Translated sequence.

    Returns
    -------
    bool
        ``True`` if the first residue is not ``M``.
    """
    return protein[0] != 'M'


def missing_stop(protein: str) -> bool:
    """Return whether a protein does not end with a stop codon.

    Parameters
    ----------
    protein : str
        Translated sequence.

    Returns
    -------
    bool
        ``True`` if the last residue is not ``*``.
    """
    return protein[-1] != '*'


def inframe_stop(protein: str) -> bool:
    """Return whether a protein contains a premature stop codon.

    Parameters
    ----------
    protein : str
        Translated sequence.

    Returns
    -------
    bool
        ``True`` if ``*`` occurs before the last residue.
    """
    return '*' in protein[:-1]


def matches_ref(
    cds_group: Sequence[Feature],
    ref_fasta: Any,
    ref_children: Iterable[Feature],
    target_protein: str,
) -> bool:
    """Return whether a lifted protein equals the reference protein.

    Parameters
    ----------
    cds_group : sequence of Feature
        Lifted CDS features of one transcript.
    ref_fasta : pyfaidx.Fasta
        Reference genome.
    ref_children : iterable of Feature
        Reference child features of the gene.
    target_protein : str
        Translation of the lifted CDS.

    Returns
    -------
    bool
        ``True`` if the translations are identical.
    """
    transcript = cds_group[0].attributes['Parent']
    ref_cds_group = [
        ref_cds
        for ref_cds in ref_children
        if ref_cds.featuretype == 'CDS' and ref_cds.attributes['Parent'] == transcript
    ]
    return target_protein == _translate(get_seq(ref_cds_group, ref_fasta))


def features_by_type(features: Iterable[Feature], featuretype: str) -> list[Feature]:
    """Select features of one type.

    Parameters
    ----------
    features : iterable of Feature
        Features to filter.
    featuretype : str
        Type to keep.

    Returns
    -------
    list of Feature
        Matching features in input order.
    """
    return [feature for feature in features if feature.featuretype == featuretype]


def features_lookup(features: Iterable[Feature]) -> Mapping[str, Feature]:
    """Index features by ID, keeping the first occurrence of duplicates.

    Parameters
    ----------
    features : iterable of Feature
        Features to index.

    Returns
    -------
    mapping of str to Feature
        Features keyed by ID.
    """
    lookup: dict[str, Feature] = {}
    for feature in features:
        lookup.setdefault(feature.id, feature)
    return lookup
