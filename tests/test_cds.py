"""Tests for ORF detection and trimming of broken CDS features."""

from __future__ import annotations

from Bio.Seq import reverse_complement
import pytest

from liftoff.cds import find_longest_orf, trim_cds_to_orf
from liftoff.models import Feature


def _cds(fid: str, start: int, end: int, strand: str = '+') -> Feature:
    return Feature(fid, 'CDS', 'chr1', 'Liftoff', strand, start, end, '0', {'Parent': ['t1']})


def _spliced(genome: str, cds_group: list[Feature]) -> str:
    seq = ''.join(genome[c.start - 1 : c.end] for c in sorted(cds_group, key=lambda c: c.start))
    return str(reverse_complement(seq)) if cds_group[0].strand == '-' else seq


@pytest.mark.parametrize(
    ('seq', 'min_length', 'expected'),
    [
        ('CCATGAAATAGCC', 6, (2, 11)),
        # The whole sequence is already an ORF: nothing to trim.
        ('ATGAAATAG', 6, None),
        # Too short.
        ('CCATGAAATAGCC', 12, None),
        # No stop codon after the ATG.
        ('CCATGAAAAAA', 3, None),
        # Longest ORF is in another frame than an earlier, shorter one.
        ('ATGTAACATGAAAAAAAAATAAG', 6, (7, 22)),
        # A nested ATG before the same stop does not replace the outer ORF.
        ('CATGATGAAATAA', 6, (1, 13)),
        # Equal lengths: earliest start wins.
        ('ATGTAAATGTAA', 3, (0, 6)),
    ],
)
def test_find_longest_orf(seq: str, min_length: int, expected: tuple[int, int] | None) -> None:
    assert find_longest_orf(seq, min_length=min_length) == expected


def test_find_longest_orf_matches_brute_force() -> None:
    import random

    rng = random.Random(3)
    stops = ('TAA', 'TAG', 'TGA')
    for _ in range(300):
        seq = ''.join(rng.choice('ACGT') for _ in range(rng.randint(3, 120)))
        best = None
        for i in range(len(seq) - 2):
            if seq[i : i + 3] != 'ATG':
                continue
            for j in range(i, len(seq) - 2, 3):
                if seq[j : j + 3] in stops:
                    orf = (i, j + 3)
                    if best is None or orf[1] - orf[0] > best[1] - best[0]:
                        best = orf
                    break
        if best is not None and best[1] - best[0] == len(seq):
            best = None
        assert find_longest_orf(seq, min_length=0) == best, seq


def test_trim_plus_strand_multi_exon() -> None:
    # Transcript offsets: exon1 0-9, exon2 10-19, exon3 20-29.
    cds = [_cds('c1', 101, 110), _cds('c2', 201, 210), _cds('c3', 301, 310)]
    gene = [Feature('g', 'gene', 'chr1', 'L', '+', 101, 310, '.', {}), *cds]
    trim_cds_to_orf(cds, (12, 18), gene)  # ORF covers transcript bases 12..17
    assert [(c.id, c.start, c.end, c.frame) for c in cds] == [('c2', 203, 208, '0')]
    assert [f.id for f in gene] == ['g', 'c2']


def test_trim_minus_strand_spanning_exons_recomputes_phase() -> None:
    # Minus strand: transcript starts at the highest coordinate (c3 end).
    cds = [_cds('c1', 101, 110, '-'), _cds('c2', 201, 210, '-'), _cds('c3', 301, 310, '-')]
    gene = list(cds)
    trim_cds_to_orf(cds, (4, 26), gene)  # transcript bases 4..25
    by_id = {c.id: (c.start, c.end, c.frame) for c in cds}
    # c3 (offsets 0-9) loses 4 bases at its 5' (high) end; c1 (20-29) keeps offsets 20-25.
    assert by_id == {'c3': (301, 306, '0'), 'c2': (201, 210, '0'), 'c1': (105, 110, '2')}


def test_trimmed_cds_sequence_is_the_orf() -> None:
    import random

    rng = random.Random(11)
    genome = ''.join(rng.choice('ACGT') for _ in range(400))
    for strand in '+-':
        cds = [
            _cds('c1', 11, 60, strand),
            _cds('c2', 101, 180, strand),
            _cds('c3', 251, 330, strand),
        ]
        before = _spliced(genome, cds)
        orf = (37, 151)
        trim_cds_to_orf(cds, orf, list(cds))
        assert _spliced(genome, cds) == before[orf[0] : orf[1]]
