"""Tests for ORF detection, CDS trimming and CDS status annotation."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
import random
from typing import Any

from Bio.Seq import reverse_complement
import pytest

from liftoff.cds import (
    find_and_check_cds,
    find_longest_orf,
    inframe_stop,
    missing_start,
    missing_stop,
    trim_cds_to_orf,
)
from liftoff.io.fasta import open_fasta
from liftoff.models import Feature

from ..helpers import FeatureFactory

STOP_CODONS = ('TAA', 'TAG', 'TGA')


def brute_force_longest_orf(seq: str) -> tuple[int, int] | None:
    """Reference implementation: longest ATG..stop ORF, earliest start on ties."""
    best = None
    for start in range(len(seq) - 2):
        if seq[start : start + 3] != 'ATG':
            continue
        for codon in range(start, len(seq) - 2, 3):
            if seq[codon : codon + 3] in STOP_CODONS:
                if best is None or codon + 3 - start > best[1] - best[0]:
                    best = (start, codon + 3)
                break
    if best is not None and best[1] - best[0] == len(seq):
        return None
    return best


def spliced(genome: str, cds_group: list[Feature]) -> str:
    """Concatenate CDS sequence in transcript orientation."""
    seq = ''.join(genome[c.start - 1 : c.end] for c in sorted(cds_group, key=lambda c: c.start))
    return str(reverse_complement(seq)) if cds_group[0].strand == '-' else seq


class TestFindLongestOrf:
    @pytest.mark.parametrize(
        ('seq', 'min_length', 'expected'),
        [
            pytest.param('CCATGAAATAGCC', 6, (2, 11), id='internal-orf'),
            pytest.param('ATGAAATAG', 6, None, id='whole-sequence-is-orf'),
            pytest.param('CCATGAAATAGCC', 12, None, id='shorter-than-minimum'),
            pytest.param('CCATGAAAAAA', 3, None, id='no-stop-codon'),
            pytest.param('ATGTAACATGAAAAAAAAATAAG', 6, (7, 22), id='longest-in-other-frame'),
            pytest.param('CATGATGAAATAA', 6, (1, 13), id='nested-atg-ignored'),
            pytest.param('ATGTAAATGTAA', 3, (0, 6), id='tie-prefers-earliest'),
        ],
    )
    def test_examples(self, seq: str, min_length: int, expected: tuple[int, int] | None) -> None:
        assert find_longest_orf(seq, min_length=min_length) == expected

    def test_default_minimum_is_180_bases(self) -> None:
        orf = 'ATG' + 'AAA' * 58 + 'TAA'  # 180 bases
        assert (find_longest_orf('CC' + orf), find_longest_orf('CC' + orf[3:])) == ((2, 182), None)

    @pytest.mark.parametrize('seed', range(20))
    def test_agrees_with_brute_force(self, seed: int) -> None:
        generator = random.Random(seed)
        for _ in range(20):
            seq = ''.join(generator.choice('ACGT') for _ in range(generator.randint(3, 150)))
            assert find_longest_orf(seq, min_length=0) == brute_force_longest_orf(seq), seq


class TestTrimCdsToOrf:
    @pytest.fixture
    def plus_cds(self, make_feature: FeatureFactory) -> list[Feature]:
        """Three 10-base CDS pieces: transcript offsets 0-9, 10-19 and 20-29."""
        return [
            make_feature('c1', 'CDS', 101, 110, frame='0'),
            make_feature('c2', 'CDS', 201, 210, frame='0'),
            make_feature('c3', 'CDS', 301, 310, frame='0'),
        ]

    @pytest.fixture
    def minus_cds(self, make_feature: FeatureFactory) -> list[Feature]:
        """As ``plus_cds`` on the minus strand; transcript offset 0 is at position 310."""
        return [
            make_feature('c1', 'CDS', 101, 110, strand='-', frame='0'),
            make_feature('c2', 'CDS', 201, 210, strand='-', frame='0'),
            make_feature('c3', 'CDS', 301, 310, strand='-', frame='0'),
        ]

    def test_trims_plus_strand_boundaries(self, plus_cds: list[Feature]) -> None:
        trim_cds_to_orf(plus_cds, (12, 18), list(plus_cds))
        assert [(c.id, c.start, c.end) for c in plus_cds] == [('c2', 203, 208)]

    def test_removes_cds_outside_orf_from_gene(
        self, plus_cds: list[Feature], make_feature: FeatureFactory
    ) -> None:
        gene = make_feature('g', 'gene', 101, 310)
        gene_features = [gene, *plus_cds]
        trim_cds_to_orf(plus_cds, (12, 18), gene_features)
        assert [f.id for f in gene_features] == ['g', 'c2']

    def test_trims_minus_strand_boundaries(self, minus_cds: list[Feature]) -> None:
        trim_cds_to_orf(minus_cds, (4, 26), list(minus_cds))
        assert {c.id: (c.start, c.end) for c in minus_cds} == {
            'c1': (105, 110),
            'c2': (201, 210),
            'c3': (301, 306),
        }

    def test_recomputes_phases_in_transcript_order(self, minus_cds: list[Feature]) -> None:
        trim_cds_to_orf(minus_cds, (4, 26), list(minus_cds))
        # c3 (6 bases) then c2 (10 bases): c1 starts 16 bases in, so phase 2.
        assert {c.id: c.frame for c in minus_cds} == {'c3': '0', 'c2': '0', 'c1': '2'}

    @pytest.mark.parametrize('strand', ['+', '-'])
    def test_trimmed_sequence_equals_orf(
        self, strand: str, rng: random.Random, make_feature: FeatureFactory
    ) -> None:
        genome = ''.join(rng.choice('ACGT') for _ in range(400))
        cds_group = [
            make_feature('c1', 'CDS', 11, 60, strand=strand),
            make_feature('c2', 'CDS', 101, 180, strand=strand),
            make_feature('c3', 'CDS', 251, 330, strand=strand),
        ]
        before = spliced(genome, cds_group)

        trim_cds_to_orf(cds_group, (37, 151), list(cds_group))

        assert spliced(genome, cds_group) == before[37:151]


@pytest.mark.parametrize(
    ('check', 'protein', 'expected'),
    [
        pytest.param(missing_start, 'MKL*', False, id='start-present'),
        pytest.param(missing_start, 'KL*', True, id='start-missing'),
        pytest.param(missing_stop, 'MKL*', False, id='stop-present'),
        pytest.param(missing_stop, 'MKL', True, id='stop-missing'),
        pytest.param(inframe_stop, 'MK*L*', True, id='premature-stop'),
        pytest.param(inframe_stop, 'MKL*', False, id='no-premature-stop'),
    ],
)
def test_protein_status_checks(check: Callable[[str], bool], protein: str, expected: bool) -> None:
    assert check(protein) is expected


class TestFindAndCheckCds:
    """Status annotation on a tiny genome: ATG AAA AAA TAA split over two CDS pieces."""

    @pytest.fixture
    def fasta(self, tmp_path: Path) -> Any:
        path = tmp_path / 'genome.fa'
        path.write_text('>chr1\nCCATGAAAGGGGGAAATAACC\n')
        return open_fasta(str(path))

    @pytest.fixture
    def gene_features(self, make_feature: FeatureFactory) -> list[Feature]:
        return [
            make_feature('g', 'gene', 3, 19),
            make_feature('t', 'mRNA', 3, 19, parent='g'),
            make_feature('c1', 'CDS', 3, 8, parent='t'),
            make_feature('c2', 'CDS', 14, 19, parent='t'),
        ]

    def test_counts_valid_orfs_on_gene(self, fasta: Any, gene_features: list[Feature]) -> None:
        find_and_check_cds(gene_features, gene_features[1:], fasta, fasta)
        assert gene_features[0].attributes['valid_ORFs'] == ['1']

    def test_flags_transcript_as_valid_and_matching(
        self, fasta: Any, gene_features: list[Feature]
    ) -> None:
        find_and_check_cds(gene_features, gene_features[1:], fasta, fasta)
        transcript = gene_features[1].attributes
        assert (transcript['valid_ORF'], transcript['matches_ref_protein']) == (['True'], ['True'])

    def test_gene_without_cds_is_not_annotated(
        self, fasta: Any, make_feature: FeatureFactory
    ) -> None:
        features = [make_feature('g', 'gene', 3, 19)]
        assert find_and_check_cds(features, [], fasta, fasta) == (0, 0)
        assert 'valid_ORFs' not in features[0].attributes
