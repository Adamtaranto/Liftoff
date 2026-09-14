"""Tests for polishing helpers and polished-gene selection."""

from __future__ import annotations

import pytest

from liftoff.models import Attributes, Feature
from liftoff.polish import (
    add_splice_sites,
    cds_and_splice_sites_to_upper,
    condense_cigar_string,
    find_overlapping_exon_groups,
    make_cigar,
    polished_is_better,
)

from ..helpers import FeatureFactory


class TestAddSpliceSites:
    @pytest.fixture
    def gene(self, make_feature: FeatureFactory) -> Feature:
        return make_feature('g', 'gene', 100, 200)

    @pytest.fixture
    def exons(self, make_feature: FeatureFactory) -> list[Feature]:
        return [
            make_feature('e1', 'exon', 101, 120),  # 5' site would extend beyond the gene
            make_feature('e2', 'exon', 150, 170),
            make_feature('e3', 'exon', 180, 199),  # 3' site would extend beyond the gene
        ]

    def test_extends_exons_within_the_gene(self, gene: Feature, exons: list[Feature]) -> None:
        spans, _ = add_splice_sites(exons, gene)
        assert spans == [('e1', 101, 122), ('e2', 148, 172), ('e3', 178, 199)]

    def test_reports_added_splice_sites(self, gene: Feature, exons: list[Feature]) -> None:
        _, sites = add_splice_sites(exons, gene)
        assert sites == [[121, 122], [148, 149], [171, 172], [178, 179]]

    def test_does_not_modify_exons(self, gene: Feature, exons: list[Feature]) -> None:
        before = [(exon.start, exon.end) for exon in exons]
        add_splice_sites(exons, gene)
        assert [(exon.start, exon.end) for exon in exons] == before

    def test_is_repeatable(self, gene: Feature, exons: list[Feature]) -> None:
        assert add_splice_sites(exons, gene) == add_splice_sites(exons, gene)


def test_exon_groups_are_assigned_by_extended_start() -> None:
    spans = [('e1', 101, 122), ('e2', 120, 140), ('e3', 178, 199)]
    assert find_overlapping_exon_groups([[101, 140], [178, 199]], spans) == [['e1', 'e2'], ['e3']]


class TestCdsAndSpliceSitesToUpper:
    def test_capitalises_segments_inside_interval(self) -> None:
        result = cds_and_splice_sites_to_upper([10, 19], [[10, 11]], [[14, 16]], 'a' * 10, False)
        assert result == 'AAaaAAAaaa'

    def test_ignores_segments_crossing_interval_edges(self) -> None:
        assert cds_and_splice_sites_to_upper([10, 19], [[8, 11]], [], 'a' * 10, False) == 'a' * 10

    def test_reverse_complements_when_strand_changes(self) -> None:
        assert cds_and_splice_sites_to_upper([1, 4], [], [[1, 2]], 'acgt', True) == 'acGT'


class TestCigar:
    @pytest.mark.parametrize(
        ('ref', 'target', 'cigar', 'offset'),
        [
            pytest.param('ACGT', 'acgt', '4=', 0, id='all-match'),
            pytest.param('--ACGT', 'ttacgt', '4=', 2, id='leading-target-gap-offsets'),
            pytest.param('ACGT', 'aggt', '1=1X2=', 0, id='mismatch'),
            pytest.param('AC-GT', 'acggt', '2=1D2=', 0, id='deletion-from-reference'),
            pytest.param('ACGT', 'a-gt', '1=1I2=', 0, id='insertion-in-reference'),
        ],
    )
    def test_make_cigar(self, ref: str, target: str, cigar: str, offset: int) -> None:
        assert make_cigar(ref, target, 0) == (cigar, offset)

    @pytest.mark.parametrize(
        ('expanded', 'hard_clip', 'expected'),
        [
            pytest.param('===XX=', 3, '3H3=2X1=', id='with-hard-clip'),
            pytest.param('==', 0, '2=', id='no-hard-clip'),
            pytest.param('DD==', 0, '2=', id='leading-deletion-dropped'),
        ],
    )
    def test_condense_cigar_string(self, expanded: str, hard_clip: int, expected: str) -> None:
        assert condense_cigar_string(expanded, hard_clip) == expected


def attributes(orfs: str | None, identity: str, coverage: str) -> Attributes:
    result = {'sequence_ID': [identity], 'coverage': [coverage]}
    if orfs is not None:
        result['valid_ORFs'] = [orfs]
    return result


@pytest.mark.parametrize(
    ('original', 'polished', 'expected'),
    [
        pytest.param(
            attributes('0', '1.0', '1.0'), attributes('1', '0.9', '0.9'), True, id='more-valid-orfs'
        ),
        pytest.param(
            attributes('10', '0.9', '1.0'),
            attributes('9', '1.0', '1.0'),
            False,
            id='orf-count-numeric',
        ),
        pytest.param(
            attributes('1', '0.99', '1.0'),
            attributes('1', '1.0', '1.0'),
            True,
            id='identity-numeric',
        ),
        pytest.param(
            attributes('1', '1.0', '0.9'),
            attributes('1', '0.99', '1.0'),
            False,
            id='identity-before-coverage',
        ),
        pytest.param(
            attributes('1', '0.99', '0.9'),
            attributes('1', '0.99', '1.0'),
            True,
            id='coverage-breaks-tie',
        ),
        pytest.param(
            attributes('1', '0.99', '1.0'),
            attributes('1', '0.99', '1.0'),
            False,
            id='equal-keeps-original',
        ),
        pytest.param(
            attributes('1', '0.9', '0.9'),
            attributes(None, '1.0', '1.0'),
            False,
            id='cds-lost-when-polishing',
        ),
    ],
)
def test_polished_is_better(original: Attributes, polished: Attributes, expected: bool) -> None:
    assert polished_is_better(original, polished) is expected
