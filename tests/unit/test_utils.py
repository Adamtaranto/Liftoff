"""Tests for coordinate, identifier, ordering and model helpers."""

from __future__ import annotations

import pytest

from liftoff.models import AlignedSegment, Feature
from liftoff.utils import (
    ParentOrder,
    clear_scores,
    convert_id_to_original,
    count_overlap,
    find_nonoverlapping_upstream_neighbor,
    get_copy_tag,
    get_parent_list,
    get_relative_child_coord,
    get_strand,
    merge_children_intervals,
    merge_intervals,
    overlaps_in_ref_annotation,
)

from ..helpers import FeatureFactory


@pytest.mark.parametrize(
    ('feature_id', 'original', 'tag'),
    [
        pytest.param('gene1_0', 'gene1', '_0', id='copy-zero'),
        pytest.param(
            'gene_with_underscores_12', 'gene_with_underscores', '_12', id='underscored-id'
        ),
        pytest.param('gene1_3_frag2', 'gene1', '_frag2', id='fragment'),
    ],
)
class TestCopyIds:
    def test_original_id(self, feature_id: str, original: str, tag: str) -> None:
        assert convert_id_to_original(feature_id) == original

    def test_copy_tag(self, feature_id: str, original: str, tag: str) -> None:
        assert get_copy_tag(feature_id) == tag


@pytest.mark.parametrize(
    ('interval_a', 'interval_b', 'overlap'),
    [
        pytest.param((1, 10), (5, 20), 6, id='partial'),
        pytest.param((1, 10), (10, 20), 1, id='single-base'),
        pytest.param((1, 10), (11, 20), 0, id='adjacent'),
        pytest.param((1, 100), (20, 30), 11, id='contained'),
    ],
)
def test_count_overlap(
    interval_a: tuple[int, int], interval_b: tuple[int, int], overlap: int
) -> None:
    assert count_overlap(*interval_a, *interval_b) == overlap


@pytest.mark.parametrize(('is_reverse', 'expected'), [(False, 10), (True, 90)])
def test_relative_child_coordinate(
    make_feature: FeatureFactory, is_reverse: bool, expected: int
) -> None:
    assert get_relative_child_coord(make_feature(start=100, end=200), 110, is_reverse) == expected


class TestMergeIntervals:
    def test_merges_overlapping_and_sorts(self) -> None:
        assert merge_intervals([(50, 60), (1, 10), (5, 20)]) == [[1, 20], [50, 60]]

    def test_keeps_adjacent_intervals_separate(self) -> None:
        assert merge_intervals([(1, 10), (11, 20)]) == [[1, 10], [11, 20]]

    def test_empty(self) -> None:
        assert merge_intervals([]) == []

    def test_children_use_annotated_bounds(self, make_feature: FeatureFactory) -> None:
        flanked = make_feature(start=90, end=210)
        flanked.original_bounds = (100, 200)
        assert merge_children_intervals([flanked]) == [[100, 200]]


@pytest.mark.parametrize(
    ('parent_strand', 'is_reverse', 'expected'),
    [
        pytest.param('+', False, '+', id='forward-plus'),
        pytest.param('-', False, '-', id='forward-minus'),
        pytest.param('+', True, '-', id='reverse-plus'),
        pytest.param('-', True, '+', id='reverse-minus'),
    ],
)
def test_get_strand(
    make_feature: FeatureFactory, parent_strand: str, is_reverse: bool, expected: str
) -> None:
    block = AlignedSegment(1, 'gene1_0', 'chr1', 0, 9, 0, 9, is_reverse)
    assert get_strand(block, make_feature(strand=parent_strand)) == expected


def test_parent_list_skips_empty_groups(make_feature: FeatureFactory) -> None:
    gene = make_feature()
    assert get_parent_list({'gene1_0': [gene], 'gene2_0': []}) == [gene]


def test_clear_scores_resets_reference_features_only(make_feature: FeatureFactory) -> None:
    gene, exon = make_feature('gene1'), make_feature('exon1', 'exon')
    gene.score = exon.score = 0.2
    clear_scores({'gene1_0': [gene, exon]}, {'gene1': gene})
    assert (gene.score, exon.score) == (-1, 0.2)


class TestOverlapsInRefAnnotation:
    @pytest.mark.parametrize(
        ('changes', 'expected'),
        [
            pytest.param({}, True, id='overlapping'),
            pytest.param({'seqid': 'chr2'}, False, id='other-sequence'),
            pytest.param({'strand': '-'}, False, id='other-strand'),
            pytest.param({'start': 101, 'end': 200}, False, id='disjoint'),
            pytest.param({'feature_id': 'gene1'}, False, id='same-feature'),
        ],
    )
    def test_overlap(
        self, make_feature: FeatureFactory, changes: dict[str, object], expected: bool
    ) -> None:
        gene = make_feature('gene1', start=1, end=100)
        options: dict[str, object] = {'feature_id': 'gene2', 'start': 50, 'end': 150} | changes
        assert overlaps_in_ref_annotation(gene, make_feature(**options)) is expected


class TestParentOrder:
    @pytest.fixture
    def order(self, make_feature: FeatureFactory) -> ParentOrder:
        return ParentOrder.from_features(
            [
                make_feature('b', start=50),
                make_feature('a', start=1),
                make_feature('c', start=5, seqid='chr2'),
            ]
        )

    def test_sorts_by_sequence_then_start(self, order: ParentOrder) -> None:
        assert order.ids == ['a', 'b', 'c']

    def test_duplicate_ids_resolve_to_first_occurrence(self, make_feature: FeatureFactory) -> None:
        order = ParentOrder.from_features([make_feature('x', start=5), make_feature('x', start=1)])
        assert order.index('x') == 0

    def test_missing_id(self, order: ParentOrder) -> None:
        assert order.index('missing') is None

    @pytest.mark.parametrize(
        ('feature_id', 'neighbour'),
        [
            pytest.param('a', None, id='first-feature'),
            pytest.param('b', 'a', id='same-sequence'),
            pytest.param('c', None, id='first-on-new-sequence'),
        ],
    )
    def test_upstream_neighbour(
        self, order: ParentOrder, feature_id: str, neighbour: str | None
    ) -> None:
        assert find_nonoverlapping_upstream_neighbor(order, feature_id) == neighbour

    def test_upstream_neighbour_of_unknown_feature_raises(self, order: ParentOrder) -> None:
        with pytest.raises(KeyError):
            find_nonoverlapping_upstream_neighbor(order, 'missing')


class TestFeatureModel:
    def test_equality_is_identity(self, make_feature: FeatureFactory) -> None:
        assert make_feature() != make_feature()

    def test_annotated_bounds_default_to_coordinates(self, make_feature: FeatureFactory) -> None:
        assert make_feature(start=5, end=9).annotated_bounds == (5, 9)

    def test_annotated_bounds_prefer_original_bounds(self, make_feature: FeatureFactory) -> None:
        feature: Feature = make_feature(start=1, end=20)
        feature.original_bounds = (5, 9)
        assert feature.annotated_bounds == (5, 9)
