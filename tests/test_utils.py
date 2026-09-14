"""Tests for coordinate/identifier helpers and the interval index."""

from __future__ import annotations

import random

import pytest

from liftoff.intervals import Interval, IntervalIndex
from liftoff.models import AlignedSegment, Feature
from liftoff.utils import (
    ParentOrder,
    convert_id_to_original,
    count_overlap,
    find_nonoverlapping_upstream_neighbor,
    get_copy_tag,
    get_relative_child_coord,
    get_strand,
    merge_children_intervals,
)


def _feature(fid: str, start: int, end: int, seqid: str = 'chr1', strand: str = '+') -> Feature:
    return Feature(fid, 'gene', seqid, 'test', strand, start, end, '.', {'ID': [fid]})


@pytest.mark.parametrize(
    ('feature_id', 'original', 'tag'),
    [
        ('gene1_0', 'gene1', '_0'),
        ('gene_with_underscores_12', 'gene_with_underscores', '_12'),
        ('gene1_3_frag2', 'gene1', '_frag2'),
    ],
)
def test_copy_ids(feature_id: str, original: str, tag: str) -> None:
    assert convert_id_to_original(feature_id) == original
    assert get_copy_tag(feature_id) == tag


def test_count_overlap() -> None:
    assert count_overlap(1, 10, 5, 20) == 6
    assert count_overlap(1, 10, 10, 20) == 1
    assert count_overlap(1, 10, 11, 20) == 0


def test_relative_coordinates() -> None:
    parent = _feature('g', 100, 200)
    assert get_relative_child_coord(parent, 110, is_reverse=False) == 10
    assert get_relative_child_coord(parent, 110, is_reverse=True) == 90


def test_merge_children_intervals() -> None:
    children = [_feature('a', 50, 60), _feature('b', 1, 10), _feature('c', 5, 20)]
    assert merge_children_intervals(children) == [[1, 20], [50, 60]]
    assert merge_children_intervals([]) == []


def test_get_strand_flips_for_reverse_alignments() -> None:
    parent = _feature('g', 1, 10, strand='-')
    forward = AlignedSegment(1, 'g_0', 'chr1', 0, 9, 0, 9, False)
    reverse = AlignedSegment(1, 'g_0', 'chr1', 0, 9, 0, 9, True)
    assert get_strand(forward, parent) == '-'
    assert get_strand(reverse, parent) == '+'


def test_upstream_neighbor() -> None:
    order = ParentOrder.from_features(
        [_feature('b', 50, 60), _feature('a', 1, 10), _feature('c', 5, 10, seqid='chr2')]
    )
    assert order.ids == ['a', 'b', 'c']
    assert find_nonoverlapping_upstream_neighbor(order, 'a') is None
    assert find_nonoverlapping_upstream_neighbor(order, 'b') == 'a'
    # First feature on a new sequence has no upstream neighbour.
    assert find_nonoverlapping_upstream_neighbor(order, 'c') is None


def test_parent_order_duplicate_ids_use_first_occurrence() -> None:
    order = ParentOrder.from_features([_feature('x', 1, 2), _feature('x', 5, 6)])
    assert order.index('x') == 0
    assert order.index('missing') is None


def test_interval_index_matches_brute_force() -> None:
    rng = random.Random(1)
    intervals = []
    for i in range(300):
        start = rng.randint(0, 10_000)
        intervals.append(Interval(start, start + rng.randint(0, 500), f'k{i}', i))
    index = IntervalIndex(intervals)
    ordered = sorted(intervals, key=lambda iv: (iv.start, iv.end, iv.key))
    for _ in range(200):
        start = rng.randint(0, 10_500)
        end = start + rng.randint(0, 300)
        expected = [iv for iv in ordered if iv.start <= end and iv.end >= start]
        assert index.find(start, end) == expected


def test_interval_index_edges_are_inclusive() -> None:
    index = IntervalIndex([Interval(2, 3, 'a', None)])
    assert len(index.find(3, 3)) == 1
    assert len(index.find(0, 2)) == 1
    assert index.find(4, 9) == []
    assert IntervalIndex[None]().find(0, 10) == []
