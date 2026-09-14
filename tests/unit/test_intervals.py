"""Tests for the interval index that replaces ``interlap``."""

from __future__ import annotations

import random

import pytest

from liftoff.intervals import Interval, IntervalIndex


@pytest.fixture
def random_intervals(rng: random.Random) -> list[Interval[int]]:
    intervals = []
    for index in range(300):
        start = rng.randint(0, 10_000)
        intervals.append(Interval(start, start + rng.randint(0, 500), f'k{index}', index))
    return intervals


def test_find_matches_brute_force_in_sorted_order(
    random_intervals: list[Interval[int]], rng: random.Random
) -> None:
    index = IntervalIndex(random_intervals)
    ordered = sorted(random_intervals, key=lambda iv: (iv.start, iv.end, iv.key))
    for _ in range(200):
        start = rng.randint(0, 10_500)
        end = start + rng.randint(0, 300)
        expected = [iv for iv in ordered if iv.start <= end and iv.end >= start]
        assert index.find(start, end) == expected


@pytest.mark.parametrize(
    ('query', 'hits'),
    [
        pytest.param((3, 3), 1, id='touches-end'),
        pytest.param((0, 2), 1, id='touches-start'),
        pytest.param((4, 9), 0, id='after'),
        pytest.param((0, 1), 0, id='before'),
    ],
)
def test_bounds_are_inclusive(query: tuple[int, int], hits: int) -> None:
    index = IntervalIndex([Interval(2, 3, 'a', None)])
    assert len(index.find(*query)) == hits


def test_ties_are_ordered_by_key() -> None:
    index = IntervalIndex([Interval(1, 5, 'b', 2), Interval(1, 5, 'a', 1)])
    assert [iv.value for iv in index.find(1, 1)] == [1, 2]


def test_empty_index() -> None:
    index = IntervalIndex[None]()
    assert (len(index), index.find(0, 10)) == (0, [])


def test_iterates_in_sorted_order() -> None:
    index = IntervalIndex([Interval(9, 10, 'b', None), Interval(1, 2, 'a', None)])
    assert [iv.key for iv in index] == ['a', 'b']
