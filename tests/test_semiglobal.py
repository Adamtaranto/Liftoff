"""The NumPy port must reproduce parasail's sg_dx_trace_scan_sat exactly."""

from __future__ import annotations

import gzip
import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from liftoff.align.semiglobal import ScoringMatrix, sg_dx_trace
from liftoff.polish import make_scoring_matrix


def _load_cases(path: Path) -> list[dict[str, Any]]:
    with gzip.open(path, 'rt') as handle:
        return json.load(handle)  # type: ignore[no-any-return]


@pytest.mark.parametrize('corpus', ['sg_dx_liftoff_chr1.json.gz', 'sg_dx_synthetic.json.gz'])
def test_matches_parasail_corpus(data_dir: Path, corpus: str) -> None:
    matrix = make_scoring_matrix(3)
    cases = _load_cases(data_dir / 'parasail' / corpus)
    assert cases
    for index, case in enumerate(cases):
        result = sg_dx_trace(case['s1'], case['s2'], case['open'], case['extend'], matrix)
        observed = (
            result.score,
            result.end_query,
            result.end_ref,
            result.query_traceback,
            result.ref_traceback,
        )
        expected = (case['score'], case['end_query'], case['end_ref'], case['query'], case['ref'])
        assert observed == expected, f'case {index} ({case.get("tag", "liftoff")}) differs'


def test_scoring_matrix_matches_parasail_layout() -> None:
    matrix = make_scoring_matrix(3)
    # Alphabet plus a catch-all index for unknown characters.
    assert matrix.scores.shape == (10, 10)
    a, lower_a, star, unknown = 0, 4, 8, 9
    assert matrix.scores[a, a] == 2
    assert matrix.scores[a, lower_a] == 3
    assert matrix.scores[lower_a, a] == 3
    assert matrix.scores[a, 1] == -4
    assert matrix.scores[star, star] == 2
    assert (matrix.scores[unknown] == 0).all()
    np.testing.assert_array_equal(matrix.encode('ACgtN*'), [0, 1, 6, 7, unknown, star])


def test_case_insensitive_matrix_shares_indices() -> None:
    matrix = ScoringMatrix.create('ACGT', 1, -1)
    np.testing.assert_array_equal(matrix.encode('acgtACGT'), [0, 1, 2, 3, 0, 1, 2, 3])


def test_query_must_align_end_to_end_but_database_ends_are_free() -> None:
    matrix = ScoringMatrix.create('ACGT', 2, -4)
    result = sg_dx_trace('ACGT', 'TTTACGTTTT', 10, 1, matrix)
    assert result.score == 8
    assert result.query_traceback == '---ACGT---'
    assert result.ref_traceback == 'TTTACGTTTT'
    assert result.end_ref == 6


def test_gap_costs_open_plus_extensions() -> None:
    matrix = ScoringMatrix.create('ACGT', 5, -4)
    # One 2-base deletion from the database: 8 matches * 5 - (open + extend).
    result = sg_dx_trace('AAAACCGGGG', 'AAAAGGGG', 3, 1, matrix)
    assert result.score == 8 * 5 - 4
    assert result.query_traceback == 'AAAACCGGGG'
    assert result.ref_traceback == 'AAAA--GGGG'


@pytest.mark.parametrize(('s1', 's2'), [('', 'ACGT'), ('ACGT', '')])
def test_empty_sequences_rejected(s1: str, s2: str) -> None:
    with pytest.raises(ValueError, match='empty'):
        sg_dx_trace(s1, s2, 1, 1, ScoringMatrix.create('ACGT', 1, -1))


def test_negative_penalties_rejected() -> None:
    with pytest.raises(ValueError, match='negative'):
        sg_dx_trace('A', 'A', -1, 1, ScoringMatrix.create('ACGT', 1, -1))
