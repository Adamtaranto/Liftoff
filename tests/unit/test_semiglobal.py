"""The NumPy semi-global aligner must reproduce parasail's sg_dx_trace_scan_sat."""

from __future__ import annotations

import gzip
import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from liftoff.align.semiglobal import ScoringMatrix, sg_dx_trace
from liftoff.polish import make_scoring_matrix

CORPORA = ['sg_dx_liftoff_chr1.json.gz', 'sg_dx_synthetic.json.gz']


def load_corpus(path: Path) -> list[dict[str, Any]]:
    with gzip.open(path, 'rt') as handle:
        cases: list[dict[str, Any]] = json.load(handle)
    return cases


def pytest_generate_tests(metafunc: pytest.Metafunc) -> None:
    """Parametrize ``parasail_case`` with every recorded alignment."""
    if 'parasail_case' not in metafunc.fixturenames:
        return
    data = Path(__file__).parents[1] / 'data' / 'parasail'
    params = []
    for corpus in CORPORA:
        for index, case in enumerate(load_corpus(data / corpus)):
            tag = case.get('tag', 'liftoff')
            params.append(pytest.param(case, id=f'{corpus.split(".")[0]}-{index}-{tag}'))
    metafunc.parametrize('parasail_case', params)


@pytest.fixture(scope='module')
def polish_matrix() -> ScoringMatrix:
    return make_scoring_matrix(3)


def test_matches_parasail(parasail_case: dict[str, Any], polish_matrix: ScoringMatrix) -> None:
    result = sg_dx_trace(
        parasail_case['s1'],
        parasail_case['s2'],
        parasail_case['open'],
        parasail_case['extend'],
        polish_matrix,
    )

    assert (
        result.score,
        result.end_query,
        result.end_ref,
        result.query_traceback,
        result.ref_traceback,
    ) == (
        parasail_case['score'],
        parasail_case['end_query'],
        parasail_case['end_ref'],
        parasail_case['query'],
        parasail_case['ref'],
    )


class TestScoringMatrix:
    def test_has_catch_all_row_and_column(self, polish_matrix: ScoringMatrix) -> None:
        assert polish_matrix.scores.shape == (10, 10)

    @pytest.mark.parametrize(
        ('row', 'column', 'score'),
        [
            pytest.param(0, 0, 2, id='A-A-match'),
            pytest.param(0, 4, 3, id='A-a-coding-match'),
            pytest.param(4, 0, 3, id='a-A-coding-match-symmetric'),
            pytest.param(0, 1, -4, id='A-C-mismatch'),
            pytest.param(8, 8, 2, id='star-star'),
            pytest.param(9, 3, 0, id='unknown-scores-zero'),
        ],
    )
    def test_polish_matrix_scores(
        self, polish_matrix: ScoringMatrix, row: int, column: int, score: int
    ) -> None:
        assert polish_matrix.scores[row, column] == score

    def test_case_sensitive_encoding_maps_unknown_characters_last(
        self, polish_matrix: ScoringMatrix
    ) -> None:
        np.testing.assert_array_equal(polish_matrix.encode('ACgtN*'), [0, 1, 6, 7, 9, 8])

    def test_case_insensitive_encoding_shares_indices(self) -> None:
        matrix = ScoringMatrix.create('ACGT', 1, -1)
        np.testing.assert_array_equal(matrix.encode('acgtACGT'), [0, 1, 2, 3, 0, 1, 2, 3])

    def test_with_scores_returns_modified_copy(self) -> None:
        matrix = ScoringMatrix.create('ACGT', 1, -1)
        updated = matrix.with_scores([(0, 1, 7)])
        assert (updated.scores[0, 1], matrix.scores[0, 1]) == (7, -1)


class TestSemiglobalAlignment:
    @pytest.fixture
    def matrix(self) -> ScoringMatrix:
        return ScoringMatrix.create('ACGT', 2, -4)

    def test_database_end_gaps_are_free(self, matrix: ScoringMatrix) -> None:
        result = sg_dx_trace('ACGT', 'TTTACGTTTT', 10, 1, matrix)
        assert (result.score, result.end_ref) == (8, 6)

    def test_traceback_includes_unaligned_database_ends(self, matrix: ScoringMatrix) -> None:
        result = sg_dx_trace('ACGT', 'TTTACGTTTT', 10, 1, matrix)
        assert (result.query_traceback, result.ref_traceback) == ('---ACGT---', 'TTTACGTTTT')

    def test_gap_costs_open_plus_extensions(self) -> None:
        matrix = ScoringMatrix.create('ACGT', 5, -4)
        result = sg_dx_trace('AAAACCGGGG', 'AAAAGGGG', 3, 1, matrix)
        # 8 matches, one 2-base gap costing open + 1 extension.
        assert result.score == 8 * 5 - (3 + 1)

    def test_gap_is_placed_in_database_traceback(self) -> None:
        matrix = ScoringMatrix.create('ACGT', 5, -4)
        result = sg_dx_trace('AAAACCGGGG', 'AAAAGGGG', 3, 1, matrix)
        assert result.ref_traceback == 'AAAA--GGGG'

    @pytest.mark.parametrize(
        ('s1', 's2'),
        [pytest.param('', 'ACGT', id='empty-query'), pytest.param('ACGT', '', id='empty-database')],
    )
    def test_rejects_empty_sequences(self, matrix: ScoringMatrix, s1: str, s2: str) -> None:
        with pytest.raises(ValueError, match='empty'):
            sg_dx_trace(s1, s2, 1, 1, matrix)

    @pytest.mark.parametrize(
        ('open_penalty', 'extend_penalty'),
        [pytest.param(-1, 1, id='open'), pytest.param(1, -1, id='extend')],
    )
    def test_rejects_negative_penalties(
        self, matrix: ScoringMatrix, open_penalty: int, extend_penalty: int
    ) -> None:
        with pytest.raises(ValueError, match='negative'):
            sg_dx_trace('A', 'A', open_penalty, extend_penalty, matrix)
