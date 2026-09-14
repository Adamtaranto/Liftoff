"""Semi-global affine-gap alignment with traceback (a port of parasail).

Liftoff previously used ``parasail.sg_dx_trace_scan_sat`` to re-align exons
while polishing. This module reimplements exactly that routine in NumPy so that
Liftoff has no compiled dependencies. Results — score, end position and the
gapped traceback strings, including how ties between equally scoring
alignments are broken — are identical to parasail's vectorised ("scan")
implementation.

Terminology follows parasail:

* ``s1`` is the *query* (rows ``i``), ``s2`` the *database* (columns ``j``).
* ``sg_dx``: gaps at the beginning and end of the **database** are free, i.e.
  the query must align end to end but may start and stop anywhere in the
  database. Leading/trailing query gaps are penalised.
* A gap of length ``L`` costs ``open + (L - 1) * extend``.
* ``E`` holds scores of alignments ending in a gap in the query (an insertion,
  moving along ``s2``), ``F`` those ending in a gap in the database (a
  deletion, moving along ``s1``), and ``H`` the best of ``E``, ``F`` and the
  diagonal.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

# Trace flags, bit-compatible with parasail's PARASAIL_* constants.
_INS = 1  # H came from E
_DEL = 2  # H came from F
_DIAG = 4  # H came from the diagonal
_DIAG_E = 8  # E was opened from H
_INS_E = 16  # E was extended from E
_DIAG_F = 32  # F was opened from H
_DEL_F = 64  # F was extended from F

# Far below any reachable score, but safely away from integer overflow.
_NEG_INF = np.int64(-(2**40))


@dataclass(frozen=True, slots=True)
class ScoringMatrix:
    """A substitution matrix over a byte alphabet.

    Mirrors ``parasail.matrix_create``: the matrix has one row/column per
    alphabet character plus a final catch-all row/column (scored 0) for
    characters not in the alphabet.

    Attributes
    ----------
    alphabet : str
        Characters with explicit scores.
    scores : numpy.ndarray
        ``(len(alphabet) + 1, len(alphabet) + 1)`` integer score table.
    mapper : numpy.ndarray
        Length-256 array mapping byte values to matrix indices.
    """

    alphabet: str
    scores: npt.NDArray[np.int64]
    mapper: npt.NDArray[np.intp]

    @classmethod
    def create(
        cls, alphabet: str, match: int, mismatch: int, *, case_sensitive: bool = False
    ) -> ScoringMatrix:
        """Build a simple match/mismatch matrix.

        Parameters
        ----------
        alphabet : str
            Characters with explicit scores.
        match : int
            Score on the diagonal.
        mismatch : int
            Score off the diagonal.
        case_sensitive : bool, default False
            If false, upper- and lower-case forms of a letter share an index.

        Returns
        -------
        ScoringMatrix
            The new matrix.
        """
        size = len(alphabet)
        scores = np.full((size + 1, size + 1), mismatch, dtype=np.int64)
        np.fill_diagonal(scores, match)
        # Unknown characters score zero against everything.
        scores[size, :] = 0
        scores[:, size] = 0
        mapper = np.full(256, size, dtype=np.intp)
        for index, char in enumerate(alphabet):
            if case_sensitive:
                mapper[ord(char)] = index
            else:
                mapper[ord(char.upper())] = index
                mapper[ord(char.lower())] = index
        return cls(alphabet, scores, mapper)

    def with_scores(self, updates: Iterable[tuple[int, int, int]]) -> ScoringMatrix:
        """Return a copy with individual entries replaced.

        Parameters
        ----------
        updates : iterable of tuple of int
            ``(row, column, score)`` triples.

        Returns
        -------
        ScoringMatrix
            The modified matrix.
        """
        scores = self.scores.copy()
        for row, column, score in updates:
            scores[row, column] = score
        return ScoringMatrix(self.alphabet, scores, self.mapper)

    def encode(self, sequence: str) -> npt.NDArray[np.intp]:
        """Map a sequence to matrix indices.

        Parameters
        ----------
        sequence : str
            Sequence to encode (non-Latin-1 characters are treated as unknown).

        Returns
        -------
        numpy.ndarray
            Matrix index for each character.
        """
        raw = np.frombuffer(sequence.encode('latin-1', errors='replace'), dtype=np.uint8)
        return self.mapper[raw]


@dataclass(frozen=True, slots=True)
class Alignment:
    """Result of :func:`sg_dx_trace`.

    Attributes
    ----------
    score : int
        Optimal alignment score.
    end_query : int
        Zero-based query index of the alignment end.
    end_ref : int
        Zero-based database index of the alignment end.
    query_traceback : str
        Gapped query, covering the whole query and database.
    ref_traceback : str
        Gapped database, aligned column-for-column with ``query_traceback``.
    """

    score: int
    end_query: int
    end_ref: int
    query_traceback: str
    ref_traceback: str


def sg_dx_trace(
    s1: str, s2: str, open_penalty: int, extend_penalty: int, matrix: ScoringMatrix
) -> Alignment:
    """Semi-global alignment with free end gaps in ``s2``, with traceback.

    Equivalent to ``parasail.sg_dx_trace_scan_sat(s1, s2, open, extend,
    matrix)`` followed by reading ``result.traceback.query`` and
    ``result.traceback.ref``.

    Parameters
    ----------
    s1 : str
        Query sequence; aligned end to end.
    s2 : str
        Database sequence; unaligned leading and trailing bases are free.
    open_penalty : int
        Gap opening penalty (non-negative); a length-1 gap costs this much.
    extend_penalty : int
        Gap extension penalty (non-negative).
    matrix : ScoringMatrix
        Substitution scores.

    Returns
    -------
    Alignment
        Score, end coordinates and traceback strings.

    Raises
    ------
    ValueError
        If either sequence is empty or a penalty is negative.
    """
    if not s1 or not s2:
        raise ValueError('sequences to align must not be empty')
    if open_penalty < 0 or extend_penalty < 0:
        raise ValueError('gap penalties must not be negative')
    trace, score, end_ref = _fill(s1, s2, open_penalty, extend_penalty, matrix)
    query_tb, ref_tb = _traceback(s1, s2, trace, len(s1) - 1, end_ref)
    return Alignment(score, len(s1) - 1, end_ref, query_tb, ref_tb)


def _fill(
    s1: str, s2: str, open_: int, extend: int, matrix: ScoringMatrix
) -> tuple[npt.NDArray[np.uint8], int, int]:
    """Fill the DP matrices column by column and record trace flags.

    The outer loop walks the database; each column is computed for all query
    positions at once. ``E`` depends only on the previous column, and the
    in-column ``F`` recurrence is solved with a prefix maximum (the "scan"
    trick), so no per-cell Python loop is needed.
    """
    m, n = len(s1), len(s2)
    profile = matrix.scores[matrix.encode(s1)]  # (m, alphabet) query profile
    db = matrix.encode(s2)
    rows = np.arange(m, dtype=np.int64)
    # Offsets for the prefix-max formulation of F:
    #   F[i] = max_{k < i} (Ht[k] - open - (i - 1 - k) * extend)
    # with the top boundary H[-1] = 0 (free leading database gaps) as k = -1.
    offsets = rows * extend

    trace = np.empty((m, n), dtype=np.uint8)
    # Left boundary (column -1): leading query gaps are penalised.
    h = -open_ - extend * rows
    e = h - open_
    diag_src = np.empty(m, dtype=np.int64)
    scan = np.empty(m + 1, dtype=np.int64)
    h_up = np.empty(m, dtype=np.int64)
    f_up = np.empty(m, dtype=np.int64)

    best_score = int(_NEG_INF)
    end_ref = n - 1
    for j in range(n):
        # E: gap in the query, coming from the previous column.
        e_open = h - open_
        e_ext = e - extend
        e = np.maximum(e_open, e_ext)
        e_flags = np.where(e_open > e_ext, _DIAG_E, _INS_E)

        # Diagonal: H[i-1][j-1] + score; H[-1][*] = 0 on the top boundary.
        diag_src[0] = 0
        diag_src[1:] = h[:-1]
        h_diag = diag_src + profile[:, db[j]]
        h_tmp = np.maximum(e, h_diag)

        # F: gap in the database, propagated down the column via prefix max.
        scan[0] = -extend  # boundary term for k = -1
        np.add(h_tmp, offsets, out=scan[1:])
        np.maximum.accumulate(scan, out=scan)
        f = scan[:-1] - open_ - offsets + extend
        h_new = np.maximum(h_tmp, f)

        # F flags compare opening from the final H above with extending F above.
        h_up[0] = 0
        h_up[1:] = h_new[:-1]
        f_up[0] = _NEG_INF
        f_up[1:] = f[:-1]
        f_flags = np.where(h_up - open_ > f_up - extend, _DIAG_F, _DEL_F)

        # H flags: prefer diagonal, then deletion (F), then insertion (E).
        h_flags = np.where(h_new == h_diag, _DIAG, np.where(h_new == f, _DEL, _INS))
        trace[:, j] = h_flags | e_flags | f_flags
        h = h_new

        # Free trailing database gaps: best score in the last query row,
        # earliest column wins ties.
        if h[-1] > best_score:
            best_score = int(h[-1])
            end_ref = j
    return trace, best_score, end_ref


def _traceback(s1: str, s2: str, trace: npt.NDArray[np.uint8], i: int, j: int) -> tuple[str, str]:
    """Walk trace flags back from ``(i, j)`` and build gapped strings."""
    query: list[str] = []
    ref: list[str] = []
    # Semi-global alignments include the free trailing database bases.
    for k in range(len(s2) - 1, j, -1):
        query.append('-')
        ref.append(s2[k])
    state = _DIAG
    while i >= 0 or j >= 0:
        if i < 0:
            # Query exhausted: remaining database bases are leading gaps.
            query.extend('-' * (j + 1))
            ref.extend(reversed(s2[: j + 1]))
            break
        if j < 0:
            # Database exhausted: remaining query bases are deletions.
            query.extend(reversed(s1[: i + 1]))
            ref.extend('-' * (i + 1))
            break
        flags = int(trace[i, j])
        if state == _DIAG:
            if flags & _DIAG:
                query.append(s1[i])
                ref.append(s2[j])
                i -= 1
                j -= 1
            elif flags & _INS:
                state = _INS
            elif flags & _DEL:
                state = _DEL
            else:  # pragma: no cover - every cell carries an H flag
                break
        elif state == _INS:
            query.append('-')
            ref.append(s2[j])
            j -= 1
            state = _DIAG if flags & _DIAG_E else _INS
        else:
            query.append(s1[i])
            ref.append('-')
            i -= 1
            state = _DIAG if flags & _DIAG_F else _DEL
    query.reverse()
    ref.reverse()
    return ''.join(query), ''.join(ref)
