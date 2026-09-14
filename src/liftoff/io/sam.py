"""Minimal pure-Python SAM reader.

Liftoff only needs a handful of fields from each alignment record, so this
module replaces the ``pysam`` dependency (a C extension that is unavailable in
Pyodide). Property semantics — in particular :attr:`SamRecord.query_alignment_start`
and :attr:`SamRecord.query_alignment_end` — mirror ``pysam.AlignedSegment`` so
that downstream coordinate arithmetic is unchanged.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from enum import IntEnum
from pathlib import Path
import re

from liftoff.errors import InputError


class CigarOp(IntEnum):
    """Numeric CIGAR operation codes as defined by the SAM specification."""

    MATCH = 0  # M
    INSERTION = 1  # I
    DELETION = 2  # D
    REF_SKIP = 3  # N
    SOFT_CLIP = 4  # S
    HARD_CLIP = 5  # H
    PAD = 6  # P
    EQUAL = 7  # =
    DIFF = 8  # X
    BACK = 9  # B


_CIGAR_CODES = {code: int(op) for code, op in zip('MIDNSHP=XB', CigarOp, strict=True)}
_CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=XB])')

#: SAM FLAG bits used by Liftoff.
FLAG_UNMAPPED = 0x4
FLAG_REVERSE = 0x10


def parse_cigar(cigar: str) -> list[tuple[int, int]]:
    """Convert a CIGAR string into ``(operation, length)`` tuples.

    Parameters
    ----------
    cigar : str
        CIGAR string, or ``"*"`` when unavailable.

    Returns
    -------
    list of tuple of int
        Operations as ``(CigarOp, length)`` pairs, like ``pysam``'s
        ``cigartuples``.

    Raises
    ------
    InputError
        If the string is not a valid CIGAR.

    Examples
    --------
    >>> parse_cigar('5H10=1X2I')
    [(5, 5), (7, 10), (8, 1), (1, 2)]
    """
    if cigar == '*':
        return []
    pairs = _CIGAR_RE.findall(cigar)
    # Every character must belong to a matched "<length><op>" pair.
    if sum(len(length) + 1 for length, _ in pairs) != len(cigar):
        raise InputError(f'invalid CIGAR string: {cigar!r}')
    return [(_CIGAR_CODES[code], int(length)) for length, code in pairs]


@dataclass(slots=True)
class SamRecord:
    """The subset of a SAM alignment record used by Liftoff.

    Attributes
    ----------
    query_name : str
        QNAME field. Mutable, because Liftoff appends copy suffixes.
    flag : int
        Bitwise FLAG field.
    reference_name : str or None
        RNAME field, ``None`` when ``"*"``.
    reference_start : int
        Zero-based leftmost mapping position (``POS - 1``).
    cigartuples : list of tuple of int
        Parsed CIGAR operations.
    query_length : int
        Length of the SEQ field (0 when SEQ is ``"*"``).
    """

    query_name: str
    flag: int
    reference_name: str | None
    reference_start: int
    cigartuples: list[tuple[int, int]]
    query_length: int

    @property
    def is_unmapped(self) -> bool:
        """Whether the read is flagged as unmapped.

        Returns
        -------
        bool
            ``True`` if FLAG bit 0x4 is set.
        """
        return bool(self.flag & FLAG_UNMAPPED)

    @property
    def is_reverse(self) -> bool:
        """Whether the read aligned to the reverse strand.

        Returns
        -------
        bool
            ``True`` if FLAG bit 0x10 is set.
        """
        return bool(self.flag & FLAG_REVERSE)

    @property
    def query_alignment_start(self) -> int:
        """Start of the aligned query portion, excluding soft clips.

        Leading hard clips are skipped (they are not part of SEQ), and leading
        soft clips are counted, matching ``pysam``.

        Returns
        -------
        int
            Zero-based offset into SEQ of the first aligned base.
        """
        start = 0
        for op, length in self.cigartuples:
            if op == CigarOp.HARD_CLIP:
                continue
            if op == CigarOp.SOFT_CLIP:
                start += length
            else:
                break
        return start

    @property
    def query_alignment_end(self) -> int:
        """End (exclusive) of the aligned query portion, excluding soft clips.

        When SEQ is absent (common for minimap2 secondary alignments) the query
        length is derived from the CIGAR, again matching ``pysam``.

        Returns
        -------
        int
            Zero-based exclusive offset into SEQ after the last aligned base.
        """
        cigar = self.cigartuples
        end = self.query_length
        if end == 0:
            # No sequence: sum query-consuming operations. Soft clips only
            # count while nothing else has been counted yet (i.e. leading).
            for op, length in cigar:
                if op in (CigarOp.MATCH, CigarOp.INSERTION, CigarOp.EQUAL, CigarOp.DIFF) or (
                    op == CigarOp.SOFT_CLIP and end == 0
                ):
                    end += length
            return end
        # Walk backwards over trailing clips (never examining the first op).
        for op, length in reversed(cigar[1:]):
            if op == CigarOp.SOFT_CLIP:
                end -= length
            elif op != CigarOp.HARD_CLIP:
                break
        return end


def parse_sam_line(line: str) -> SamRecord:
    """Parse one tab-separated SAM alignment line.

    Parameters
    ----------
    line : str
        A non-header SAM line.

    Returns
    -------
    SamRecord
        The parsed record.

    Raises
    ------
    InputError
        If the line has fewer than 11 fields.
    """
    fields = line.rstrip('\r\n').split('\t', 11)
    if len(fields) < 11:
        raise InputError(f'SAM record has {len(fields)} fields, expected at least 11: {line!r}')
    qname, flag, rname, pos, _mapq, cigar, _rnext, _pnext, _tlen, seq = fields[:10]
    return SamRecord(
        query_name=qname,
        flag=int(flag),
        reference_name=None if rname == '*' else rname,
        reference_start=int(pos) - 1,
        cigartuples=parse_cigar(cigar),
        query_length=0 if seq == '*' else len(seq),
    )


def read_sam(path: str | Path) -> Iterator[SamRecord]:
    """Iterate over alignment records of a SAM file in file order.

    Parameters
    ----------
    path : str or pathlib.Path
        SAM file to read. Header lines (starting with ``@``) are skipped.

    Yields
    ------
    SamRecord
        Each alignment record.

    Raises
    ------
    InputError
        If the file does not exist.
    """
    path = Path(path)
    if not path.is_file():
        raise InputError(f'SAM file not found: {path}')
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith('@'):
                continue
            yield parse_sam_line(line)
