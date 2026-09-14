"""Tests for the pure-Python SAM reader (pysam-compatible semantics)."""

from __future__ import annotations

from pathlib import Path

import pytest

from liftoff.errors import InputError
from liftoff.io.sam import CigarOp, parse_cigar, parse_sam_line, read_sam


def _line(
    flag: int = 0, rname: str = 'chr1', pos: int = 101, cigar: str = '10=', seq: str = '*'
) -> str:
    return '\t'.join(['gene_1', str(flag), rname, str(pos), '60', cigar, '*', '0', '0', seq, '*'])


def test_parse_cigar() -> None:
    assert parse_cigar('5H10=1X2I3D4S') == [
        (CigarOp.HARD_CLIP, 5),
        (CigarOp.EQUAL, 10),
        (CigarOp.DIFF, 1),
        (CigarOp.INSERTION, 2),
        (CigarOp.DELETION, 3),
        (CigarOp.SOFT_CLIP, 4),
    ]
    assert parse_cigar('*') == []


@pytest.mark.parametrize('bad', ['10', '10Z', '=10', '5M3'])
def test_invalid_cigar(bad: str) -> None:
    with pytest.raises(InputError):
        parse_cigar(bad)


def test_basic_fields() -> None:
    record = parse_sam_line(_line(flag=16, pos=101))
    assert record.query_name == 'gene_1'
    assert record.reference_name == 'chr1'
    assert record.reference_start == 100
    assert record.is_reverse
    assert not record.is_unmapped


def test_unmapped_record() -> None:
    record = parse_sam_line(_line(flag=4, rname='*', pos=0, cigar='*'))
    assert record.is_unmapped
    assert record.reference_name is None
    assert record.cigartuples == []


@pytest.mark.parametrize(
    ('cigar', 'seq', 'start', 'end'),
    [
        ('10=', 'A' * 10, 0, 10),
        ('3S10=2S', 'A' * 15, 3, 13),
        ('5H3S10=2S', 'A' * 15, 3, 13),  # hard clips are not part of SEQ
        ('3S10=2S7H', 'A' * 15, 3, 13),
        # Without SEQ the end is derived from the CIGAR; only leading soft
        # clips count (matching pysam).
        ('3S10=2I2S', '*', 3, 15),
        ('5H10=1X', '*', 0, 11),
    ],
)
def test_query_alignment_bounds(cigar: str, seq: str, start: int, end: int) -> None:
    record = parse_sam_line(_line(cigar=cigar, seq=seq))
    assert record.query_alignment_start == start
    assert record.query_alignment_end == end


def test_read_sam_skips_headers(tmp_path: Path) -> None:
    sam = tmp_path / 'x.sam'
    sam.write_text(
        '@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n' + _line() + '\n\n' + _line(flag=16) + '\n'
    )
    records = list(read_sam(sam))
    assert [r.is_reverse for r in records] == [False, True]


def test_read_sam_missing_file(tmp_path: Path) -> None:
    with pytest.raises(InputError, match='not found'):
        list(read_sam(tmp_path / 'missing.sam'))


def test_short_line_rejected() -> None:
    with pytest.raises(InputError, match='fields'):
        parse_sam_line('a\tb\tc')
