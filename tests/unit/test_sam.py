"""Tests for the pure-Python SAM reader, which mirrors pysam's semantics."""

from __future__ import annotations

from pathlib import Path

import pytest

from liftoff.errors import InputError
from liftoff.io.sam import CigarOp, parse_cigar, parse_sam_line, read_sam

from ..helpers import SamLineFactory


class TestParseCigar:
    def test_parses_all_operation_codes(self) -> None:
        assert parse_cigar('1M2I3D4N5S6H7P8=9X10B') == [
            (CigarOp.MATCH, 1),
            (CigarOp.INSERTION, 2),
            (CigarOp.DELETION, 3),
            (CigarOp.REF_SKIP, 4),
            (CigarOp.SOFT_CLIP, 5),
            (CigarOp.HARD_CLIP, 6),
            (CigarOp.PAD, 7),
            (CigarOp.EQUAL, 8),
            (CigarOp.DIFF, 9),
            (CigarOp.BACK, 10),
        ]

    def test_star_means_no_operations(self) -> None:
        assert parse_cigar('*') == []

    @pytest.mark.parametrize('cigar', ['10', '10Z', '=10', '5M3', 'M'])
    def test_rejects_malformed_cigar(self, cigar: str) -> None:
        with pytest.raises(InputError, match='invalid CIGAR'):
            parse_cigar(cigar)


class TestParseSamLine:
    def test_reads_query_name(self, make_sam_line: SamLineFactory) -> None:
        assert parse_sam_line(make_sam_line(qname='geneX')).query_name == 'geneX'

    def test_reference_start_is_zero_based(self, make_sam_line: SamLineFactory) -> None:
        assert parse_sam_line(make_sam_line(pos=101)).reference_start == 100

    @pytest.mark.parametrize(
        ('flag', 'is_reverse', 'is_unmapped'),
        [
            pytest.param(0, False, False, id='forward'),
            pytest.param(16, True, False, id='reverse'),
            pytest.param(4, False, True, id='unmapped'),
            pytest.param(256 | 16, True, False, id='secondary-reverse'),
        ],
    )
    def test_flag_properties(
        self, make_sam_line: SamLineFactory, flag: int, is_reverse: bool, is_unmapped: bool
    ) -> None:
        record = parse_sam_line(make_sam_line(flag=flag))
        assert (record.is_reverse, record.is_unmapped) == (is_reverse, is_unmapped)

    def test_star_reference_name_is_none(self, make_sam_line: SamLineFactory) -> None:
        record = parse_sam_line(make_sam_line(flag=4, rname='*', pos=0, cigar='*'))
        assert record.reference_name is None

    def test_rejects_short_lines(self) -> None:
        with pytest.raises(InputError, match='fields'):
            parse_sam_line('a\tb\tc')

    @pytest.mark.parametrize(
        ('cigar', 'seq', 'start', 'end'),
        [
            pytest.param('10=', 'A' * 10, 0, 10, id='no-clipping'),
            pytest.param('3S10=2S', 'A' * 15, 3, 13, id='soft-clips'),
            pytest.param('5H3S10=2S', 'A' * 15, 3, 13, id='leading-hard-clip-ignored'),
            pytest.param('3S10=2S7H', 'A' * 15, 3, 13, id='trailing-hard-clip-ignored'),
            pytest.param('3S10=2I2S', '*', 3, 15, id='no-seq-counts-leading-soft-clip-only'),
            pytest.param('5H10=1X', '*', 0, 11, id='no-seq-with-hard-clip'),
        ],
    )
    def test_query_alignment_bounds(
        self, make_sam_line: SamLineFactory, cigar: str, seq: str, start: int, end: int
    ) -> None:
        record = parse_sam_line(make_sam_line(cigar=cigar, seq=seq))
        assert (record.query_alignment_start, record.query_alignment_end) == (start, end)


class TestReadSam:
    def test_skips_header_and_blank_lines(
        self, tmp_path: Path, make_sam_line: SamLineFactory
    ) -> None:
        sam = tmp_path / 'alignments.sam'
        lines = [
            '@HD\tVN:1.6',
            '@SQ\tSN:chr1\tLN:1000',
            make_sam_line(),
            '',
            make_sam_line(flag=16),
        ]
        sam.write_text('\n'.join(lines) + '\n')

        records = list(read_sam(sam))

        assert [record.is_reverse for record in records] == [False, True]

    def test_missing_file_raises_input_error(self, tmp_path: Path) -> None:
        with pytest.raises(InputError, match='not found'):
            list(read_sam(tmp_path / 'missing.sam'))
