"""Tests for splitting SAM alignments into gap-free aligned blocks."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import pytest

from liftoff.align.blocks import get_aligned_blocks, parse_alignment
from liftoff.io.sam import SamRecord, parse_sam_line
from liftoff.models import Feature, FeatureHierarchy, LiftoverType

from ..helpers import FeatureFactory, SamLineFactory

RecordFactory = Callable[..., SamRecord]


@pytest.fixture
def hierarchy(make_feature: FeatureFactory) -> FeatureHierarchy:
    """gene1 (101-150) with exons at query offsets 0-9 and 30-49."""
    return FeatureHierarchy(
        parents={'gene1': make_feature('gene1', 'gene', 101, 150)},
        intermediates={},
        children={
            'gene1': [
                make_feature('exon1', 'exon', 101, 110, parent='gene1'),
                make_feature('exon2', 'exon', 131, 150, parent='gene1'),
            ]
        },
    )


@pytest.fixture
def record(make_sam_line: SamLineFactory) -> RecordFactory:
    def factory(cigar: str, qname: str = 'gene1_0') -> SamRecord:
        return parse_sam_line(make_sam_line(qname=qname, rname='chrT', pos=1001, cigar=cigar))

    return factory


class TestGetAlignedBlocks:
    def test_splits_blocks_at_insertions(
        self, record: RecordFactory, hierarchy: FeatureHierarchy
    ) -> None:
        blocks = get_aligned_blocks(record('10=20I20='), 7, hierarchy, LiftoverType.CHROM_BY_CHROM)
        assert [(b.query_block_start, b.query_block_end) for b in blocks] == [(0, 9), (30, 49)]

    def test_target_coordinates_do_not_advance_over_insertions(
        self, record: RecordFactory, hierarchy: FeatureHierarchy
    ) -> None:
        blocks = get_aligned_blocks(record('10=20I20='), 7, hierarchy, LiftoverType.CHROM_BY_CHROM)
        assert [(b.reference_block_start, b.reference_block_end) for b in blocks] == [
            (1000, 1009),
            (1010, 1029),
        ]

    def test_records_mismatch_positions_per_block(
        self, record: RecordFactory, hierarchy: FeatureHierarchy
    ) -> None:
        blocks = get_aligned_blocks(
            record('5=1X4=20I20='), 7, hierarchy, LiftoverType.CHROM_BY_CHROM
        )
        assert [b.mismatches for b in blocks] == [[5], []]

    def test_blocks_share_the_alignment_id(
        self, record: RecordFactory, hierarchy: FeatureHierarchy
    ) -> None:
        blocks = get_aligned_blocks(record('10=20I20='), 7, hierarchy, LiftoverType.CHROM_BY_CHROM)
        assert {b.aln_id for b in blocks} == {7}

    def test_drops_blocks_outside_child_features(
        self, record: RecordFactory, hierarchy: FeatureHierarchy
    ) -> None:
        # After the 10-base hard clip, the first block covers query 10-24 (intron only).
        blocks = get_aligned_blocks(
            record('10H15=5D25='), 1, hierarchy, LiftoverType.CHROM_BY_CHROM
        )
        assert [(b.query_block_start, b.query_block_end) for b in blocks] == [(25, 49)]

    @pytest.mark.parametrize(
        ('cigar', 'expected_blocks'),
        [
            pytest.param('40=', 0, id='partial-alignment-rejected'),
            pytest.param('50=', 1, id='end-to-end-kept'),
        ],
    )
    def test_copies_require_end_to_end_alignment(
        self, record: RecordFactory, hierarchy: FeatureHierarchy, cigar: str, expected_blocks: int
    ) -> None:
        blocks = get_aligned_blocks(
            record(cigar, qname='gene1_1'), 1, hierarchy, LiftoverType.COPIES
        )
        assert len(blocks) == expected_blocks


class TestParseAlignment:
    @pytest.fixture
    def sam_file(
        self,
        tmp_path: Path,
        make_sam_line: SamLineFactory,
        hierarchy: FeatureHierarchy,
        make_feature: FeatureFactory,
    ) -> Path:
        """Two alignments of gene1 and an unmapped gene2."""
        hierarchy.parents['gene2'] = make_feature('gene2', 'gene', 1, 50)
        hierarchy.children['gene2'] = [make_feature('exon3', 'exon', 1, 50, parent='gene2')]
        lines = [
            make_sam_line('gene1', 0, 'chrT', 1001, '50='),
            make_sam_line('gene1', 256, 'chrT', 5001, '50='),
            make_sam_line('gene2', 4, '*', 0, '*'),
        ]
        path = tmp_path / 'alignments.sam'
        path.write_text('\n'.join(lines) + '\n')
        return path

    def test_copies_stage_numbers_each_alignment(
        self, sam_file: Path, hierarchy: FeatureHierarchy
    ) -> None:
        blocks = parse_alignment(sam_file, hierarchy, [], LiftoverType.COPIES)
        assert list(blocks) == ['gene1_1', 'gene1_2']

    def test_other_stages_group_alignments_as_copy_zero(
        self, sam_file: Path, hierarchy: FeatureHierarchy
    ) -> None:
        blocks = parse_alignment(sam_file, hierarchy, [], LiftoverType.CHROM_BY_CHROM)
        assert [b.aln_id for b in blocks['gene1_0']] == [1, 2]

    def test_unmapped_records_are_reported(
        self, sam_file: Path, hierarchy: FeatureHierarchy
    ) -> None:
        unmapped: list[Feature] = []
        parse_alignment(sam_file, hierarchy, unmapped, LiftoverType.CHROM_BY_CHROM)
        assert [feature.id for feature in unmapped] == ['gene2']

    def test_alignments_without_child_overlap_are_reported_unmapped(
        self, tmp_path: Path, make_sam_line: SamLineFactory, hierarchy: FeatureHierarchy
    ) -> None:
        sam = tmp_path / 'intron_only.sam'
        sam.write_text(make_sam_line('gene1', 0, 'chrT', 1, '10H15=') + '\n')
        unmapped: list[Feature] = []

        blocks = parse_alignment(sam, hierarchy, unmapped, LiftoverType.CHROM_BY_CHROM)

        assert (blocks, [f.id for f in unmapped]) == ({}, ['gene1'])
