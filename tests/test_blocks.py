"""Tests for splitting alignments into gap-free blocks."""

from __future__ import annotations

from pathlib import Path

from liftoff.align.blocks import get_aligned_blocks, parse_alignment
from liftoff.io.sam import parse_sam_line
from liftoff.models import Feature, FeatureHierarchy, LiftoverType


def _feature(fid: str, ftype: str, start: int, end: int, parent: str | None = None) -> Feature:
    attributes = {'ID': [fid]}
    if parent:
        attributes['Parent'] = [parent]
    return Feature(fid, ftype, 'chr1', 'test', '+', start, end, '.', attributes)


def _hierarchy() -> FeatureHierarchy:
    gene = _feature('gene1', 'gene', 101, 150)
    exons = [
        _feature('exon1', 'exon', 101, 110, 'gene1'),
        _feature('exon2', 'exon', 131, 150, 'gene1'),
    ]
    return FeatureHierarchy({'gene1': gene}, {}, {'gene1': exons})


def _record(cigar: str, flag: int = 0, pos: int = 1001, name: str = 'gene1_0'):  # type: ignore[no-untyped-def]
    fields = [name, str(flag), 'chrT', str(pos), '60', cigar, '*', '0', '0', '*', '*']
    return parse_sam_line('\t'.join(fields))


def test_blocks_split_at_gaps_and_record_mismatches() -> None:
    # 10 aligned (1 mismatch), 20-base deletion from the query (intron), 20 aligned.
    record = _record('5=1X4=20I20=')
    blocks = get_aligned_blocks(record, 7, _hierarchy(), LiftoverType.CHROM_BY_CHROM)
    assert [(b.query_block_start, b.query_block_end) for b in blocks] == [(0, 9), (30, 49)]
    assert [(b.reference_block_start, b.reference_block_end) for b in blocks] == [
        (1000, 1009),
        (1010, 1029),
    ]
    assert blocks[0].mismatches == [5]
    assert blocks[1].mismatches == []
    assert {b.aln_id for b in blocks} == {7}


def test_blocks_without_children_are_dropped() -> None:
    # The first block lies entirely in the intron (query 10-29) and is removed.
    record = _record('10H15=5D25=')
    blocks = get_aligned_blocks(record, 1, _hierarchy(), LiftoverType.CHROM_BY_CHROM)
    assert [(b.query_block_start, b.query_block_end) for b in blocks] == [(25, 49)]


def test_copies_require_end_to_end_alignments() -> None:
    partial = _record('40=', name='gene1_1')
    assert get_aligned_blocks(partial, 1, _hierarchy(), LiftoverType.COPIES) == []
    full = _record('50=', name='gene1_1')
    assert len(get_aligned_blocks(full, 1, _hierarchy(), LiftoverType.COPIES)) == 1


def test_parse_alignment_names_copies_and_collects_unmapped(tmp_path: Path) -> None:
    hierarchy = _hierarchy()
    gene2 = _feature('gene2', 'gene', 1, 50)
    hierarchy.parents['gene2'] = gene2
    hierarchy.children['gene2'] = [_feature('exon3', 'exon', 1, 50, 'gene2')]
    sam = tmp_path / 'aln.sam'
    rows = [
        ['gene1', '0', 'chrT', '1001', '60', '50=', '*', '0', '0', '*', '*'],
        ['gene1', '256', 'chrT', '5001', '60', '50=', '*', '0', '0', '*', '*'],
        ['gene2', '4', '*', '0', '0', '*', '*', '0', '0', '*', '*'],
    ]
    sam.write_text('\n'.join('\t'.join(r) for r in rows) + '\n')

    unmapped: list[Feature] = []
    copies = parse_alignment(sam, hierarchy, unmapped, LiftoverType.COPIES)
    assert list(copies) == ['gene1_1', 'gene1_2']
    assert [f.id for f in unmapped] == ['gene2']

    unmapped.clear()
    normal = parse_alignment(sam, hierarchy, unmapped, LiftoverType.CHROM_BY_CHROM)
    assert list(normal) == ['gene1_0']
    assert [b.aln_id for b in normal['gene1_0']] == [1, 2]
