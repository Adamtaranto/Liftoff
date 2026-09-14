"""Tests for polishing helpers and polished-gene selection."""

from __future__ import annotations

from pathlib import Path

import pytest

from liftoff.io.fasta import write_gene_sequences
from liftoff.models import Feature, LiftoverType
from liftoff.pipeline import _polished_is_better
from liftoff.polish import add_splice_sites, find_overlapping_exon_groups


def _feature(fid: str, ftype: str, start: int, end: int) -> Feature:
    return Feature(fid, ftype, 'chr1', 'test', '+', start, end, '.', {'ID': [fid]})


def test_add_splice_sites_does_not_modify_exons() -> None:
    gene = _feature('g', 'gene', 100, 200)
    exons = [
        _feature('e1', 'exon', 101, 120),  # one base inside the gene: no 5' splice site
        _feature('e2', 'exon', 150, 170),
        _feature('e3', 'exon', 180, 199),  # one base inside the gene end: no 3' site
    ]
    before = [(e.start, e.end) for e in exons]
    spans, sites = add_splice_sites(exons, gene)
    assert [(e.start, e.end) for e in exons] == before
    assert spans == [('e1', 101, 122), ('e2', 148, 172), ('e3', 178, 199)]
    assert sites == [[121, 122], [148, 149], [171, 172], [178, 179]]
    # Repeated polishing of the same gene gives the same result.
    assert add_splice_sites(exons, gene) == (spans, sites)


def test_exon_groups_use_extended_spans() -> None:
    spans = [('e1', 101, 122), ('e2', 120, 140), ('e3', 178, 199)]
    assert find_overlapping_exon_groups([[101, 140], [178, 199]], spans) == [['e1', 'e2'], ['e3']]


def _attrs(orfs: str | None, identity: str, coverage: str) -> dict[str, list[str]]:
    attributes = {'sequence_ID': [identity], 'coverage': [coverage]}
    if orfs is not None:
        attributes['valid_ORFs'] = [orfs]
    return attributes


@pytest.mark.parametrize(
    ('original', 'polished', 'better'),
    [
        (_attrs('0', '1.0', '1.0'), _attrs('1', '0.9', '0.9'), True),  # more valid ORFs
        (_attrs('10', '0.9', '1.0'), _attrs('9', '1.0', '1.0'), False),  # numeric, not '9' > '10'
        (_attrs('1', '0.99', '1.0'), _attrs('1', '1.0', '1.0'), True),  # 1.0 > 0.99 numerically
        (_attrs('1', '1.0', '0.9'), _attrs('1', '0.99', '1.0'), False),  # identity before coverage
        (_attrs('1', '0.99', '0.9'), _attrs('1', '0.99', '1.0'), True),  # coverage breaks ties
        (_attrs('1', '0.99', '1.0'), _attrs('1', '0.99', '1.0'), False),  # equal: keep original
        (_attrs('1', '0.9', '0.9'), _attrs(None, '1.0', '1.0'), False),  # CDS lost in polishing
    ],
)
def test_polished_is_better(
    original: dict[str, list[str]], polished: dict[str, list[str]], better: bool
) -> None:
    assert _polished_is_better(original, polished) is better


def test_flank_is_not_applied_cumulatively(tmp_path: Path) -> None:
    fasta = tmp_path / 'ref.fa'
    fasta.write_text('>chr1\n' + 'ACGT' * 100 + '\n')
    gene = _feature('g', 'gene', 101, 200)
    for stage in (LiftoverType.CHROM_BY_CHROM, LiftoverType.UNMAPPED, LiftoverType.COPIES):
        write_gene_sequences({'g': gene}, [str(fasta)], str(fasta), tmp_path, stage, 0.1)
        assert (gene.start, gene.end) == (91, 210)
    assert gene.original_bounds == (101, 200)
    assert gene.annotated_bounds == (101, 200)
