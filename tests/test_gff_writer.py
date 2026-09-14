"""Tests for GFF3/GTF output formatting."""

from __future__ import annotations

from pathlib import Path

from liftoff.io.gff_writer import edit_copy_ids, make_gff_line, make_gtf_line, write_new_gff
from liftoff.models import Feature


def _gene(copy_id: str, start: int, coverage: str = '1.0', identity: str = '1.0') -> Feature:
    return Feature(
        'gene1',
        'gene',
        'chr1',
        'Liftoff',
        '+',
        start,
        start + 99,
        '.',
        {
            'ID': ['gene1'],
            'Name': ['ABC1'],
            'copy_num_ID': [copy_id],
            'coverage': [coverage],
            'sequence_ID': [identity],
        },
    )


def _exon(start: int) -> Feature:
    return Feature(
        'exon1',
        'exon',
        'chr1',
        'Liftoff',
        '+',
        start,
        start + 49,
        '.',
        {'ID': ['exon1'], 'Parent': ['gene1'], 'transcript_id': ['t1']},
    )


def test_gff_line_puts_id_first_and_hides_copy_id() -> None:
    feature = Feature(
        'cds1',
        'CDS',
        'chr1',
        'Liftoff',
        '-',
        10,
        20,
        '2',
        {'Parent': ['t1'], 'ID': ['cds1'], 'Dbxref': ['a', 'b'], 'copy_id': ['x'], 'empty': []},
    )
    line = make_gff_line(feature.attributes, feature)
    assert line == 'chr1\tLiftoff\tCDS\t10\t20\t.\t-\t2\tID=cds1;Parent=t1;Dbxref=a,b;empty='


def test_gtf_line() -> None:
    feature = _exon(1)
    line = make_gtf_line(feature.attributes, feature)
    assert line == (
        'chr1\tLiftoff\texon\t1\t50\t.\t+\t.\tID "exon1"; Parent "gene1"; transcript_id "t1"; '
    )


def test_edit_copy_ids_suffixes_identifiers() -> None:
    feature = _exon(1)
    feature.attributes['extra_copy_number'] = ['2']
    edited = edit_copy_ids(feature)
    assert edited['ID'] == ['exon1_2']
    assert edited['Parent'] == ['gene1_2']
    assert edited['transcript_id'] == ['t1_2']
    assert feature.attributes['ID'] == ['exon1']  # original untouched


def test_write_new_gff_numbers_copies_and_flags_partial(tmp_path: Path) -> None:
    lifted = {
        'gene1_0': [_gene('gene1_0', 1000), _exon(1000)],
        'gene1_1': [_gene('gene1_1', 1, coverage='0.4', identity='0.3'), _exon(1)],
    }
    output = tmp_path / 'out.gff3'
    write_new_gff(lifted, str(output), 'gff3', 0.5, 0.5, command_line='liftoff test')
    lines = output.read_text().splitlines()
    assert lines[0] == '##gff-version 3'
    assert lines[2] == '# liftoff test'
    body = lines[3:]
    # Sorted by position; the copy at position 1 was numbered second (ID order).
    assert body[0].split('\t')[8] == (
        'ID=gene1_1;Name=ABC1;coverage=0.4;sequence_ID=0.3;'
        'extra_copy_number=1;copy_num_ID=gene1_1;partial_mapping=True;low_identity=True'
    )
    assert body[1].split('\t')[8] == (
        'ID=exon1_1;Parent=gene1_1;transcript_id=t1_1;extra_copy_number=1'
    )
    assert body[2].split('\t')[8].endswith('extra_copy_number=0;copy_num_ID=gene1_0')
