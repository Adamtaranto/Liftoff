"""Tests for GFF3/GTF output formatting."""

from __future__ import annotations

from pathlib import Path

import pytest

from liftoff.io.gff_writer import edit_copy_ids, make_gff_line, make_gtf_line, write_new_gff
from liftoff.models import Feature

from ..helpers import FeatureFactory


@pytest.fixture
def cds(make_feature: FeatureFactory) -> Feature:
    feature = make_feature('cds1', 'CDS', 10, 20, strand='-', frame='2', parent='t1')
    feature.attributes.update({'Dbxref': ['a', 'b'], 'copy_id': ['x'], 'empty': []})
    return feature


class TestMakeGffLine:
    def test_columns(self, cds: Feature) -> None:
        columns = make_gff_line(cds.attributes, cds).split('\t')
        assert columns[:8] == ['chr1', 'test', 'CDS', '10', '20', '.', '-', '2']

    def test_id_is_first_attribute(self, cds: Feature) -> None:
        attributes = make_gff_line(cds.attributes, cds).split('\t')[8]
        assert attributes.startswith('ID=cds1;')

    def test_multiple_values_are_comma_separated(self, cds: Feature) -> None:
        assert 'Dbxref=a,b' in make_gff_line(cds.attributes, cds)

    def test_internal_copy_id_is_hidden(self, cds: Feature) -> None:
        assert 'copy_id' not in make_gff_line(cds.attributes, cds)

    def test_empty_attribute_has_no_value(self, cds: Feature) -> None:
        assert make_gff_line(cds.attributes, cds).endswith(';empty=')


def test_gtf_line(make_feature: FeatureFactory) -> None:
    exon = make_feature('exon1', 'exon', 1, 50, parent='gene1', transcript_id=['t1'])
    assert make_gtf_line(exon.attributes, exon) == (
        'chr1\ttest\texon\t1\t50\t.\t+\t.\tID "exon1"; Parent "gene1"; transcript_id "t1"; '
    )


class TestEditCopyIds:
    @pytest.fixture
    def exon_copy(self, make_feature: FeatureFactory) -> Feature:
        return make_feature(
            'exon1', 'exon', 1, 50, parent='gene1', transcript_id=['t1'], extra_copy_number=['2']
        )

    @pytest.mark.parametrize(
        ('key', 'value'),
        [
            pytest.param('ID', 'exon1_2', id='ID'),
            pytest.param('Parent', 'gene1_2', id='Parent'),
            pytest.param('transcript_id', 't1_2', id='suffix-_id'),
        ],
    )
    def test_suffixes_identifiers(self, exon_copy: Feature, key: str, value: str) -> None:
        assert edit_copy_ids(exon_copy)[key] == [value]

    def test_leaves_original_attributes_untouched(self, exon_copy: Feature) -> None:
        edit_copy_ids(exon_copy)
        assert exon_copy.attributes['ID'] == ['exon1']


class TestWriteNewGff:
    @pytest.fixture
    def lifted(self, make_feature: FeatureFactory) -> dict[str, list[Feature]]:
        """Two copies of gene1; the extra copy at position 1 maps poorly."""

        def gene_copy(copy_id: str, start: int, coverage: str, identity: str) -> list[Feature]:
            gene = make_feature(
                'gene1',
                'gene',
                start,
                start + 99,
                copy_num_ID=[copy_id],
                coverage=[coverage],
                sequence_ID=[identity],
            )
            exon = make_feature('exon1', 'exon', start, start + 49, parent='gene1')
            return [gene, exon]

        return {
            'gene1_0': gene_copy('gene1_0', 1000, '1.0', '1.0'),
            'gene1_1': gene_copy('gene1_1', 1, '0.4', '0.3'),
        }

    @pytest.fixture
    def written(self, tmp_path: Path, lifted: dict[str, list[Feature]]) -> list[str]:
        output = tmp_path / 'out.gff3'
        write_new_gff(lifted, str(output), 'gff3', 0.5, 0.5, command_line='liftoff test')
        return output.read_text().splitlines()

    def test_writes_header(self, written: list[str]) -> None:
        assert written[0] == '##gff-version 3'
        assert written[2] == '# liftoff test'

    def test_features_are_sorted_by_position(self, written: list[str]) -> None:
        starts = [int(line.split('\t')[3]) for line in written[3:]]
        assert starts == [1, 1, 1000, 1000]

    def test_copy_numbers_follow_copy_id_order(self, written: list[str]) -> None:
        copy_numbers = [
            line.split('extra_copy_number=')[1].split(';')[0]
            for line in written[3:]
            if '\tgene\t' in line
        ]
        # gene1_1 (position 1) is numbered after gene1_0 despite being written first.
        assert copy_numbers == ['1', '0']

    def test_low_quality_copy_is_flagged(self, written: list[str]) -> None:
        assert written[3].endswith('partial_mapping=True;low_identity=True')

    def test_extra_copy_children_have_suffixed_ids(self, written: list[str]) -> None:
        assert written[4].split('\t')[8] == 'ID=exon1_1;Parent=gene1_1;extra_copy_number=1'

    def test_writes_to_stdout(
        self, lifted: dict[str, list[Feature]], capsys: pytest.CaptureFixture[str]
    ) -> None:
        write_new_gff(lifted, 'stdout', 'gff3', 0.5, 0.5, command_line='liftoff test')
        assert capsys.readouterr().out.startswith('##gff-version 3\n')
