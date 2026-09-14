"""Tests for building the feature database and classifying reference features."""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path

import gffutils
import pytest

from liftoff.errors import GffSyntaxError, InputError
from liftoff.io.gff_db import (
    annotation_format,
    build_database,
    get_feature_order,
    separate_parents_and_children,
)
from liftoff.models import FeatureHierarchy
from liftoff.utils import ParentOrder

ANNOTATION = """\
##gff-version 3
chr1\tsrc\tgene\t500\t800\t.\t+\t.\tID=gene2;Name=B
chr1\tsrc\tmRNA\t500\t800\t.\t+\t.\tID=tx2;Parent=gene2
chr1\tsrc\texon\t500\t800\t.\t+\t.\tID=exon3;Parent=tx2
chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=gene1;Name=A%3Bescaped
chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=tx1;Parent=gene1
chr1\tsrc\texon\t1\t40\t.\t+\t.\tID=exon1;Parent=tx1
chr1\tsrc\tCDS\t10\t40\t.\t+\t0\tID=cds1;Parent=tx1
chr1\tsrc\texon\t60\t100\t.\t+\t.\tID=exon2;Parent=tx1
chr1\tsrc\tintron\t41\t59\t.\t+\t.\tID=intron1;Parent=tx1
chr1\tsrc\trepeat_region\t300\t400\t.\t+\t.\tID=repeat1
"""


@pytest.fixture
def gff(tmp_path: Path) -> Path:
    path = tmp_path / 'annotation.gff3'
    path.write_text(ANNOTATION)
    return path


@pytest.fixture
def feature_db(gff: Path) -> Iterator[gffutils.FeatureDB]:
    database = build_database(str(gff), None)
    yield database
    database.conn.close()


def classify(
    feature_db: gffutils.FeatureDB, types: list[str]
) -> tuple[FeatureHierarchy, ParentOrder]:
    return separate_parents_and_children(feature_db, types)


class TestBuildDatabase:
    def test_writes_database_next_to_annotation(
        self, feature_db: gffutils.FeatureDB, gff: Path
    ) -> None:
        assert Path(f'{gff}_db').is_file()

    def test_existing_database_can_be_reopened(
        self, feature_db: gffutils.FeatureDB, gff: Path
    ) -> None:
        reopened = build_database(None, f'{gff}_db')
        try:
            assert reopened.count_features_of_type('gene') == 2
        finally:
            reopened.conn.close()

    def test_percent_encoding_is_preserved(self, feature_db: gffutils.FeatureDB) -> None:
        assert feature_db['gene1'].attributes['Name'] == ['A%3Bescaped']

    @pytest.mark.parametrize(
        ('gff_name', 'db_name'),
        [
            pytest.param('missing.gff3', None, id='annotation'),
            pytest.param(None, 'missing_db', id='database'),
        ],
    )
    def test_missing_input_raises(
        self, tmp_path: Path, gff_name: str | None, db_name: str | None
    ) -> None:
        gff_path = str(tmp_path / gff_name) if gff_name else None
        db_path = str(tmp_path / db_name) if db_name else None
        with pytest.raises(InputError, match='not found'):
            build_database(gff_path, db_path)

    def test_syntax_error_reports_line_number(self, tmp_path: Path) -> None:
        bad = tmp_path / 'bad.gff3'
        bad.write_text(ANNOTATION + 'chr1\tsrc\tgene\tnot-a-number\t10\t.\t+\t.\tID=broken\n')
        with pytest.raises(GffSyntaxError) as excinfo:
            build_database(str(bad), None)
        assert excinfo.value.line_number == len(ANNOTATION.splitlines()) + 1


def test_annotation_format_is_detected(feature_db: gffutils.FeatureDB) -> None:
    assert annotation_format(feature_db) == 'gff3'


def test_feature_order_ranks_exon_then_cds_first(feature_db: gffutils.FeatureDB) -> None:
    order = get_feature_order(feature_db)
    assert (order['exon'], order['CDS']) == (0, 1)


class TestSeparateParentsAndChildren:
    def test_only_requested_types_are_parents(self, feature_db: gffutils.FeatureDB) -> None:
        hierarchy, _ = classify(feature_db, ['gene'])
        assert list(hierarchy.parents) == ['gene2', 'gene1']

    def test_parents_follow_annotation_order(self, feature_db: gffutils.FeatureDB) -> None:
        hierarchy, _ = classify(feature_db, ['gene', 'repeat_region'])
        assert list(hierarchy.parents) == ['gene2', 'gene1', 'repeat1']

    def test_parent_order_is_positional(self, feature_db: gffutils.FeatureDB) -> None:
        _, order = classify(feature_db, ['gene', 'repeat_region'])
        assert order.ids == ['gene1', 'repeat1', 'gene2']

    def test_children_are_leaf_features_without_introns(
        self, feature_db: gffutils.FeatureDB
    ) -> None:
        hierarchy, _ = classify(feature_db, ['gene'])
        assert [child.id for child in hierarchy.children['gene1']] == ['exon1', 'cds1', 'exon2']

    def test_children_keep_their_direct_parent(self, feature_db: gffutils.FeatureDB) -> None:
        hierarchy, _ = classify(feature_db, ['gene'])
        assert {child.attributes['Parent'][0] for child in hierarchy.children['gene1']} == {'tx1'}

    def test_transcripts_are_intermediates(self, feature_db: gffutils.FeatureDB) -> None:
        hierarchy, _ = classify(feature_db, ['gene'])
        assert set(hierarchy.intermediates) == {'tx1', 'tx2'}

    def test_single_level_feature_is_its_own_child(self, feature_db: gffutils.FeatureDB) -> None:
        hierarchy, _ = classify(feature_db, ['gene', 'repeat_region'])
        assert hierarchy.children['repeat1'] == [hierarchy.parents['repeat1']]
