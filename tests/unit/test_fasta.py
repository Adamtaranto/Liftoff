"""Tests for FASTA helpers used to extract gene and chromosome sequences."""

from __future__ import annotations

from pathlib import Path

import pytest

from liftoff.errors import InputError
from liftoff.io.fasta import (
    features_file_name,
    get_genome_size,
    open_fasta,
    split_target_sequence,
    write_gene_sequences,
)
from liftoff.models import LiftoverType

from ..helpers import FeatureFactory

SEQUENCE = 'ACGT' * 100  # 400 bases


@pytest.fixture
def genome(tmp_path: Path) -> Path:
    path = tmp_path / 'genome.fa'
    path.write_text(f'>chr1 description\n{SEQUENCE}\n>chr2\n{SEQUENCE[:40]}\n')
    return path


class TestFeaturesFileName:
    @pytest.mark.parametrize(
        ('stage', 'expected'),
        [
            pytest.param(LiftoverType.CHROM_BY_CHROM, 'reference_all', id='chrom-by-chrom'),
            pytest.param(LiftoverType.COPIES, 'reference_all_copies', id='copies'),
            pytest.param(LiftoverType.UNMAPPED, 'unmapped_to_expected_chrom', id='unmapped'),
            pytest.param(LiftoverType.UNPLACED, 'unplaced', id='unplaced'),
        ],
    )
    def test_whole_genome_names_are_unique_per_stage(
        self, stage: LiftoverType, expected: str
    ) -> None:
        assert features_file_name('ref.fa', 'ref.fa', stage) == expected

    def test_chromosome_stage_uses_chromosome_name(self) -> None:
        assert features_file_name('chrX', 'ref.fa', LiftoverType.CHROM_BY_CHROM) == 'chrX'


def test_open_fasta_missing_file(tmp_path: Path) -> None:
    with pytest.raises(InputError, match='not found'):
        open_fasta(str(tmp_path / 'missing.fa'))


def test_genome_size_sums_all_sequences(genome: Path) -> None:
    assert get_genome_size(open_fasta(str(genome))) == 440


def test_split_target_writes_requested_chromosomes(genome: Path, tmp_path: Path) -> None:
    split_target_sequence(['chr2'], str(genome), tmp_path)
    assert (tmp_path / 'chr2.fa').read_text() == f'>chr2\n{SEQUENCE[:40]}'


class TestWriteGeneSequences:
    def test_writes_feature_sequence(
        self, genome: Path, tmp_path: Path, make_feature: FeatureFactory
    ) -> None:
        gene = make_feature('g', 'gene', 5, 12)
        write_gene_sequences(
            {'g': gene}, [str(genome)], str(genome), tmp_path, LiftoverType.CHROM_BY_CHROM, 0
        )
        assert (tmp_path / 'reference_all_genes.fa').read_text() == f'>g\n{SEQUENCE[4:12]}\n'

    def test_flank_extends_extracted_region(
        self, genome: Path, tmp_path: Path, make_feature: FeatureFactory
    ) -> None:
        gene = make_feature('g', 'gene', 101, 200)
        write_gene_sequences(
            {'g': gene}, [str(genome)], str(genome), tmp_path, LiftoverType.CHROM_BY_CHROM, 0.1
        )
        assert (gene.start, gene.end) == (91, 210)

    def test_flank_is_clipped_to_the_chromosome(
        self, genome: Path, tmp_path: Path, make_feature: FeatureFactory
    ) -> None:
        gene = make_feature('g', 'gene', 1, 400)
        write_gene_sequences(
            {'g': gene}, [str(genome)], str(genome), tmp_path, LiftoverType.CHROM_BY_CHROM, 0.5
        )
        assert (gene.start, gene.end) == (1, 400)

    def test_flank_is_not_applied_cumulatively(
        self, genome: Path, tmp_path: Path, make_feature: FeatureFactory
    ) -> None:
        gene = make_feature('g', 'gene', 101, 200)
        for stage in (LiftoverType.CHROM_BY_CHROM, LiftoverType.UNMAPPED, LiftoverType.COPIES):
            write_gene_sequences({'g': gene}, [str(genome)], str(genome), tmp_path, stage, 0.1)
        assert (gene.start, gene.end) == (91, 210)

    def test_annotated_bounds_are_preserved(
        self, genome: Path, tmp_path: Path, make_feature: FeatureFactory
    ) -> None:
        gene = make_feature('g', 'gene', 101, 200)
        write_gene_sequences(
            {'g': gene}, [str(genome)], str(genome), tmp_path, LiftoverType.CHROM_BY_CHROM, 0.1
        )
        assert gene.annotated_bounds == (101, 200)

    def test_rejects_empty_feature_set(self, genome: Path, tmp_path: Path) -> None:
        with pytest.raises(InputError, match='does not contain any gene features'):
            write_gene_sequences(
                {}, [str(genome)], str(genome), tmp_path, LiftoverType.CHROM_BY_CHROM, 0
            )
