"""Whole-genome lift-over of S. cerevisiae R64 onto itself.

Expected outputs were recorded with the original Liftoff 1.6.3.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
import shutil

import pytest

from .helpers import RunPaths

pytestmark = [pytest.mark.minimap2, pytest.mark.slow]

GENOME = 'GCA_000146045.2_R64_genomic.fna.gz'
ANNOTATION = 'GCA_000146045.2_R64_genomic.gff.gz'


@dataclass(frozen=True)
class YeastRun:
    scenario: str
    paths: RunPaths


def _scenario_options(scenario: str, inputs: Path) -> list[str]:
    if scenario == 'basic':
        return []
    return [
        '--chroms',
        str(inputs / 'chroms.txt'),
        '--unplaced',
        str(inputs / 'unplaced.txt'),
        '--copies',
        '--threads',
        '2',
    ]


@pytest.fixture(scope='module', params=['basic', 'advanced'])
def yeast_run(
    request: pytest.FixtureRequest,
    tmp_path_factory: pytest.TempPathFactory,
    data_dir: Path,
    run_cli: Callable[..., RunPaths],
) -> YeastRun:
    """Run one whole-genome scenario (``advanced`` adds chroms, unplaced and copies)."""
    scenario = str(request.param)
    workdir = tmp_path_factory.mktemp(f'yeast_{scenario}')
    inputs = Path(shutil.copytree(data_dir / 'yeast', workdir / 'inputs'))
    genome = inputs / GENOME
    extra = _scenario_options(scenario, inputs)
    paths = run_cli(workdir, genome, genome, inputs / ANNOTATION, extra)
    return YeastRun(scenario, paths)


def test_annotation_matches_expected(
    yeast_run: YeastRun,
    expected_dir: Path,
    assert_same_annotation: Callable[[Path, Path], None],
) -> None:
    expected = expected_dir / f'yeast_{yeast_run.scenario}.gff3.gz'
    assert_same_annotation(yeast_run.paths.output, expected)


def test_unmapped_features_match_expected(yeast_run: YeastRun, expected_dir: Path) -> None:
    expected = expected_dir / f'yeast_{yeast_run.scenario}_unmapped.txt'
    assert yeast_run.paths.unmapped.read_text() == expected.read_text()
