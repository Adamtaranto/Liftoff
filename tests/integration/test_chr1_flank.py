"""Lifting chromosome I onto itself with flanking sequence."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import pytest

from .helpers import Chr1Inputs, RunPaths

pytestmark = pytest.mark.minimap2

#: Top-level feature types lifted in this test; all but ``gene`` are
#: single-level features (their own only child).
TOP_LEVEL_TYPES = ('gene', 'origin_of_replication', 'telomere', 'long_terminal_repeat')


def top_level_coordinates(path: Path) -> dict[str, tuple[int, int]]:
    """Map feature ID to ``(start, end)`` for original (non-copy) top-level features."""
    coordinates = {}
    for line in path.read_text().splitlines():
        fields = line.split('\t')
        if line.startswith('#') or len(fields) < 9 or fields[2] not in TOP_LEVEL_TYPES:
            continue
        attributes = dict(item.split('=', 1) for item in fields[8].split(';'))
        if attributes.get('extra_copy_number', '0') == '0':
            coordinates[attributes['ID']] = (int(fields[3]), int(fields[4]))
    return coordinates


@pytest.fixture(scope='module')
def flank_run(
    tmp_path_factory: pytest.TempPathFactory,
    copy_chr1_inputs: Callable[[Path], Chr1Inputs],
    run_cli: Callable[..., RunPaths],
) -> tuple[Chr1Inputs, RunPaths]:
    """Self lift-over with 20% flank, extra copies and single-level feature types."""
    workdir = tmp_path_factory.mktemp('chr1_flank')
    inputs = copy_chr1_inputs(workdir / 'inputs')
    feature_types = workdir / 'feature_types.txt'
    feature_types.write_text('\n'.join(TOP_LEVEL_TYPES[1:]) + '\n')
    extra = ['--flank', '0.2', '--copies', '--feature-types', str(feature_types)]
    paths = run_cli(workdir, inputs.reference, inputs.reference, inputs.gff, extra)
    return inputs, paths


def test_all_top_level_features_are_lifted(flank_run: tuple[Chr1Inputs, RunPaths]) -> None:
    inputs, paths = flank_run
    assert top_level_coordinates(paths.output).keys() == top_level_coordinates(inputs.gff).keys()


def test_lifted_coordinates_exclude_the_flank(flank_run: tuple[Chr1Inputs, RunPaths]) -> None:
    inputs, paths = flank_run
    assert top_level_coordinates(paths.output) == top_level_coordinates(inputs.gff)
