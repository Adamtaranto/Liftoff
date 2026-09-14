"""Fixtures for end-to-end lift-over tests.

Pipeline runs are comparatively slow, so each scenario is run once by a
module-scoped fixture and inspected by several focused tests. Inputs are
always copied to a temporary directory because Liftoff writes indices
(``.fai``, ``.mmi``, ``_db``) next to its input files.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from pathlib import Path
import shutil

import pytest

from liftoff.cli import main

from .helpers import CHR1_SAM_NAMES, Chr1Inputs, RunPaths, read_feature_lines


@pytest.fixture(scope='session')
def assert_same_annotation() -> Callable[[Path, Path], None]:
    """Return an assertion comparing two annotations, ignoring header lines.

    Reports the first differing line rather than a diff of the whole file,
    which can be tens of thousands of lines long.
    """

    def check(observed: Path, expected: Path) -> None:
        observed_lines = read_feature_lines(observed)
        expected_lines = read_feature_lines(expected)
        for number, (got, want) in enumerate(
            zip(observed_lines, expected_lines, strict=False), start=1
        ):
            assert got == want, f'line {number} of {observed.name} differs from {expected.name}'
        assert len(observed_lines) == len(expected_lines), (
            f'{observed.name} has {len(observed_lines)} feature lines, '
            f'{expected.name} has {len(expected_lines)}'
        )

    return check


@pytest.fixture(scope='session')
def copy_chr1_inputs(data_dir: Path) -> Callable[[Path], Chr1Inputs]:
    """Return a function copying the chromosome I inputs into a directory."""

    def copy(destination: Path) -> Chr1Inputs:
        destination.mkdir(parents=True, exist_ok=True)
        source = data_dir / 'chr1'
        return Chr1Inputs(
            reference=Path(shutil.copy(source / 'ref.fa', destination)),
            target=Path(shutil.copy(source / 'target_mutated.fa', destination)),
            gff=Path(shutil.copy(source / 'annotation.gff3', destination)),
        )

    return copy


@pytest.fixture(scope='session')
def chr1_alignments(tmp_path_factory: pytest.TempPathFactory, expected_dir: Path) -> Path:
    """Directory of pre-computed SAM files for chromosome I runs.

    The copies stage aligns the same gene sequences to the same target with the
    same options, so its alignment is identical to the first stage's.
    """
    directory = tmp_path_factory.mktemp('chr1_alignments')
    recorded = expected_dir / f'chr1_basic_{CHR1_SAM_NAMES[0]}'
    for name in CHR1_SAM_NAMES:
        shutil.copy(recorded, directory / name)
    return directory


@pytest.fixture(scope='session')
def run_cli() -> Callable[..., RunPaths]:
    """Return a function running the ``liftoff`` CLI with standard output paths."""

    def run(
        workdir: Path,
        target: Path,
        reference: Path,
        gff: Path,
        extra: Sequence[str] = (),
    ) -> RunPaths:
        paths = RunPaths(workdir)
        argv = [
            '--quiet',
            '--gff',
            str(gff),
            '--output',
            str(paths.output),
            '--unmapped',
            str(paths.unmapped),
            '--intermediate-dir',
            str(paths.intermediate),
            *extra,
            str(target),
            str(reference),
        ]
        exit_status = main(argv)
        assert exit_status == 0, f'liftoff exited with status {exit_status}'
        return paths

    return run
