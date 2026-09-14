"""Shared fixtures for the Liftoff test-suite.

Test data layout (see ``tests/data/README.md``):

* ``yeast/``    - S. cerevisiae R64 genome and annotation (whole-genome tests)
* ``chr1/``     - chromosome I subset and a CDS-mutated copy (fast tests)
* ``expected/`` - reference outputs recorded with Liftoff 1.6.3 (parasail/pysam)
* ``parasail/`` - parasail alignment corpora and polishing traces
"""

from __future__ import annotations

from collections.abc import Callable
import gzip
from pathlib import Path
import shutil

import pytest

DATA_DIR = Path(__file__).parent / 'data'


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    """Skip tests marked ``minimap2`` when the executable is unavailable."""
    if shutil.which('minimap2') is not None:
        return
    skip = pytest.mark.skip(reason='minimap2 executable not found on PATH')
    for item in items:
        if 'minimap2' in item.keywords:
            item.add_marker(skip)


@pytest.fixture(scope='session')
def data_dir() -> Path:
    """Root of the test data directory."""
    return DATA_DIR


@pytest.fixture
def chr1_inputs(tmp_path: Path) -> dict[str, Path]:
    """Copy the chromosome I inputs into a scratch directory.

    Liftoff writes indices (``.fai``, ``.mmi``, ``_db``) next to its inputs, so
    tests never run directly against files in ``tests/data``.
    """
    work = tmp_path / 'inputs'
    work.mkdir()
    paths = {}
    for key, name in [
        ('reference', 'ref.fa'),
        ('target', 'target_mutated.fa'),
        ('gff', 'annotation.gff3'),
    ]:
        paths[key] = Path(shutil.copy(DATA_DIR / 'chr1' / name, work / name))
    return paths


@pytest.fixture
def yeast_inputs(tmp_path: Path) -> dict[str, Path]:
    """Copy the whole-genome yeast inputs into a scratch directory."""
    work = tmp_path / 'inputs'
    shutil.copytree(DATA_DIR / 'yeast', work)
    return {
        'genome': work / 'GCA_000146045.2_R64_genomic.fna.gz',
        'gff': work / 'GCA_000146045.2_R64_genomic.gff.gz',
        'chroms': work / 'chroms.txt',
        'unplaced': work / 'unplaced.txt',
    }


def _read_features(path: Path) -> list[str]:
    opener: Callable[..., object] = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as handle:  # type: ignore[operator]
        return [line.rstrip('\n') for line in handle if not line.startswith('#')]


@pytest.fixture(scope='session')
def read_features() -> Callable[[Path], list[str]]:
    """Return a reader yielding annotation lines without header comments."""
    return _read_features


def assert_same_annotation(observed: Path, expected: Path) -> None:
    """Assert two annotation files are identical, ignoring ``#`` header lines."""
    observed_lines = _read_features(observed)
    expected_lines = _read_features(expected)
    assert len(observed_lines) == len(expected_lines), (
        f'{observed} has {len(observed_lines)} feature lines, expected {len(expected_lines)}'
    )
    for number, (got, want) in enumerate(zip(observed_lines, expected_lines, strict=True), 1):
        assert got == want, f'line {number} differs:\n got: {got}\nwant: {want}'


@pytest.fixture(scope='session')
def same_annotation() -> Callable[[Path, Path], None]:
    """Return the annotation comparison helper."""
    return assert_same_annotation
