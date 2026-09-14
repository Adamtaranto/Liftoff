"""Shared configuration and fixtures for the Liftoff test-suite.

Layout
------
``tests/unit/``
    Fast, isolated tests of individual modules (marked ``unit``).
``tests/integration/``
    End-to-end lift-overs on real data (marked ``integration``).
``tests/data/``
    Input genomes, expected outputs and parasail corpora; see its README.

Markers ``unit`` and ``integration`` are applied automatically from the
directory. Tests marked ``minimap2`` are skipped when the executable is not on
``PATH`` (e.g. under Pyodide).
"""

from __future__ import annotations

from collections.abc import Iterator
import logging
from pathlib import Path
import random
import shutil

import pytest

from liftoff.log import ROOT_LOGGER_NAME
from liftoff.models import Feature

from .helpers import FeatureFactory, SamLineFactory

DATA_DIR = Path(__file__).parent / 'data'


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    """Apply directory-based markers and skip minimap2 tests without minimap2."""
    minimap2_missing = shutil.which('minimap2') is None
    skip_minimap2 = pytest.mark.skip(reason='minimap2 executable not found on PATH')
    for item in items:
        parts = item.path.parts
        if 'integration' in parts:
            item.add_marker(pytest.mark.integration)
        elif 'unit' in parts:
            item.add_marker(pytest.mark.unit)
        if minimap2_missing and item.get_closest_marker('minimap2') is not None:
            item.add_marker(skip_minimap2)


@pytest.fixture(autouse=True)
def restore_package_logger() -> Iterator[logging.Logger]:
    """Undo logging configuration made by a test (e.g. via ``liftoff.cli.main``)."""
    logger = logging.getLogger(ROOT_LOGGER_NAME)
    handlers, level, propagate = list(logger.handlers), logger.level, logger.propagate
    yield logger
    logger.handlers[:] = handlers
    logger.setLevel(level)
    logger.propagate = propagate


# ---------------------------------------------------------------------------
# Test data locations
# ---------------------------------------------------------------------------
@pytest.fixture(scope='session')
def data_dir() -> Path:
    """Root of the test data directory."""
    return DATA_DIR


@pytest.fixture(scope='session')
def expected_dir(data_dir: Path) -> Path:
    """Directory of expected lift-over outputs."""
    return data_dir / 'expected'


@pytest.fixture(scope='session')
def parasail_dir(data_dir: Path) -> Path:
    """Directory of parasail reference corpora and traces."""
    return data_dir / 'parasail'


# ---------------------------------------------------------------------------
# Object factories
# ---------------------------------------------------------------------------
@pytest.fixture
def make_feature() -> FeatureFactory:
    """Return a factory building :class:`~liftoff.models.Feature` objects.

    Only the fields relevant to a test need to be given; ``ID`` (and
    ``Parent`` when ``parent`` is set) attributes are filled in automatically.
    """

    def factory(
        feature_id: str = 'gene1',
        featuretype: str = 'gene',
        start: int = 1,
        end: int = 100,
        *,
        strand: str = '+',
        seqid: str = 'chr1',
        frame: str = '.',
        parent: str | None = None,
        source: str = 'test',
        **attributes: list[str],
    ) -> Feature:
        attrs: dict[str, list[str]] = {'ID': [feature_id]}
        if parent is not None:
            attrs['Parent'] = [parent]
        attrs.update(attributes)
        return Feature(feature_id, featuretype, seqid, source, strand, start, end, frame, attrs)

    return factory


@pytest.fixture
def make_sam_line() -> SamLineFactory:
    """Return a factory building tab-separated SAM alignment lines."""

    def factory(
        qname: str = 'gene1',
        flag: int = 0,
        rname: str = 'chr1',
        pos: int = 1,
        cigar: str = '10=',
        seq: str = '*',
    ) -> str:
        fields = [qname, str(flag), rname, str(pos), '60', cigar, '*', '0', '0', seq, '*']
        return '\t'.join(fields)

    return factory


@pytest.fixture
def rng() -> random.Random:
    """A deterministically seeded random number generator."""
    return random.Random(20260915)
