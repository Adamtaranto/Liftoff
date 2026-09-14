"""Liftoff: accurate mapping of gene annotations between genome assemblies.

The command line interface lives in :mod:`liftoff.cli`. For programmatic use,
build a :class:`~liftoff.config.LiftoffConfig` and call
:func:`~liftoff.pipeline.run_liftoff`, optionally supplying a custom
:class:`~liftoff.align.base.Aligner` (for example when running under Pyodide,
where the minimap2 executable is unavailable).
"""

from __future__ import annotations

try:
    # Written by hatch-vcs at build/install time.
    from liftoff._version import __version__
except ImportError:  # pragma: no cover - only hit when running from a raw checkout
    __version__ = '0.0.0+unknown'

__all__ = ['__version__']
