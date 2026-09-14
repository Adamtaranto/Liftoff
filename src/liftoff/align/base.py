"""Pluggable alignment backends.

Liftoff needs gene sequences aligned to the target genome in SAM format. The
pipeline describes each required alignment as an :class:`AlignmentJob` and
hands it to an :class:`Aligner`, which returns the path of the resulting SAM
file. Backends:

* :class:`~liftoff.align.minimap2.Minimap2Aligner` runs the ``minimap2``
  executable (the default for native installs).
* :class:`~liftoff.align.external.PrecomputedSamAligner` reads SAM files that
  were produced elsewhere.
* :class:`~liftoff.align.external.CallableAligner` delegates to a Python
  function — for example one that calls a WebAssembly build of minimap2 from a
  Pyodide host.
"""

from __future__ import annotations

import abc
from dataclasses import dataclass, field
from pathlib import Path

from liftoff.models import LiftoverType


@dataclass(frozen=True, slots=True)
class AlignmentJob:
    """Description of one gene-to-genome alignment.

    Attributes
    ----------
    features_fasta : pathlib.Path
        FASTA of extracted reference gene sequences (the queries).
    target_fasta : pathlib.Path
        Target genome, or a single target chromosome, to align to.
    output_sam : pathlib.Path
        Where the SAM output is expected. Its file name is unique within one
        pipeline stage and is used to look up pre-computed alignments.
    minimap2_options : tuple of str
        Arguments passed to minimap2 (always includes ``-a --eqx``).
    threads : int
        Threads available to this alignment.
    liftover_type : LiftoverType
        Pipeline stage requesting the alignment.
    split_prefix : pathlib.Path or None, default None
        When set, the target is too large for a single minimap2 index and
        ``--split-prefix`` should be used instead of a pre-built index.
    """

    features_fasta: Path
    target_fasta: Path
    output_sam: Path
    minimap2_options: tuple[str, ...]
    threads: int
    liftover_type: LiftoverType
    split_prefix: Path | None = field(default=None)


class Aligner(abc.ABC):
    """Base class for alignment backends."""

    #: Whether jobs may be dispatched to worker processes. Backends holding
    #: unpicklable state (such as arbitrary callables) must leave this False.
    parallel_safe: bool = False

    @abc.abstractmethod
    def align(self, job: AlignmentJob) -> Path:
        """Produce the SAM file for ``job``.

        Parameters
        ----------
        job : AlignmentJob
            The alignment to perform.

        Returns
        -------
        pathlib.Path
            Path to a SAM file with one record per alignment (unmapped queries
            may be reported with FLAG 4).
        """
