"""Alignment backends that obtain SAM files from outside Liftoff.

These backends make it possible to run the full pipeline where the minimap2
executable is unavailable, such as in Pyodide/WebAssembly.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from liftoff.align.base import Aligner, AlignmentJob
from liftoff.errors import AlignmentError


@dataclass(slots=True)
class PrecomputedSamAligner(Aligner):
    """Look up SAM files produced ahead of time.

    For each job the file ``<directory>/<job.output_sam.name>`` is used. Run
    Liftoff once with ``--verbose`` (or read the error message raised for a
    missing file) to discover which alignments are required; each lists
    the query FASTA written to the intermediate directory, the target, and the
    minimap2 options to use.

    Attributes
    ----------
    directory : pathlib.Path
        Directory containing the SAM files.
    """

    directory: Path
    parallel_safe = True

    def align(self, job: AlignmentJob) -> Path:
        """Return the pre-computed SAM file for ``job``.

        Parameters
        ----------
        job : AlignmentJob
            The alignment being requested.

        Returns
        -------
        pathlib.Path
            Path of the matching SAM file.

        Raises
        ------
        AlignmentError
            If no SAM file with the expected name exists.
        """
        candidate = Path(self.directory) / job.output_sam.name
        if not candidate.is_file():
            raise AlignmentError(
                f'pre-computed alignment {candidate} not found. Produce it with:\n'
                f'  minimap2 {" ".join(job.minimap2_options)} {job.target_fasta} '
                f'{job.features_fasta} > {candidate}'
            )
        return candidate


@dataclass(slots=True)
class CallableAligner(Aligner):
    """Delegate alignment to a Python callable.

    Attributes
    ----------
    function : callable
        Called with the :class:`~liftoff.align.base.AlignmentJob`. It must
        write SAM output to ``job.output_sam`` and may return ``None``, or
        return the path of the SAM file it wrote.
    """

    function: Callable[[AlignmentJob], str | Path | None]

    def align(self, job: AlignmentJob) -> Path:
        """Invoke the wrapped callable.

        Parameters
        ----------
        job : AlignmentJob
            The alignment to perform.

        Returns
        -------
        pathlib.Path
            Path of the SAM file produced.

        Raises
        ------
        AlignmentError
            If the callable did not produce a SAM file.
        """
        result = self.function(job)
        path = job.output_sam if result is None else Path(result)
        if not path.is_file():
            raise AlignmentError(f'aligner callable did not produce a SAM file at {path}')
        return path
