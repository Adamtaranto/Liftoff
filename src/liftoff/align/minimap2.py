"""Alignment backend that runs the minimap2 executable."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import shutil
import subprocess

from liftoff.align.base import Aligner, AlignmentJob
from liftoff.errors import AlignmentError
from liftoff.log import get_logger

logger = get_logger(__name__)


@dataclass(slots=True)
class Minimap2Aligner(Aligner):
    """Run ``minimap2`` in a subprocess.

    Attributes
    ----------
    executable : str, default "minimap2"
        Path to, or name of, the minimap2 executable.
    """

    executable: str = 'minimap2'
    parallel_safe = True

    def resolve_executable(self) -> str:
        """Return the absolute path of the executable.

        Returns
        -------
        str
            Resolved executable path.

        Raises
        ------
        AlignmentError
            If minimap2 cannot be found.
        """
        resolved = shutil.which(self.executable)
        if resolved is None:
            raise AlignmentError(
                f'minimap2 executable {self.executable!r} was not found. Install minimap2, '
                'pass its location with --minimap2, or supply pre-computed alignments '
                'with --alignments.'
            )
        return resolved

    def align(self, job: AlignmentJob) -> Path:
        """Align ``job.features_fasta`` to ``job.target_fasta``.

        Small targets are indexed once (``<target>.mmi``, reused across
        stages); targets above minimap2's single-index size are aligned with
        ``--split-prefix`` instead.

        Parameters
        ----------
        job : AlignmentJob
            The alignment to perform.

        Returns
        -------
        pathlib.Path
            The SAM file written to ``job.output_sam``.
        """
        executable = self.resolve_executable()
        options = list(job.minimap2_options)
        threads = ['-t', str(job.threads)]
        if job.split_prefix is not None:
            target = str(job.target_fasta)
            extra = ['--split-prefix', str(job.split_prefix), *threads]
        else:
            target = str(self._build_index(executable, job))
            extra = threads
        command = [
            executable,
            '-o',
            str(job.output_sam),
            target,
            str(job.features_fasta),
            *options,
            *extra,
        ]
        _run(command)
        return job.output_sam

    @staticmethod
    def _build_index(executable: str, job: AlignmentJob) -> Path:
        """Build ``<target>.mmi`` unless it already exists."""
        index = Path(f'{job.target_fasta}.mmi')
        if not index.exists():
            _run(
                [
                    executable,
                    '-d',
                    str(index),
                    str(job.target_fasta),
                    *job.minimap2_options,
                    '-t',
                    str(job.threads),
                ]
            )
        return index


def _run(command: list[str]) -> None:
    """Run a minimap2 command, logging its stderr at debug level."""
    logger.debug('running: %s', ' '.join(command))
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    for line in completed.stderr.splitlines():
        logger.debug('minimap2: %s', line)
    if completed.returncode != 0:
        tail = '\n'.join(completed.stderr.splitlines()[-20:])
        raise AlignmentError(
            f'minimap2 exited with status {completed.returncode}: {" ".join(command)}\n{tail}'
        )
