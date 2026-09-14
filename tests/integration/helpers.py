"""Shared types and helpers for integration tests."""

from __future__ import annotations

from dataclasses import dataclass
import gzip
from pathlib import Path

#: SAM file names requested by the whole-genome stages of a chromosome I run.
CHR1_SAM_NAMES = ('reference_all_to_target_all.sam', 'reference_all_copies_to_target_all.sam')


@dataclass(frozen=True)
class Chr1Inputs:
    """Scratch copies of the chromosome I test inputs."""

    reference: Path
    target: Path
    gff: Path


@dataclass(frozen=True)
class RunPaths:
    """Output locations of one Liftoff run rooted at ``workdir``."""

    workdir: Path

    @property
    def output(self) -> Path:
        return self.workdir / 'out.gff3'

    @property
    def polished(self) -> Path:
        return self.workdir / 'out.gff3_polished'

    @property
    def unmapped(self) -> Path:
        return self.workdir / 'unmapped.txt'

    @property
    def intermediate(self) -> Path:
        return self.workdir / 'intermediate'


def read_feature_lines(path: Path) -> list[str]:
    """Return annotation lines, excluding ``#`` header/comment lines."""
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as handle:
        return [line.rstrip('\n') for line in handle if not line.startswith('#')]
