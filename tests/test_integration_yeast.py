"""Whole-genome lift-over of S. cerevisiae onto itself (requires minimap2)."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import pytest

from liftoff.cli import main

pytestmark = [pytest.mark.minimap2, pytest.mark.slow]


def _run(args: list[str]) -> None:
    assert main(['-q', *args]) == 0


def test_yeast_basic(
    tmp_path: Path,
    data_dir: Path,
    yeast_inputs: dict[str, Path],
    same_annotation: Callable[[Path, Path], None],
) -> None:
    output = tmp_path / 'basic.gff3'
    unmapped = tmp_path / 'unmapped.txt'
    genome = str(yeast_inputs['genome'])
    _run(
        [
            '--gff',
            str(yeast_inputs['gff']),
            '--output',
            str(output),
            '--unmapped',
            str(unmapped),
            '--intermediate-dir',
            str(tmp_path / 'intermediate'),
            genome,
            genome,
        ]
    )
    same_annotation(output, data_dir / 'expected' / 'yeast_basic.gff3.gz')
    expected_unmapped = (data_dir / 'expected' / 'yeast_basic_unmapped.txt').read_text()
    assert unmapped.read_text() == expected_unmapped


def test_yeast_advanced(
    tmp_path: Path,
    data_dir: Path,
    yeast_inputs: dict[str, Path],
    same_annotation: Callable[[Path, Path], None],
) -> None:
    output = tmp_path / 'advanced.gff3'
    unmapped = tmp_path / 'unmapped.txt'
    genome = str(yeast_inputs['genome'])
    _run(
        [
            '--gff',
            str(yeast_inputs['gff']),
            '--output',
            str(output),
            '--unmapped',
            str(unmapped),
            '--intermediate-dir',
            str(tmp_path / 'intermediate'),
            '--chroms',
            str(yeast_inputs['chroms']),
            '--unplaced',
            str(yeast_inputs['unplaced']),
            '--copies',
            '--threads',
            '2',
            genome,
            genome,
        ]
    )
    same_annotation(output, data_dir / 'expected' / 'yeast_advanced.gff3.gz')
    expected_unmapped = (data_dir / 'expected' / 'yeast_advanced_unmapped.txt').read_text()
    assert unmapped.read_text() == expected_unmapped
