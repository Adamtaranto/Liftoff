"""Lifting every top-level feature type of chromosome I with ``--all-feature-types``."""

from __future__ import annotations

from collections import Counter
from collections.abc import Callable
from pathlib import Path

import pytest

from .helpers import Chr1Inputs, RunPaths, read_feature_lines

pytestmark = pytest.mark.minimap2

#: Features shorter than this cannot be aligned by minimap2 and stay unmapped.
MIN_ALIGNABLE_LENGTH = 30


def top_level_features(lines: list[str]) -> dict[str, tuple[str, int]]:
    """Map top-level (Parent-less) feature IDs to ``(type, length)``."""
    features = {}
    for line in lines:
        fields = line.split('\t')
        attributes = dict(item.split('=', 1) for item in fields[8].split(';') if '=' in item)
        if 'Parent' not in attributes and attributes.get('extra_copy_number', '0') == '0':
            features[attributes['ID']] = (fields[2], int(fields[4]) - int(fields[3]) + 1)
    return features


def gene_lines(lines: list[str]) -> list[str]:
    """Return lines belonging to gene records (each gene followed by its descendants)."""
    selected, in_gene = [], False
    for line in lines:
        fields = line.split('\t')
        if 'Parent=' not in fields[8]:
            in_gene = fields[2] == 'gene'
        if in_gene:
            selected.append(line)
    return selected


@pytest.fixture(scope='module')
def all_types_run(
    tmp_path_factory: pytest.TempPathFactory,
    copy_chr1_inputs: Callable[[Path], Chr1Inputs],
    run_cli: Callable[..., RunPaths],
) -> tuple[Chr1Inputs, RunPaths]:
    workdir = tmp_path_factory.mktemp('chr1_all_feature_types')
    inputs = copy_chr1_inputs(workdir / 'inputs')
    paths = run_cli(workdir, inputs.target, inputs.reference, inputs.gff, ['--all-feature-types'])
    return inputs, paths


def test_every_alignable_top_level_feature_is_lifted(
    all_types_run: tuple[Chr1Inputs, RunPaths],
) -> None:
    inputs, paths = all_types_run
    reference = top_level_features(read_feature_lines(inputs.gff))
    lifted = top_level_features(read_feature_lines(paths.output))
    alignable = {k for k, (_, length) in reference.items() if length >= MIN_ALIGNABLE_LENGTH}
    assert alignable <= lifted.keys()


def test_all_reference_feature_types_are_present(
    all_types_run: tuple[Chr1Inputs, RunPaths],
) -> None:
    inputs, paths = all_types_run
    reference_types = Counter(
        t for t, _ in top_level_features(read_feature_lines(inputs.gff)).values()
    )
    lifted_types = Counter(
        t for t, _ in top_level_features(read_feature_lines(paths.output)).values()
    )
    assert set(lifted_types) == set(reference_types)


def test_only_unalignable_features_are_unmapped(all_types_run: tuple[Chr1Inputs, RunPaths]) -> None:
    inputs, paths = all_types_run
    reference = top_level_features(read_feature_lines(inputs.gff))
    unmapped = paths.unmapped.read_text().split()
    assert all(reference[feature_id][1] < MIN_ALIGNABLE_LENGTH for feature_id in unmapped)


def test_genes_are_lifted_as_without_the_option(
    all_types_run: tuple[Chr1Inputs, RunPaths], expected_dir: Path
) -> None:
    _, paths = all_types_run
    expected = gene_lines(read_feature_lines(expected_dir / 'chr1_basic.gff3'))
    assert gene_lines(read_feature_lines(paths.output)) == expected
