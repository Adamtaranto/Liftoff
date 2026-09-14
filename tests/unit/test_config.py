"""Tests for run configuration validation and minimap2 option handling."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

from liftoff.config import DEFAULT_MM2_OPTIONS, LiftoffConfig, merge_minimap2_options
from liftoff.errors import ConfigError


def make_config(**overrides: Any) -> LiftoffConfig:
    options: dict[str, Any] = {'target': 't.fa', 'reference': 'r.fa', 'gff': 'a.gff3'}
    options.update(overrides)
    return LiftoffConfig(**options)


def test_defaults_are_valid() -> None:
    config = make_config()
    assert (config.output, config.cds, config.threads) == ('stdout', True, 1)


@pytest.mark.parametrize(
    ('overrides', 'message'),
    [
        pytest.param({'gff': None}, 'exactly one', id='no-annotation'),
        pytest.param({'db': 'a.db'}, 'exactly one', id='annotation-and-database'),
        pytest.param(
            {'min_coverage': 1.5}, 'min_coverage must be between 0 and 1', id='coverage-range'
        ),
        pytest.param(
            {'min_identity': -0.1}, 'min_identity must be between 0 and 1', id='identity-range'
        ),
        pytest.param({'flank': 2.0}, 'flank must be between 0 and 1', id='flank-range'),
        pytest.param(
            {'max_overlap': 1.1}, 'max_overlap must be between 0 and 1', id='overlap-range'
        ),
        pytest.param(
            {'min_identity': 0.9, 'copy_identity': 0.5}, 'copy identity', id='copy-identity'
        ),
        pytest.param({'unplaced': 'u.txt'}, 'chroms', id='unplaced-without-chroms'),
        pytest.param({'threads': 0}, 'threads', id='threads'),
        pytest.param({'distance_factor': 0}, 'distance factor', id='distance-factor'),
        pytest.param({'gap_open': -1}, 'gap_open penalty', id='negative-penalty'),
    ],
)
def test_invalid_options_raise_config_error(overrides: dict[str, Any], message: str) -> None:
    with pytest.raises(ConfigError, match=message):
        make_config(**overrides)


def test_intermediate_path_is_a_path() -> None:
    assert make_config(intermediate_dir='work/tmp').intermediate_path == Path('work/tmp')


class TestMergeMinimap2Options:
    def test_default_options_are_unchanged(self) -> None:
        assert merge_minimap2_options(DEFAULT_MM2_OPTIONS) == DEFAULT_MM2_OPTIONS.split()

    def test_required_options_are_appended(self) -> None:
        assert merge_minimap2_options('-r 2k') == [
            '-r',
            '2k',
            '-a',
            '--end-bonus',
            '5',
            '--eqx',
            '-N',
            '50',
            '-p',
            '0.5',
        ]

    def test_user_values_take_precedence(self) -> None:
        tokens = merge_minimap2_options('-N 5 --end-bonus 10')
        assert tokens[:4] == ['-N', '5', '--end-bonus', '10']

    def test_flags_are_matched_as_whole_tokens(self) -> None:
        # '--end-bonus' contains '-N'-like substrings but must not satisfy '-N'.
        assert '-N' in merge_minimap2_options('--end-bonus 10')

    def test_shell_quoting_is_honoured(self) -> None:
        assert merge_minimap2_options("--rg 'a b'")[:2] == ['--rg', 'a b']

    def test_config_exposes_merged_arguments(self) -> None:
        assert make_config(mm2_options='-r 2k').minimap2_arguments[:2] == ['-r', '2k']
