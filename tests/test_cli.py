"""Tests for command line parsing and configuration validation."""

from __future__ import annotations

from pathlib import Path

import pytest

from liftoff import __version__
from liftoff.cli import build_parser, config_from_args, main
from liftoff.config import DEFAULT_MM2_OPTIONS, LiftoffConfig, merge_minimap2_options
from liftoff.errors import ConfigError


def _parse(*argv: str) -> LiftoffConfig:
    args = build_parser().parse_args(argv)
    return config_from_args(args, argv)


def test_defaults() -> None:
    config = _parse('-g', 'a.gff', 'target.fa', 'ref.fa')
    assert config.target == 'target.fa'
    assert config.reference == 'ref.fa'
    assert config.output == 'stdout'
    assert config.cds is True
    assert config.minimap2_arguments == DEFAULT_MM2_OPTIONS.split()
    assert config.command_line == 'liftoff -g a.gff target.fa ref.fa'


def test_long_options() -> None:
    config = _parse(
        '--db',
        'a.db',
        '--min-coverage',
        '0.9',
        '--min-identity',
        '0.8',
        '--copy-identity',
        '0.95',
        '--copies',
        '--no-cds',
        '--threads',
        '4',
        '--mm2-options=-r 2k -z 5000',
        't.fa',
        'r.fa',
    )
    assert (config.db, config.gff) == ('a.db', None)
    assert (config.min_coverage, config.min_identity, config.copy_identity) == (0.9, 0.8, 0.95)
    assert config.copies
    assert config.cds is False
    assert config.threads == 4
    assert config.minimap2_arguments[:4] == ['-r', '2k', '-z', '5000']


def test_gff_and_db_are_mutually_exclusive(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit):
        build_parser().parse_args(['-g', 'a.gff', '--db', 'a.db', 't.fa', 'r.fa'])
    assert 'not allowed' in capsys.readouterr().err


@pytest.mark.parametrize(
    ('extra', 'message'),
    [
        (['--min-identity', '0.9', '--copy-identity', '0.5'], 'copy identity'),
        (['--unplaced', 'u.txt'], 'chroms'),
        (['--min-coverage', '1.5'], 'between 0 and 1'),
        (['--threads', '0'], 'threads'),
    ],
)
def test_invalid_combinations_exit_with_parser_error(
    extra: list[str], message: str, capsys: pytest.CaptureFixture[str]
) -> None:
    with pytest.raises(SystemExit) as excinfo:
        main(['-g', 'a.gff', *extra, 't.fa', 'r.fa'])
    assert excinfo.value.code == 2
    assert message in capsys.readouterr().err


def test_config_requires_exactly_one_annotation_source() -> None:
    with pytest.raises(ConfigError):
        LiftoffConfig(target='t', reference='r')


def test_version(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit):
        main(['--version'])
    assert __version__ in capsys.readouterr().out


def test_missing_input_returns_error_status(tmp_path: Path) -> None:
    status = main(
        [
            '-q',
            '-g',
            str(tmp_path / 'missing.gff'),
            '--intermediate-dir',
            str(tmp_path / 'inter'),
            '--unmapped',
            str(tmp_path / 'u.txt'),
            str(tmp_path / 't.fa'),
            str(tmp_path / 'r.fa'),
        ]
    )
    assert status == 1


def test_merge_minimap2_options_uses_whole_tokens() -> None:
    # "--end-bonus" must not be mistaken for "-N"; user values are kept.
    tokens = merge_minimap2_options('--end-bonus 10 -N 5')
    assert tokens == ['--end-bonus', '10', '-N', '5', '-a', '--eqx', '-p', '0.5']
