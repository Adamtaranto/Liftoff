"""Tests for command line parsing and the ``liftoff`` entry points."""

from __future__ import annotations

import logging
from pathlib import Path
import runpy
import sys

import pytest

from liftoff import __version__
from liftoff.cli import build_parser, config_from_args, main
from liftoff.config import LiftoffConfig
import liftoff.pipeline


def parse(*argv: str) -> LiftoffConfig:
    return config_from_args(build_parser().parse_args(argv), argv)


class TestParsing:
    def test_positional_sequences(self) -> None:
        config = parse('-g', 'a.gff', 'target.fa', 'ref.fa')
        assert (config.target, config.reference) == ('target.fa', 'ref.fa')

    def test_command_line_is_recorded(self) -> None:
        assert (
            parse('-g', 'a b.gff', 't.fa', 'r.fa').command_line == "liftoff -g 'a b.gff' t.fa r.fa"
        )

    @pytest.mark.parametrize(
        ('argv', 'attribute', 'value'),
        [
            pytest.param(['--db', 'a.db'], 'db', 'a.db', id='db'),
            pytest.param(['-o', 'out.gff3'], 'output', 'out.gff3', id='output'),
            pytest.param(['-a', '0.9'], 'min_coverage', 0.9, id='min-coverage'),
            pytest.param(['-s', '0.8'], 'min_identity', 0.8, id='min-identity'),
            pytest.param(['-d', '3'], 'distance_factor', 3.0, id='distance-factor'),
            pytest.param(['-p', '4'], 'threads', 4, id='threads'),
            pytest.param(['--copy-identity', '0.95'], 'copy_identity', 0.95, id='copy-identity'),
            pytest.param(['--max-overlap', '0.2'], 'max_overlap', 0.2, id='max-overlap'),
            pytest.param(
                ['--intermediate-dir', 'work'], 'intermediate_dir', 'work', id='intermediate-dir'
            ),
            pytest.param(['--alignments', 'sams'], 'alignments', 'sams', id='alignments'),
            pytest.param(['--copies'], 'copies', True, id='copies'),
            pytest.param(['--polish'], 'polish', True, id='polish'),
            pytest.param(['--no-cds'], 'cds', False, id='no-cds'),
            pytest.param(['--exclude-partial'], 'exclude_partial', True, id='exclude-partial'),
            pytest.param(['--mm2-options=-r 2k'], 'mm2_options', '-r 2k', id='mm2-options'),
        ],
    )
    def test_option_is_mapped_to_config(
        self, argv: list[str], attribute: str, value: object
    ) -> None:
        source = [] if argv[0] == '--db' else ['-g', 'a.gff']
        config = parse(*source, *argv, 't.fa', 'r.fa')
        assert getattr(config, attribute) == value

    def test_annotation_sources_are_mutually_exclusive(
        self, capsys: pytest.CaptureFixture[str]
    ) -> None:
        with pytest.raises(SystemExit):
            build_parser().parse_args(['-g', 'a.gff', '--db', 'a.db', 't.fa', 'r.fa'])
        assert 'not allowed with argument' in capsys.readouterr().err

    def test_verbose_and_quiet_are_mutually_exclusive(self) -> None:
        with pytest.raises(SystemExit):
            build_parser().parse_args(['-g', 'a.gff', '-v', '-q', 't.fa', 'r.fa'])


class TestMain:
    @pytest.fixture
    def pipeline_calls(self, monkeypatch: pytest.MonkeyPatch) -> list[LiftoffConfig]:
        """Replace the pipeline so CLI behaviour is tested in isolation."""
        calls: list[LiftoffConfig] = []
        monkeypatch.setattr(liftoff.pipeline, 'run_liftoff', lambda config: calls.append(config))
        return calls

    def test_runs_pipeline_and_returns_zero(self, pipeline_calls: list[LiftoffConfig]) -> None:
        assert main(['-g', 'a.gff', 't.fa', 'r.fa']) == 0
        assert [config.gff for config in pipeline_calls] == ['a.gff']

    @pytest.mark.parametrize(
        ('flag', 'level'),
        [
            pytest.param([], logging.INFO, id='default'),
            pytest.param(['-v'], logging.DEBUG, id='verbose'),
            pytest.param(['-q'], logging.WARNING, id='quiet'),
        ],
    )
    def test_sets_log_level(
        self, pipeline_calls: list[LiftoffConfig], flag: list[str], level: int
    ) -> None:
        main([*flag, '-g', 'a.gff', 't.fa', 'r.fa'])
        assert logging.getLogger('liftoff').level == level

    @pytest.mark.parametrize(
        ('extra', 'message'),
        [
            pytest.param(
                ['--min-identity', '0.9', '--copy-identity', '0.5'],
                'copy identity',
                id='copy-identity',
            ),
            pytest.param(['--unplaced', 'u.txt'], 'chroms', id='unplaced'),
        ],
    )
    def test_invalid_configuration_is_a_usage_error(
        self,
        pipeline_calls: list[LiftoffConfig],
        capsys: pytest.CaptureFixture[str],
        extra: list[str],
        message: str,
    ) -> None:
        with pytest.raises(SystemExit) as excinfo:
            main(['-g', 'a.gff', *extra, 't.fa', 'r.fa'])
        assert (excinfo.value.code, message in capsys.readouterr().err) == (2, True)

    def test_liftoff_errors_return_status_one(self, tmp_path: Path) -> None:
        argv = [
            '-q',
            '-g',
            str(tmp_path / 'missing.gff'),
            '--intermediate-dir',
            str(tmp_path / 'intermediate'),
            '--unmapped',
            str(tmp_path / 'unmapped.txt'),
            str(tmp_path / 't.fa'),
            str(tmp_path / 'r.fa'),
        ]
        assert main(argv) == 1

    def test_version(self, capsys: pytest.CaptureFixture[str]) -> None:
        with pytest.raises(SystemExit):
            main(['--version'])
        assert capsys.readouterr().out.strip() == f'liftoff {__version__}'


def test_module_entry_point(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    monkeypatch.setattr(sys, 'argv', ['liftoff', '--version'])
    with pytest.raises(SystemExit) as excinfo:
        runpy.run_module('liftoff', run_name='__main__')
    assert (excinfo.value.code, __version__ in capsys.readouterr().out) == (0, True)
