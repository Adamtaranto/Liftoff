"""Tests for logging helpers."""

from __future__ import annotations

import logging

import pytest

from liftoff.log import configure_logging, get_logger


@pytest.mark.parametrize(
    ('name', 'expected'),
    [
        pytest.param('liftoff', 'liftoff', id='package'),
        pytest.param('liftoff.cli', 'liftoff.cli', id='submodule'),
        pytest.param('plugin', 'liftoff.plugin', id='prefixed'),
    ],
)
def test_get_logger_nests_under_package(name: str, expected: str) -> None:
    assert get_logger(name).name == expected


def test_configure_logging_sets_level(restore_package_logger: logging.Logger) -> None:
    configure_logging(logging.WARNING)
    assert restore_package_logger.level == logging.WARNING


def test_configure_logging_does_not_duplicate_handlers(
    restore_package_logger: logging.Logger,
) -> None:
    configure_logging()
    configure_logging()
    installed = [
        h for h in restore_package_logger.handlers if getattr(h, '_liftoff_handler', False)
    ]
    assert len(installed) == 1


def test_messages_go_to_stderr(
    restore_package_logger: logging.Logger, capsys: pytest.CaptureFixture[str]
) -> None:
    configure_logging(logging.INFO)
    get_logger('liftoff.test').info('hello')
    captured = capsys.readouterr()
    assert ('hello' in captured.err, captured.out) == (True, '')
