"""Logging helpers.

Liftoff reports progress through the standard :mod:`logging` module. All
loggers are children of the ``"liftoff"`` logger, and log records are written to
standard error so that annotation output written to standard output is never
interleaved with progress messages.
"""

from __future__ import annotations

import logging
import sys

#: Name of the package-level logger that all module loggers descend from.
ROOT_LOGGER_NAME = 'liftoff'

_LOG_FORMAT = '%(asctime)s %(levelname)-7s %(name)s: %(message)s'
_DATE_FORMAT = '%H:%M:%S'


def get_logger(name: str) -> logging.Logger:
    """Return a logger nested under the package logger.

    Parameters
    ----------
    name : str
        Usually ``__name__`` of the calling module.

    Returns
    -------
    logging.Logger
        The logger for ``name``. Names that do not already start with
        ``"liftoff"`` are prefixed so that configuration of the package logger
        applies to them.
    """
    if name != ROOT_LOGGER_NAME and not name.startswith(ROOT_LOGGER_NAME + '.'):
        name = f'{ROOT_LOGGER_NAME}.{name}'
    return logging.getLogger(name)


def configure_logging(level: int = logging.INFO) -> None:
    """Attach a standard-error handler to the package logger.

    Calling this more than once replaces the previously installed handler
    rather than duplicating output.

    Parameters
    ----------
    level : int, default logging.INFO
        Minimum severity that is emitted.
    """
    logger = logging.getLogger(ROOT_LOGGER_NAME)
    for handler in list(logger.handlers):
        if getattr(handler, '_liftoff_handler', False):
            logger.removeHandler(handler)
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(logging.Formatter(_LOG_FORMAT, _DATE_FORMAT))
    handler._liftoff_handler = True  # type: ignore[attr-defined]
    logger.addHandler(handler)
    logger.setLevel(level)
    # Keep records from propagating to the root logger to avoid duplicates when
    # an application has configured logging itself.
    logger.propagate = False
