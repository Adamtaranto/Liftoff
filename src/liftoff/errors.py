"""Exception hierarchy for Liftoff.

Library code raises these exceptions instead of calling :func:`sys.exit`, so
that callers embedding Liftoff (including the command line interface) decide how
failures are reported.
"""

from __future__ import annotations


class LiftoffError(Exception):
    """Base class for all errors raised by Liftoff."""


class ConfigError(LiftoffError, ValueError):
    """Raised when a combination of options is invalid."""


class InputError(LiftoffError):
    """Raised when an input file is missing, empty or malformed."""


class GffSyntaxError(InputError):
    """Raised when the annotation file cannot be parsed by gffutils.

    Parameters
    ----------
    line_number : int or None
        One-based line number of the first offending record, if it could be
        located.
    message : str
        Human readable description of the problem.
    """

    def __init__(self, line_number: int | None, message: str) -> None:
        self.line_number = line_number
        super().__init__(message)


class DuplicateFeatureIdError(InputError):
    """Raised when distinct features in an annotation share an ``ID``.

    Parameters
    ----------
    duplicates : list of tuple
        ``(feature_id, first_line, duplicate_line)`` for each conflict found.
    message : str
        Human readable description of the problem.
    """

    def __init__(self, duplicates: list[tuple[str, int, int]], message: str) -> None:
        self.duplicates = duplicates
        super().__init__(message)


class AlignmentError(LiftoffError):
    """Raised when an alignment backend fails to produce a SAM file."""
