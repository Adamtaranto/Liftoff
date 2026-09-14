"""Allow ``python -m liftoff``."""

from __future__ import annotations

import sys

from liftoff.cli import main

if __name__ == '__main__':
    sys.exit(main())
