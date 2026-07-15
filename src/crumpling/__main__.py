"""Run the command-line interface with ``python -m crumpling``."""

import sys

from crumpling import cli


if __name__ == "__main__":
    sys.exit(cli.main())
