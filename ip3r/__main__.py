"""Entry point: ``python -m ip3r`` (GUI) or ``python -m ip3r <command>``."""

import sys

from .cli import main

if __name__ == "__main__":
    sys.exit(main())
