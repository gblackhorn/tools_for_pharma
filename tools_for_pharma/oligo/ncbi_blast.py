"""Compatibility facade and executable module for transcript scanning.

Implementation now lives in the ``transcript_scan`` package.  Public imports
remain available here until the compatibility surface is finalized in Phase 11.
"""

from __future__ import annotations

import sys

from tools_for_pharma.oligo.transcript_scan.cli import main as run_cli
from tools_for_pharma.oligo.transcript_scan.gui import *  # noqa: F403
from tools_for_pharma.oligo.transcript_scan.gui import run_gui


def main() -> int:
    """Preserve the historical module entry point through the CLI interface."""
    return run_cli(gui_runner=run_gui)


if __name__ == "__main__":
    raise SystemExit(main())
