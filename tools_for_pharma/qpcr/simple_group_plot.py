"""Compatibility wrapper for the generic grouped bar plotting module."""

from __future__ import annotations

from tools_for_pharma.plotting.bar import *  # noqa: F401,F403
from tools_for_pharma.plotting.bar import main


if __name__ == "__main__":
    raise SystemExit(main())
