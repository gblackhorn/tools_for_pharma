"""Friendly entrypoint for IC50/IC75/IC90 fitted dose-response plotting."""

from __future__ import annotations

from tools_for_pharma.plotting.curve import *  # noqa: F401,F403
from tools_for_pharma.plotting.curve import main as curve_main
from tools_for_pharma.plotting.ic50_summary import add_summary_panel


def main() -> int:
    return curve_main(summary_renderer=add_summary_panel)


if __name__ == "__main__":
    raise SystemExit(main())
