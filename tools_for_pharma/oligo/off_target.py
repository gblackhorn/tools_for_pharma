"""Utility for oligo off-target analysis preparation.

Purpose:
  Extract the antisense seed/core region, usually positions 2-18, and compute
  the complementary sense sequence in 5'->3' orientation. This is the workflow
  used to prepare sequence information for in silico off-target analysis.

Examples:
  Process one antisense sequence:
    python util_oligo_offtarget.py "AUGCUACGGAUCUAGCUAGCU"

  Use a non-default antisense region:
    python util_oligo_offtarget.py "AUGCUACGGAUCUAGCUAGCU" --start 3 --end 19

  Process many antisense sequences from an Excel/CSV table:
    python util_oligo_offtarget.py --input oligos.xlsx --column antisense

  By default, table mode appends result columns to the selected input sheet.
  Use --output only when you want to write a separate file.

  Open the off-target table GUI:
    python util_oligo_offtarget.py --gui
"""

from __future__ import annotations

import argparse
import sys

from tools_for_pharma.oligo.core import (
    DEFAULT_END,
    DEFAULT_START,
    antisense_region_to_sense,
)
from tools_for_pharma.oligo.gui import run_offtarget_gui
from tools_for_pharma.oligo.table import output_column_names, process_table


def build_parser() -> argparse.ArgumentParser:
    """Build the off-target utility parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Prepare oligo antisense region information for in silico "
            "off-target analysis."
        )
    )
    parser.add_argument("antisense", nargs="?", help="Antisense sequence in 5'->3'.")
    parser.add_argument(
        "-i",
        "--input",
        help="CSV or Excel table containing many antisense sequences.",
    )
    parser.add_argument(
        "-c",
        "--column",
        default="antisense",
        help="Column header containing antisense sequences in batch mode.",
    )
    parser.add_argument("-o", "--output", help="Output .xlsx or .csv path.")
    parser.add_argument("--sheet", help="Excel sheet name. Defaults to the first sheet.")
    parser.add_argument(
        "--start",
        type=int,
        default=DEFAULT_START,
        help=f"Antisense start position. Defaults to {DEFAULT_START}.",
    )
    parser.add_argument(
        "--end",
        type=int,
        default=DEFAULT_END,
        help=f"Antisense end position. Defaults to {DEFAULT_END}.",
    )
    parser.add_argument("--gui", action="store_true", help="Open the off-target table GUI.")
    return parser


def main() -> int:
    """Run the off-target utility."""
    args = build_parser().parse_args()

    try:
        if args.gui:
            return run_offtarget_gui()

        if args.input:
            output = process_table(
                args.input,
                column=args.column,
                output_path=args.output,
                sheet_name=args.sheet,
                start=args.start,
                end=args.end,
            )
            action = "Wrote" if args.output else "Appended"
            print(f"{action} off-target preparation results to: {output}")
            return 0

        if args.antisense:
            antisense_region, sense_5to3 = antisense_region_to_sense(
                args.antisense,
                start=args.start,
                end=args.end,
            )
            antisense_column, sense_column = output_column_names(args.start, args.end)
            print(f"{antisense_column}: {antisense_region}")
            print(f"{sense_column}: {sense_5to3}")
            return 0

        build_parser().error("provide an antisense sequence, --input, or --gui")
        return 2
    except ValueError as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
