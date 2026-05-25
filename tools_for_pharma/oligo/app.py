"""Combined entrypoint for oligo sequence utilities.

Prefer these clearer user-facing utilities:
  util_oligo_offtarget.py
    Prepare antisense region information for in silico off-target analysis.

  util_transcript_sequence.py
    Extract SS/AS oligo design sequences from long transcripts.

Usage examples:
  Convert one antisense sequence region to its complementary sense strand:
    python util_oligo_offtarget.py "AUGCUACGGAUCUAGCUAGCU"
    python util_oligo_offtarget.py "AUGCUACGGAUCUAGCUAGCU" --start 3 --end 19

  Extract one transcript range and output SS/AS:
    python util_transcript_sequence.py --transcript-file transcript.fasta --start 120 --end 140

  Match many transcript ranges from an Excel/CSV table with start/end columns:
    python util_transcript_sequence.py --transcript-file transcript.fasta --range-table ranges.xlsx

  Process many antisense sequences from one Excel/CSV column:
    python util_oligo_offtarget.py --input oligos.xlsx --column antisense

  Open the combined file-picker GUI:
    python oligo_utils.py --gui

Batch mode appends result columns to the original table by default. It adds
normalized_antisense, antisense_START-END_5to3, sense_END-START_5to3, and error.

Transcript mode accepts FASTA or plain text and outputs SS/AS for the selected
1-based inclusive transcript range.

Transcript table mode reads start/end from the first sheet by default and fills
SS matched, AS matched, matched length (nt), and match error.
"""

from __future__ import annotations

import argparse
import sys

from tools_for_pharma.oligo.core import (
    DEFAULT_END,
    DEFAULT_START,
    antisense_region_to_sense,
)
from tools_for_pharma.oligo.gui import run_gui
from tools_for_pharma.oligo.table import output_column_names, process_table
from tools_for_pharma.oligo.transcript import (
    format_transcript_oligo_result,
    transcript_file_to_oligo,
    transcript_text_to_oligo,
)
from tools_for_pharma.oligo.transcript_table import process_transcript_range_table


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser."""
    parser = argparse.ArgumentParser(
        description="Extract an antisense region and compute complementary sense 5'->3'."
    )
    parser.add_argument("antisense", nargs="?", help="Antisense sequence in 5'->3'.")
    parser.add_argument(
        "-i",
        "--input",
        help=(
            "CSV or Excel table containing many antisense sequences. Use --column "
            "to choose the sequence column; appends normalized_antisense, "
            "antisense_START-END_5to3, sense_END-START_5to3, and error columns."
        ),
    )
    parser.add_argument(
        "--transcript-file",
        help="FASTA or plain text transcript file. Extract --start to --end as SS/AS.",
    )
    parser.add_argument(
        "--transcript",
        help="Pasted FASTA or plain transcript sequence. Extract --start to --end as SS/AS.",
    )
    parser.add_argument(
        "--range-table",
        help=(
            "Excel/CSV table with start and end columns. Requires --transcript-file "
            "and fills SS matched, AS matched, matched length (nt), and match error."
        ),
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
        help=f"Start position. Defaults to {DEFAULT_START}.",
    )
    parser.add_argument(
        "--end",
        type=int,
        help=f"End position. Defaults to {DEFAULT_END}.",
    )
    parser.add_argument(
        "--gui",
        action="store_true",
        help=(
            "Open a file-picker GUI for transcript range extraction, transcript "
            "range-table matching, or batch antisense-table processing."
        ),
    )
    return parser


def resolve_positions(args: argparse.Namespace) -> tuple[int, int]:
    """Return start/end positions using legacy defaults when omitted."""
    return (
        args.start if args.start is not None else DEFAULT_START,
        args.end if args.end is not None else DEFAULT_END,
    )


def require_transcript_positions(
    parser: argparse.ArgumentParser,
    args: argparse.Namespace,
) -> tuple[int, int]:
    """Return explicit transcript coordinates or stop with a CLI error."""
    if args.start is None or args.end is None:
        parser.error("transcript mode requires both --start and --end")
    return args.start, args.end


def main() -> int:
    """Run CLI, batch, or GUI mode."""
    parser = build_parser()
    args = parser.parse_args()

    try:
        if args.gui:
            return run_gui()

        if args.transcript_file and args.transcript:
            parser.error("use either --transcript-file or --transcript, not both")

        if args.range_table and not args.transcript_file:
            parser.error("--range-table requires --transcript-file")

        if args.range_table:
            output = process_transcript_range_table(
                transcript_path=args.transcript_file,
                table_path=args.range_table,
                output_path=args.output,
                sheet_name=args.sheet,
            )
            print(f"Wrote transcript range results to: {output}")
            return 0

        if args.transcript_file:
            start, end = require_transcript_positions(parser, args)
            result = transcript_file_to_oligo(
                args.transcript_file,
                start=start,
                end=end,
            )
            print(format_transcript_oligo_result(result))
            return 0

        if args.transcript:
            start, end = require_transcript_positions(parser, args)
            result = transcript_text_to_oligo(
                args.transcript,
                start=start,
                end=end,
            )
            print(format_transcript_oligo_result(result))
            return 0

        if args.input:
            start, end = resolve_positions(args)
            output = process_table(
                args.input,
                column=args.column,
                output_path=args.output,
                sheet_name=args.sheet,
                start=start,
                end=end,
            )
            action = "Wrote" if args.output else "Appended"
            print(f"{action} results to: {output}")
            return 0

        if args.antisense:
            start, end = resolve_positions(args)
            antisense_region, sense_5to3 = antisense_region_to_sense(
                args.antisense,
                start=start,
                end=end,
            )
            antisense_column, sense_column = output_column_names(start, end)
            print(f"{antisense_column}: {antisense_region}")
            print(f"{sense_column}: {sense_5to3}")
            return 0

        parser.error(
            "provide an antisense sequence, --input, --transcript-file, --transcript, or --gui"
        )
        return 2
    except ValueError as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
