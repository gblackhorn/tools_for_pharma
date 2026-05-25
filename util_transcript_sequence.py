"""Utility for extracting oligo design sequences from long transcripts.

Purpose:
  Read an NCBI-style FASTA/plain-text transcript, extract transcript positions
  by 1-based inclusive start/end coordinates, and output matched SS/AS
  sequences for oligo design.

Examples:
  Extract one transcript range:
    python util_transcript_sequence.py --transcript-file transcript.fasta --start 120 --end 140

  Match many transcript ranges from an Excel/CSV table with start/end columns:
    python util_transcript_sequence.py --transcript-file transcript.fasta --range-table ranges.xlsx

  Paste a transcript directly:
    python util_transcript_sequence.py --transcript "AGACGCCTGGGAACTGCGGCC" --start 1 --end 21

  Open the transcript sequence GUI:
    python util_transcript_sequence.py --gui
"""

from __future__ import annotations

import argparse
import sys

from oligo_gui import run_transcript_sequence_gui
from oligo_transcript import (
    format_transcript_oligo_result,
    transcript_file_to_oligo,
    transcript_text_to_oligo,
)
from oligo_transcript_table import process_transcript_range_table


def build_parser() -> argparse.ArgumentParser:
    """Build the transcript sequence utility parser."""
    parser = argparse.ArgumentParser(
        description="Extract SS/AS oligo design sequences from long transcripts."
    )
    parser.add_argument(
        "--transcript-file",
        help="FASTA or plain text transcript file.",
    )
    parser.add_argument(
        "--transcript",
        help="Pasted FASTA or plain transcript sequence.",
    )
    parser.add_argument(
        "--range-table",
        help=(
            "Excel/CSV table with start and end columns. Requires --transcript-file "
            "and fills SS matched, AS matched, matched length (nt), and match error."
        ),
    )
    parser.add_argument("-o", "--output", help="Output .xlsx, .csv, or .txt path.")
    parser.add_argument("--sheet", help="Excel sheet name. Defaults to the first sheet.")
    parser.add_argument("--start", type=int, help="Transcript start position.")
    parser.add_argument("--end", type=int, help="Transcript end position.")
    parser.add_argument(
        "--gui",
        action="store_true",
        help="Open the transcript sequence file-picker GUI.",
    )
    return parser


def require_positions(
    parser: argparse.ArgumentParser,
    args: argparse.Namespace,
) -> tuple[int, int]:
    """Return explicit transcript coordinates or stop with a CLI error."""
    if args.start is None or args.end is None:
        parser.error("single transcript range mode requires both --start and --end")
    return args.start, args.end


def main() -> int:
    """Run the transcript sequence utility."""
    parser = build_parser()
    args = parser.parse_args()

    try:
        if args.gui:
            return run_transcript_sequence_gui()

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
            start, end = require_positions(parser, args)
            result = transcript_file_to_oligo(
                args.transcript_file,
                start=start,
                end=end,
            )
            output_text = format_transcript_oligo_result(result)
            if args.output:
                with open(args.output, "w", encoding="utf-8") as output_file:
                    output_file.write(f"{output_text}\n")
                print(f"Wrote transcript range result to: {args.output}")
            else:
                print(output_text)
            return 0

        if args.transcript:
            start, end = require_positions(parser, args)
            result = transcript_text_to_oligo(
                args.transcript,
                start=start,
                end=end,
            )
            print(format_transcript_oligo_result(result))
            return 0

        parser.error("provide --transcript-file, --transcript, --range-table, or --gui")
        return 2
    except ValueError as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
