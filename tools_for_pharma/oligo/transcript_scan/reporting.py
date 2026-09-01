"""Text, CSV, and workbook reporting for transcript-scan workflows."""

from __future__ import annotations

import argparse
import csv
import io
from pathlib import Path
from typing import Iterable, Protocol

from tools_for_pharma.oligo.ncbi_transport import (
    DEFAULT_POLL_SECONDS,
    DEFAULT_REQUEST_SECONDS,
)
from tools_for_pharma.oligo.transcript_scan.models import (
    AntisenseQuery,
    AntisenseRegion,
    ComparisonResult,
    PrivatePanelScanResult,
    QueryTargetSummary,
    TranscriptMatch,
    TranscriptTargetResult,
)
from tools_for_pharma.oligo.transcript_scan.queries import (
    assign_unique_blast_query_ids,
    duplicate_sequence_groups,
    normalize_sequence_type,
)
from tools_for_pharma.oligo.transcript_scan.remote_blast import (
    CSV_COLUMNS,
    filter_blast_rows,
    parse_blast_csv,
)
from tools_for_pharma.sequence.nucleotides import normalize_rna
from tools_for_pharma.shared.excel_utils import sanitize_sheet_name


class BlastBatchResultLike(Protocol):
    """Reporting fields supplied by one completed remote BLAST batch."""

    batch_index: int
    submission: object
    queries: tuple[AntisenseQuery, ...]
    csv_text: str


def input_query_rows(records: list[AntisenseQuery]) -> list[dict[str, object]]:
    """Return input query rows with duplicate and source-field annotations."""
    prepared_records = assign_unique_blast_query_ids(records)
    duplicate_groups = duplicate_sequence_groups(prepared_records)
    rows = []
    for index, record in enumerate(prepared_records, start=1):
        sequence = normalize_rna(record.sequence_5to3)
        duplicate_names = duplicate_groups.get(sequence, [])
        output_row = {
            "input_order": index,
            "sequence_type": normalize_sequence_type(record.sequence_type),
            "antisense_name": record.name,
            "blast_query_id": record.blast_query_id,
            "antisense_5to3": sequence,
            "length_nt": len(sequence),
            "target_accession": record.target_accession,
            "target_gene": record.target_gene,
            "species": record.species,
            "notes": record.notes,
            "is_duplicate_sequence": bool(duplicate_names),
            "duplicate_group_names": ";".join(duplicate_names),
        }
        for column, value in record.source_fields.items():
            if column not in output_row:
                output_row[column] = value
        rows.append(output_row)
    return rows


def transcript_matches_to_csv(matches: Iterable[TranscriptMatch]) -> str:
    """Format local transcript matches as CSV text."""
    output = io.StringIO()
    writer = csv.writer(output, lineterminator="\n")
    writer.writerow(
        [
            "transcript_name",
            "antisense_name",
            "scan_region",
            "as_region_start",
            "as_region_end",
            "antisense_5to3",
            "antisense_region_5to3",
            "expected_target_5to3",
            "transcript_start",
            "transcript_end",
            "mismatches",
            "transcript_window_5to3",
            "transcript_match_as_5to3",
            "mismatch_positions_1based",
            "as_mismatch_positions_1based",
            "sequence_type",
        ]
    )
    for match in matches:
        writer.writerow(
            [
                match.transcript_name,
                match.antisense_name,
                match.scan_region,
                match.as_region_start,
                match.as_region_end,
                match.antisense_5to3,
                match.antisense_region_5to3,
                match.target_5to3,
                match.transcript_start,
                match.transcript_end,
                match.mismatches,
                match.transcript_window_5to3,
                match.transcript_match_as_5to3,
                ";".join(
                    str(position) for position in match.mismatch_positions_1based
                ),
                ";".join(
                    str(position) for position in match.as_mismatch_positions_1based
                ),
                match.sequence_type,
            ]
        )
    return output.getvalue()


def terminal_table(headers: list[str], rows: list[list[object]]) -> str:
    """Return a compact plain-text table for quick terminal review."""
    if not rows:
        return ""
    text_rows = [[str(value) for value in row] for row in rows]
    widths = [
        max(len(header), *(len(row[index]) for row in text_rows))
        for index, header in enumerate(headers)
    ]
    output = io.StringIO()
    output.write(
        "  ".join(header.ljust(widths[index]) for index, header in enumerate(headers))
    )
    output.write("\n")
    output.write("  ".join("-" * width for width in widths))
    output.write("\n")
    for row in text_rows:
        output.write(
            "  ".join(value.ljust(widths[index]) for index, value in enumerate(row))
        )
        output.write("\n")
    return output.getvalue().rstrip()


def transcript_match_terminal_table(matches: Iterable[TranscriptMatch]) -> str:
    match_list = list(matches)
    sequence_types = {match.sequence_type for match in match_list}
    ss_only = sequence_types == {"SS"}
    matched_header = "matched_ss_5to3" if ss_only else "matched_as_5to3"
    query_mm_header = "ss_mm_pos" if ss_only else "as_mm_pos"
    headers = [
        "query",
        "region",
        "start",
        "end",
        "mm",
        "target_5to3",
        matched_header,
        "mm_pos",
        query_mm_header,
    ]
    rows = [
        [
            match.antisense_name,
            match.scan_region,
            match.transcript_start,
            match.transcript_end,
            match.mismatches,
            match.transcript_window_5to3,
            match.transcript_match_as_5to3,
            ";".join(str(position) for position in match.mismatch_positions_1based)
            or "-",
            ";".join(str(position) for position in match.as_mismatch_positions_1based)
            or "-",
        ]
        for match in match_list
    ]
    return terminal_table(headers, rows)


def format_transcript_matches_for_terminal(
    matches: Iterable[TranscriptMatch],
    queries: list[AntisenseQuery],
    scan_regions: list[AntisenseRegion],
    max_mismatches: int,
) -> str:
    """Format local transcript matches as a readable terminal summary."""
    match_list = list(matches)
    transcript_names = sorted({match.transcript_name for match in match_list})
    sequence_types = {normalize_sequence_type(query.sequence_type) for query in queries}
    label = next(iter(sequence_types)) if len(sequence_types) == 1 else "Oligo"
    output = io.StringIO()
    output.write("Local transcript scan\n")
    output.write(f"{label} queries: {len(queries)}\n")
    if len(queries) == 1:
        query = queries[0]
        query_label = normalize_sequence_type(query.sequence_type)
        output.write(f"{query_label} name: {query.name}\n")
        output.write(f"{query_label} sequence: {normalize_rna(query.sequence_5to3)}\n")
    if transcript_names:
        if len(transcript_names) == 1:
            output.write(f"Transcript: {transcript_names[0]}\n")
        else:
            output.write(f"Transcripts: {len(transcript_names)}\n")
    output.write(f"Scan regions: {', '.join(region.name for region in scan_regions)}\n")
    output.write(f"Max mismatches: {max_mismatches}\n")
    output.write(f"Matches: {len(match_list)}\n")

    if not match_list:
        output.write("\nNo local transcript matches found within mismatch threshold.")
        return output.getvalue()

    output.write("\n")
    output.write(transcript_match_terminal_table(match_list))
    return output.getvalue()


def format_closest_transcript_matches_for_terminal(
    matches: Iterable[TranscriptMatch],
    closest: int,
    max_mismatches: int,
) -> str:
    match_list = list(matches)
    output = io.StringIO()
    output.write("\nClosest transcript windows\n")
    output.write(f"Showing: {len(match_list)}")
    if len(match_list) == closest:
        output.write(f" of top {closest}")
    output.write("\n")
    output.write(f"These are not filtered by --max-mismatches {max_mismatches}.\n")
    if not match_list:
        output.write("\nNo transcript windows available.")
        return output.getvalue()
    output.write("\n")
    output.write(transcript_match_terminal_table(match_list))
    return output.getvalue()


def transcript_match_rows(
    matches: Iterable[TranscriptMatch],
    *,
    include_as_oriented_match: bool = True,
) -> list[dict[str, object]]:
    """Return local-match rows, optionally including the AS-oriented match."""
    rows = []
    for match in matches:
        row = {
            "transcript_name": match.transcript_name,
            "sequence_type": match.sequence_type,
            "antisense_name": match.antisense_name,
            "scan_region": match.scan_region,
            "as_region_start": match.as_region_start,
            "as_region_end": match.as_region_end,
            "antisense_5to3": match.antisense_5to3,
            "antisense_region_5to3": match.antisense_region_5to3,
            "expected_target_5to3": match.target_5to3,
            "transcript_start": match.transcript_start,
            "transcript_end": match.transcript_end,
            "mismatches": match.mismatches,
            "transcript_window_5to3": match.transcript_window_5to3,
        }
        if include_as_oriented_match:
            row["transcript_match_as_5to3"] = match.transcript_match_as_5to3
        row.update(
            {
                "mismatch_positions_1based": ";".join(
                    str(position) for position in match.mismatch_positions_1based
                ),
                "as_mismatch_positions_1based": ";".join(
                    str(position) for position in match.as_mismatch_positions_1based
                ),
            }
        )
        rows.append(row)
    return rows


def comparison_result_rows(
    results: Iterable[ComparisonResult],
) -> list[dict[str, object]]:
    """Project compact comparison results into the agreed workbook schema."""
    return [
        {
            "input_order": result.input_order,
            "query_name": result.query_name,
            "target_accession": result.target_accession,
            "scan_region": result.scan_region,
            "region_start": result.region_start,
            "region_end": result.region_end,
            "result": result.result,
            "sites_within_threshold": result.sites_within_threshold,
            "best_mismatches": result.best_mismatches,
            "mismatch_positions_in_query_1based": ";".join(
                str(position)
                for position in result.mismatch_positions_in_query_1based
            ),
            "best_transcript_start": result.best_transcript_start,
            "best_transcript_end": result.best_transcript_end,
            "query_region_5to3": result.query_region_5to3,
            "best_match_in_query_orientation_5to3": (
                result.best_match_in_query_orientation_5to3
            ),
            "differences": result.differences,
        }
        for result in results
    ]


def transcript_target_rows(
    targets: Iterable[TranscriptTargetResult],
) -> list[dict[str, object]]:
    """Return public transcript retrieval metadata without duplicating sequences."""
    return [
        {
            "requested_accession": target.requested_accession,
            "retrieved_accession": target.retrieved_accession,
            "transcript_name": target.transcript_name,
            "sequence_length_nt": target.sequence_length_nt,
            "cache_path": target.cache_path,
            "cache_status": target.cache_status,
            "exact_version_match": target.exact_version_match,
            "sequence_sha256": target.sequence_sha256,
            "retrieved_at_utc": target.retrieved_at_utc,
            "status": target.status,
            "error": target.error,
        }
        for target in targets
    ]


def query_target_summary_rows(
    summaries: Iterable[QueryTargetSummary],
) -> list[dict[str, object]]:
    return [
        {
            "query_name": summary.query_name,
            "sequence_type": summary.sequence_type,
            "requested_accession": summary.requested_accession,
            "retrieved_accession": summary.retrieved_accession,
            "target_status": summary.target_status,
            "scan_status": summary.scan_status,
            "scan_regions": summary.scan_regions,
            "match_count": summary.match_count,
            "exact_match_count": summary.exact_match_count,
            "best_mismatches": summary.best_mismatches,
            "error": summary.error,
        }
        for summary in summaries
    ]


def query_length_by_blast_id(
    queries: Iterable[AntisenseQuery],
) -> dict[str, int]:
    return {
        query.blast_query_id: len(normalize_rna(query.sequence_5to3))
        for query in assign_unique_blast_query_ids(queries)
    }


def blast_raw_rows(
    batch_results: Iterable[BlastBatchResultLike],
    queries: Iterable[AntisenseQuery],
) -> list[dict[str, object]]:
    query_lengths = query_length_by_blast_id(queries)
    rows = []
    for result in batch_results:
        for row in parse_blast_csv(result.csv_text):
            query_length = query_lengths.get(row["query_id"])
            alignment_length = int(float(row["alignment_length"]))
            rows.append(
                {
                    "rid": getattr(result.submission, "rid"),
                    "batch_index": result.batch_index,
                    **row,
                    "query_length": query_length,
                    "alignment_fraction": (
                        alignment_length / query_length if query_length else None
                    ),
                }
            )
    return rows


def blast_batch_rows(
    batch_results: Iterable[BlastBatchResultLike],
) -> list[dict[str, object]]:
    rows = []
    for result in batch_results:
        sequences = [normalize_rna(query.sequence_5to3) for query in result.queries]
        rows.append(
            {
                "batch_index": result.batch_index,
                "rid": getattr(result.submission, "rid"),
                "rtoe_seconds": getattr(result.submission, "rtoe_seconds"),
                "query_count": len(result.queries),
                "total_query_bases": sum(len(sequence) for sequence in sequences),
                "query_names": ";".join(query.name for query in result.queries),
            }
        )
    return rows


def metadata_rows(
    args: argparse.Namespace,
    queries: list[AntisenseQuery],
    scan_regions: list[AntisenseRegion],
    started_at: str,
    completed_at: str,
    panel_targets: Iterable[TranscriptTargetResult] | None = None,
) -> list[dict[str, object]]:
    target_list = list(panel_targets or [])
    metadata = {
        "started_at_utc": started_at,
        "completed_at_utc": completed_at,
        "tool": args.tool,
        "email": args.email,
        "database": args.database,
        "expect": args.expect,
        "word_size": args.word_size,
        "hitlist_size": args.hitlist_size,
        "megablast": args.megablast,
        "request_seconds": max(args.request_seconds, DEFAULT_REQUEST_SECONDS),
        "poll_seconds": max(args.poll_seconds, DEFAULT_POLL_SECONDS),
        "max_batch_bases": args.max_batch_bases,
        "query_count": len(queries),
        "sequence_types": ";".join(
            sorted(
                {
                    normalize_sequence_type(query.sequence_type)
                    for query in queries
                }
            )
        ),
        "total_query_bases": sum(
            len(normalize_rna(query.sequence_5to3)) for query in queries
        ),
        "scan_regions": ";".join(region.name for region in scan_regions),
        "max_mismatches_local_scan": args.max_mismatches,
        "blast_filter_max_mismatches": args.filter_max_mismatches,
        "blast_filter_max_gap_opens": args.filter_max_gap_opens,
        "blast_filter_min_alignment_fraction": args.filter_min_alignment_fraction,
        "privacy_mode": (
            "remote_blast_query_submission"
            if args.blast or args.blast_only
            else "local_guide_scan"
        ),
        "guide_sequence_transmitted_to_ncbi": bool(args.blast or args.blast_only),
        "private_panel_mode": bool(getattr(args, "private_panel", False)),
        "offline_mode": bool(getattr(args, "offline", False)),
        "refresh_targets": bool(getattr(args, "refresh_targets", False)),
        "panel_target_count": len(target_list),
        "panel_targets_ready": sum(
            target.status == "ready" for target in target_list
        ),
        "panel_targets_error": sum(target.status == "error" for target in target_list),
    }
    return [{"key": key, "value": value} for key, value in metadata.items()]


def write_excel_workbook(
    path: Path,
    sheets: dict[str, list[dict[str, object]]],
) -> None:
    import pandas as pd

    path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(path) as writer:
        for sheet_name, rows in sheets.items():
            pd.DataFrame(rows).to_excel(
                writer,
                sheet_name=sanitize_sheet_name(sheet_name),
                index=False,
            )


def default_result_workbook(args: argparse.Namespace) -> Path | None:
    if args.result_workbook:
        return args.result_workbook
    source = (
        args.as_table
        or args.as_file
        or getattr(args, "ss_table", None)
        or getattr(args, "ss_file", None)
    )
    if source and not args.output and not args.blast_output:
        workflow = (
            "ncbi_blast" if args.blast or args.blast_only else "ncbi_transcript_scan"
        )
        return source.with_name(f"{source.stem}_{workflow}_results.xlsx")
    return None


def default_private_panel_workbook(args: argparse.Namespace) -> Path:
    if args.result_workbook:
        return args.result_workbook
    source = (
        args.as_table
        or args.as_file
        or getattr(args, "ss_table", None)
        or getattr(args, "ss_file", None)
        or getattr(args, "target_table", None)
    )
    if source:
        return source.with_name(f"{source.stem}_private_transcript_panel_results.xlsx")
    return Path("private_transcript_panel_results.xlsx")


def default_gui_result_workbook(input_file: Path) -> Path:
    return input_file.with_name(
        f"{input_file.stem}_ncbi_transcript_scan_results.xlsx"
    )


def write_result_workbook(
    path: Path,
    args: argparse.Namespace,
    queries: list[AntisenseQuery],
    scan_regions: list[AntisenseRegion],
    local_matches: list[TranscriptMatch],
    blast_results: list[BlastBatchResultLike],
    started_at: str,
    completed_at: str,
    *,
    include_blast_sheets: bool = True,
    comparison_results: list[ComparisonResult] | None = None,
    transcript_targets: list[TranscriptTargetResult] | None = None,
    query_target_summaries: list[QueryTargetSummary] | None = None,
    closest_local_matches: list[TranscriptMatch] | None = None,
) -> None:
    raw_blast_rows = blast_raw_rows(blast_results, queries)
    sheets = {"input_queries": input_query_rows(queries)}
    if comparison_results is not None:
        sheets["comparison_results"] = comparison_result_rows(comparison_results)
        sheets["local_transcript_scan"] = transcript_match_rows(
            local_matches,
            include_as_oriented_match=False,
        )
        if transcript_targets is not None:
            sheets["transcript_targets"] = transcript_target_rows(transcript_targets)
    else:
        if transcript_targets is not None:
            sheets["transcript_targets"] = transcript_target_rows(transcript_targets)
        sheets["local_transcript_scan"] = transcript_match_rows(
            local_matches,
            include_as_oriented_match=False,
        )
    if query_target_summaries is not None:
        sheets["query_target_summary"] = query_target_summary_rows(
            query_target_summaries
        )
    if closest_local_matches is not None:
        sheets["closest_transcript_windows"] = transcript_match_rows(
            closest_local_matches,
            include_as_oriented_match=False,
        )
    if include_blast_sheets:
        sheets.update(
            {
                "blast_hits_raw": raw_blast_rows,
                "blast_hits_filtered": filter_blast_rows(
                    raw_blast_rows,
                    args.filter_max_mismatches,
                    args.filter_max_gap_opens,
                    args.filter_min_alignment_fraction,
                ),
                "blast_batches": blast_batch_rows(blast_results),
            }
        )
    sheets["run_metadata"] = metadata_rows(
        args,
        queries,
        scan_regions,
        started_at,
        completed_at,
        transcript_targets,
    )
    write_excel_workbook(path, sheets)


def format_single_sequence_scan_result(
    args: argparse.Namespace,
    queries: list[AntisenseQuery],
    scan_regions: list[AntisenseRegion],
    result: PrivatePanelScanResult,
) -> str:
    """Format a compact, copyable single-sequence result for the GUI."""
    query = queries[0]
    target = result.targets[0]
    summary = result.summaries[0]
    output = io.StringIO()
    output.write("LOCAL SINGLE-SEQUENCE TRANSCRIPT SCAN\n")
    output.write("=" * 37 + "\n\n")
    output.write(f"Sequence type: {normalize_sequence_type(query.sequence_type)}\n")
    output.write(f"Sequence name: {query.name}\n")
    output.write(f"Sequence 5'->3': {normalize_rna(query.sequence_5to3)}\n")
    if target.retrieved_accession:
        output.write(f"Transcript accession: {target.retrieved_accession}\n")
    output.write(f"Transcript: {target.transcript_name}\n")
    output.write(f"Transcript length: {target.sequence_length_nt} nt\n")
    output.write(f"Transcript source: {target.cache_status}\n")
    if target.cache_path:
        source_label = "Cache file" if target.retrieved_accession else "Target file"
        output.write(f"{source_label}: {target.cache_path}\n")
    output.write(
        "NCBI transcript retrieval: "
        + (
            "Used accession/cache workflow"
            if target.retrieved_accession
            else "Not used"
        )
        + "\n"
    )
    output.write("Guide sequence sent to NCBI: No\n")
    output.write(
        "Scan regions: " + ", ".join(region.name for region in scan_regions) + "\n"
    )
    output.write(f"Maximum mismatches: {args.max_mismatches}\n")
    output.write(f"Matches within threshold: {summary.match_count}\n")
    output.write(f"Exact matches: {summary.exact_match_count}\n")
    output.write(
        "Best mismatch count across selected regions: "
        f"{summary.best_mismatches if summary.best_mismatches is not None else '-'}\n"
    )

    output.write("\nMATCHES WITHIN THRESHOLD\n")
    output.write("-" * 24 + "\n")
    if result.matches:
        output.write(transcript_match_terminal_table(result.matches))
        output.write("\n")
    else:
        output.write("No transcript windows passed the mismatch threshold.\n")

    output.write("\nCLOSEST TRANSCRIPT WINDOWS\n")
    output.write("-" * 26 + "\n")
    output.write(
        "Closest windows are not filtered by the mismatch threshold "
        "and do not model gaps.\n"
    )
    for region in scan_regions:
        region_matches = [
            match
            for match in result.closest_matches
            if match.scan_region == region.name
        ]
        output.write(f"\nRegion: {region.name}\n")
        if region_matches:
            output.write(transcript_match_terminal_table(region_matches))
            output.write("\n")
        else:
            output.write("No windows available for this region.\n")
    return output.getvalue().rstrip() + "\n"
