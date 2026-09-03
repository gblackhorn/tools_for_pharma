"""Command-line interface for local transcript scans and explicit remote BLAST."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import math
from pathlib import Path
import sys
from typing import Callable

from tools_for_pharma.oligo.ncbi_transport import (
    DEFAULT_DATABASE,
    DEFAULT_EMAIL,
    DEFAULT_EXPECT,
    DEFAULT_HITLIST_SIZE,
    DEFAULT_MEGABLAST_WORD_SIZE,
    DEFAULT_POLL_SECONDS,
    DEFAULT_REQUEST_SECONDS,
    DEFAULT_TOOL,
    DEFAULT_WORD_SIZE,
    NcbiBlastClient,
    resolve_blast_word_size,
)
from tools_for_pharma.oligo.transcript_accessions import (
    normalize_versioned_refseq_accession,
)
from tools_for_pharma.oligo.transcript_scan.models import (
    AntisenseQuery,
    AntisenseRegion,
    ComparisonResult,
    TranscriptMatch,
    TranscriptTargetResult,
)
from tools_for_pharma.oligo.transcript_scan.queries import (
    DEFAULT_BATCH_BASES,
    clean_text_for_id,
    default_query_name,
    normalize_sequence_type,
    parse_scan_regions,
    read_antisense_file,
)
from tools_for_pharma.oligo.transcript_scan.remote_blast import (
    BlastBatchResult,
    RemoteBlastConfig,
    combine_blast_csv,
    run_blast_batches as run_remote_blast_batches,
)
from tools_for_pharma.oligo.transcript_scan.reporting import (
    default_private_panel_workbook,
    default_result_workbook,
    format_closest_transcript_matches_for_terminal,
    format_transcript_matches_for_terminal,
    transcript_matches_to_csv,
    write_result_workbook,
)
from tools_for_pharma.oligo.transcript_scan.scanner import (
    DEFAULT_MAX_MISMATCHES,
    closest_transcript_matches,
)
from tools_for_pharma.oligo.transcript_scan.targets import (
    LocalFileTargetSource,
    PastedTargetSource,
    local_transcript_target,
    transcript_target_source,
)
from tools_for_pharma.oligo.transcript_scan.workflows import (
    LocalScanConfig,
    PrivatePanelWorkflowConfig,
    run_local_scan as run_local_scan_workflow,
    run_local_scan_with_comparison as run_local_comparison_workflow,
    run_private_panel_workflow as run_private_panel_domain_workflow,
)
from tools_for_pharma.sequence.nucleotides import normalize_rna


DEFAULT_CLOSEST_MATCHES = 10


def read_antisense_table(
    path: Path,
    sequence_column: str | None = None,
    name_column: str | None = None,
    target_accession_column: str | None = None,
    sheet_name: str | int | None = None,
    sequence_type: str = "AS",
) -> list[AntisenseQuery]:
    """Read AS or SS queries from an Excel or CSV table."""
    import pandas as pd

    normalized_type = normalize_sequence_type(sequence_type)
    suffix = path.suffix.lower()
    if suffix in {".xlsx", ".xlsm", ".xls"}:
        table = pd.read_excel(path, sheet_name=sheet_name or 0)
    elif suffix in {".csv", ".txt"}:
        table = pd.read_csv(path)
    else:
        raise ValueError(
            f"{normalized_type} table must be an Excel workbook or CSV/text file."
        )
    if table.empty:
        raise ValueError(f"{normalized_type} table is empty: {path}")

    columns_by_lower = {
        clean_text_for_id(column).lower(): str(column) for column in table.columns
    }
    if sequence_column is None:
        candidates = (
            ["sense", "ss", "ss_sequence", "ss sequence", "sequence"]
            if normalized_type == "SS"
            else ["antisense", "as", "as_sequence", "as sequence", "sequence"]
        )
        sequence_column = next(
            (columns_by_lower[item] for item in candidates if item in columns_by_lower),
            str(table.columns[0]),
        )
    if sequence_column not in table.columns:
        raise ValueError(
            f"{normalized_type} table is missing sequence column: {sequence_column}"
        )

    if name_column is None:
        name_column = next(
            (
                columns_by_lower[item]
                for item in [
                    "name",
                    "id",
                    "oligo",
                    "oligo id",
                    "oligo_id",
                    "as name",
                    "ss name",
                ]
                if item in columns_by_lower
            ),
            None,
        )
    elif name_column not in table.columns:
        raise ValueError(f"{normalized_type} table is missing name column: {name_column}")

    if target_accession_column is None:
        target_accession_column = next(
            (
                columns_by_lower[item]
                for item in [
                    "target_accession",
                    "target accession",
                    "refseq",
                    "refseq accession",
                ]
                if item in columns_by_lower
            ),
            None,
        )
    elif target_accession_column not in table.columns:
        raise ValueError(
            f"{normalized_type} table is missing target accession column: "
            f"{target_accession_column}"
        )

    metadata_columns = {
        "target_gene": next(
            (
                columns_by_lower[item]
                for item in ["target_gene", "target gene", "gene"]
                if item in columns_by_lower
            ),
            None,
        ),
        "species": columns_by_lower.get("species"),
        "notes": columns_by_lower.get("notes"),
    }

    records = []
    for row_index, row in table.iterrows():
        raw_sequence = row[sequence_column]
        if raw_sequence is None or pd.isna(raw_sequence):
            continue
        sequence = normalize_rna(str(raw_sequence))
        raw_name = row[name_column] if name_column else None
        name = (
            clean_text_for_id(raw_name)
            if raw_name is not None and not pd.isna(raw_name)
            else ""
        )
        raw_accession = row[target_accession_column] if target_accession_column else None
        target_accession = (
            clean_text_for_id(raw_accession)
            if raw_accession is not None and not pd.isna(raw_accession)
            else ""
        )
        metadata = {}
        for key, column in metadata_columns.items():
            raw_value = row[column] if column else None
            metadata[key] = (
                clean_text_for_id(raw_value)
                if raw_value is not None and not pd.isna(raw_value)
                else ""
            )
        source_fields = {
            str(column): (None if pd.isna(row[column]) else row[column])
            for column in table.columns
        }
        records.append(
            AntisenseQuery(
                name=name or default_query_name(normalized_type, row_index + 1),
                sequence_5to3=sequence,
                target_accession=target_accession,
                target_gene=metadata["target_gene"],
                species=metadata["species"],
                notes=metadata["notes"],
                sequence_type=normalized_type,
                source_fields=source_fields,
            )
        )
    if not records:
        raise ValueError(f"No {normalized_type} sequences found in table: {path}")
    return records


def read_antisense_queries(
    as_sequence: str | None = None,
    as_name: str | None = None,
    as_file: Path | None = None,
    as_table: Path | None = None,
    as_column: str | None = None,
    as_name_column: str | None = None,
    ss_sequence: str | None = None,
    ss_name: str | None = None,
    ss_file: Path | None = None,
    ss_table: Path | None = None,
    ss_column: str | None = None,
    ss_name_column: str | None = None,
    target_accession_column: str | None = None,
    as_sheet: str | int | None = None,
    ss_sheet: str | int | None = None,
) -> list[AntisenseQuery]:
    """Read AS or SS queries from exactly one input source."""
    provided = [
        as_sequence is not None,
        as_file is not None,
        as_table is not None,
        ss_sequence is not None,
        ss_file is not None,
        ss_table is not None,
    ]
    if sum(provided) != 1:
        raise ValueError(
            "Provide exactly one of --as-sequence, --as-file, --as-table, "
            "--ss-sequence, --ss-file, or --ss-table."
        )
    if as_sequence is not None:
        return [
            AntisenseQuery(
                as_name or default_query_name("AS"),
                normalize_rna(as_sequence),
                sequence_type="AS",
            )
        ]
    if as_file is not None:
        return read_antisense_file(as_file, sequence_type="AS")
    if as_table is not None:
        return read_antisense_table(
            as_table,
            as_column,
            as_name_column,
            target_accession_column,
            as_sheet,
            sequence_type="AS",
        )
    if ss_sequence is not None:
        return [
            AntisenseQuery(
                ss_name or default_query_name("SS"),
                normalize_rna(ss_sequence),
                sequence_type="SS",
            )
        ]
    if ss_file is not None:
        return read_antisense_file(ss_file, sequence_type="SS")
    assert ss_table is not None
    return read_antisense_table(
        ss_table,
        ss_column,
        ss_name_column,
        target_accession_column,
        ss_sheet,
        sequence_type="SS",
    )


def target_accession_values(value: object) -> list[str]:
    """Normalize a scalar or repeatable argparse accession value."""
    if value is None:
        return []
    values = [value] if isinstance(value, str) else list(value)
    return [clean_text_for_id(item) for item in values if clean_text_for_id(item)]


def read_target_accession_table(
    path: Path,
    accession_column: str | None = None,
    sheet_name: str | int | None = None,
) -> list[str]:
    """Read versioned target accessions from text, CSV, or Excel."""
    if not path.exists() or not path.is_file():
        raise ValueError(f"Target accession table does not exist: {path}")
    suffix = path.suffix.lower()
    if suffix in {".txt", ".list"}:
        values = [
            line.strip()
            for line in path.read_text(encoding="utf-8-sig").splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]
        if values and values[0].lower() in {
            "target_accession",
            "accession",
            "refseq",
            "refseq accession",
        }:
            values = values[1:]
        if not values:
            raise ValueError(f"No target accessions found in: {path}")
        return values

    import pandas as pd

    if suffix in {".xlsx", ".xlsm", ".xls"}:
        table = pd.read_excel(path, sheet_name=sheet_name or 0)
    elif suffix in {".csv", ".tsv"}:
        table = pd.read_csv(path, sep="\t" if suffix == ".tsv" else ",")
    else:
        raise ValueError(
            "--target-table must be a .txt, .list, .csv, .tsv, or Excel file."
        )
    if table.empty:
        raise ValueError(f"Target accession table is empty: {path}")
    columns_by_lower = {
        clean_text_for_id(column).lower(): str(column) for column in table.columns
    }
    if accession_column is None:
        accession_column = next(
            (
                columns_by_lower[item]
                for item in [
                    "target_accession",
                    "accession",
                    "refseq",
                    "refseq accession",
                ]
                if item in columns_by_lower
            ),
            str(table.columns[0]),
        )
    if accession_column not in table.columns:
        raise ValueError(f"Target table is missing accession column: {accession_column}")
    accessions = [
        clean_text_for_id(value)
        for value in table[accession_column]
        if value is not None and not pd.isna(value) and clean_text_for_id(value)
    ]
    if not accessions:
        raise ValueError(f"No target accessions found in: {path}")
    return accessions


def panel_accessions_from_args(args: argparse.Namespace) -> list[str]:
    """Return de-duplicated exact-version target accessions for panel mode."""
    raw_accessions = target_accession_values(args.target_accession)
    if getattr(args, "target_table", None):
        raw_accessions.extend(
            read_target_accession_table(
                args.target_table,
                getattr(args, "target_column", None),
                getattr(args, "target_sheet", None),
            )
        )
    if not raw_accessions:
        raise ValueError(
            "Private panel mode requires at least one --target-accession or --target-table."
        )
    accessions = []
    seen = set()
    for raw_accession in raw_accessions:
        accession = normalize_versioned_refseq_accession(raw_accession)
        if accession not in seen:
            seen.add(accession)
            accessions.append(accession)
    return accessions


def local_transcript_target_from_args(
    args: argparse.Namespace,
) -> TranscriptTargetResult:
    """Build a ready target record from a pasted sequence or one local file."""
    source = transcript_target_source(
        transcript_sequence=args.target_sequence,
        transcript_file=args.target_file,
    )
    if not isinstance(source, (PastedTargetSource, LocalFileTargetSource)):
        raise ValueError("Expected a pasted sequence or local transcript file.")
    return local_transcript_target(source)


def private_panel_requested(args: argparse.Namespace) -> bool:
    return bool(
        getattr(args, "private_panel", False)
        or getattr(args, "target_table", None)
        or getattr(args, "offline", False)
        or getattr(args, "refresh_targets", False)
        or getattr(args, "download_targets_only", False)
        or len(target_accession_values(args.target_accession)) > 1
    )


def private_panel_cache_dir(args: argparse.Namespace) -> Path:
    return args.cache_dir or Path(".ncbi_transcript_cache")


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compare AS or SS oligos to transcripts and optionally run NCBI BLAST."
    )
    parser.add_argument("--as-sequence", help="One AS oligo sequence in 5'->3'.")
    parser.add_argument("--as-name", help="Name for the single --as-sequence input.")
    parser.add_argument("--as-file", type=Path, help="Text or FASTA file containing AS sequences.")
    parser.add_argument("--as-table", type=Path, help="Excel/CSV table containing AS sequences.")
    parser.add_argument(
        "--as-column",
        help="AS sequence column for --as-table. Defaults to antisense/as/sequence or first column.",
    )
    parser.add_argument("--as-name-column", help="Optional AS name/id column for --as-table.")
    parser.add_argument("--as-sheet", help="Excel worksheet for --as-table. Defaults to first sheet.")
    parser.add_argument(
        "--ss-sequence",
        help="One SS/sense oligo sequence in 5'->3' transcript orientation.",
    )
    parser.add_argument("--ss-name", help="Name for the single --ss-sequence input.")
    parser.add_argument(
        "--ss-file", type=Path, help="Text or FASTA file containing SS/sense sequences."
    )
    parser.add_argument(
        "--ss-table", type=Path, help="Excel/CSV table containing SS/sense sequences."
    )
    parser.add_argument(
        "--ss-column",
        help="SS/sense sequence column for --ss-table. Defaults to sense/ss/sequence or first column.",
    )
    parser.add_argument("--ss-name-column", help="Optional SS/sense name/id column for --ss-table.")
    parser.add_argument("--ss-sheet", help="Excel worksheet for --ss-table. Defaults to first sheet.")
    parser.add_argument(
        "--target-accession",
        action="append",
        help=(
            "Versioned NM/XM/NR/XR accession to fetch from NCBI. Repeat this "
            "option to run a private local transcript-panel scan."
        ),
    )
    parser.add_argument(
        "--target-accession-column",
        help=(
            "Column in --as-table or --ss-table containing per-row NM/XM/NR/XR "
            "accessions. Defaults to target_accession if present."
        ),
    )
    parser.add_argument("--target-file", type=Path, help="FASTA/plain transcript file.")
    parser.add_argument("--target-sequence", help="Pasted FASTA/plain transcript sequence.")
    parser.add_argument(
        "--target-table",
        type=Path,
        help=(
            "Text/CSV/Excel list of versioned RefSeq transcript accessions for "
            "private local panel scanning."
        ),
    )
    parser.add_argument(
        "--target-column",
        help="Accession column for --target-table. Defaults to target_accession/accession/refseq.",
    )
    parser.add_argument(
        "--target-sheet",
        help="Excel worksheet for --target-table. Defaults to the first sheet.",
    )
    parser.add_argument(
        "--private-panel",
        action="store_true",
        help=(
            "Scan every input guide against every requested transcript locally. "
            "Guide sequences are never included in NCBI EFetch requests."
        ),
    )
    parser.add_argument(
        "--offline",
        action="store_true",
        help=(
            "Private panel mode only: do not contact NCBI and require every "
            "target transcript to exist in --cache-dir."
        ),
    )
    parser.add_argument(
        "--refresh-targets",
        action="store_true",
        help=(
            "Private panel mode only: download and revalidate requested transcript "
            "references even when exact versions already exist in --cache-dir."
        ),
    )
    parser.add_argument(
        "--download-targets-only",
        action="store_true",
        help=(
            "Retrieve and verify private-panel transcript references without "
            "reading or scanning guide sequences."
        ),
    )
    parser.add_argument(
        "--scan-region",
        action="append",
        help=(
            "Oligo region to scan locally: full, 2-18, or seed:2-8. "
            "Can be repeated. Defaults to full."
        ),
    )
    parser.add_argument(
        "--max-mismatches",
        type=int,
        default=DEFAULT_MAX_MISMATCHES,
        help=f"Maximum mismatches for local transcript scan. Defaults to {DEFAULT_MAX_MISMATCHES}.",
    )
    parser.add_argument(
        "--email",
        default=DEFAULT_EMAIL,
        help=(
            "Contact email for NCBI API usage guidelines. Required whenever "
            "a transcript or BLAST query must be submitted to NCBI."
        ),
    )
    parser.add_argument("--tool", default=DEFAULT_TOOL, help=f"NCBI tool name. Defaults to {DEFAULT_TOOL}.")
    parser.add_argument(
        "--blast",
        action="store_true",
        help=(
            "Also submit the input oligo sequence(s) to the remote NCBI BLAST "
            "URL API. This transmits query sequences outside the local computer."
        ),
    )
    parser.add_argument(
        "--blast-only",
        action="store_true",
        help=(
            "Submit input oligo sequence(s) to remote NCBI BLAST without a "
            "specific target transcript. This transmits query sequences outside "
            "the local computer."
        ),
    )
    parser.add_argument("--database", default=DEFAULT_DATABASE, help=f"BLAST database. Defaults to {DEFAULT_DATABASE}.")
    parser.add_argument("--expect", default=DEFAULT_EXPECT, help=f"BLAST expect value. Defaults to {DEFAULT_EXPECT}.")
    parser.add_argument("--word-size", type=int, default=DEFAULT_WORD_SIZE, help=f"BLAST word size. Defaults to {DEFAULT_WORD_SIZE}.")
    parser.add_argument(
        "--hitlist-size",
        type=int,
        default=DEFAULT_HITLIST_SIZE,
        help=f"Number of BLAST hits to retrieve. Defaults to {DEFAULT_HITLIST_SIZE}.",
    )
    parser.add_argument(
        "--megablast",
        action="store_true",
        help=(
            "Enable megablast for near-identical hits. When --word-size remains "
            f"{DEFAULT_WORD_SIZE}, it is changed to {DEFAULT_MEGABLAST_WORD_SIZE}."
        ),
    )
    parser.add_argument(
        "--timeout-seconds",
        type=int,
        default=1800,
        help="Maximum time to wait for BLAST completion. Defaults to 1800.",
    )
    parser.add_argument(
        "--max-batch-bases",
        type=int,
        default=DEFAULT_BATCH_BASES,
        help=(
            "Maximum total oligo bases per BLAST multi-FASTA submission. "
            f"Defaults to {DEFAULT_BATCH_BASES}."
        ),
    )
    parser.add_argument(
        "--request-seconds",
        type=int,
        default=DEFAULT_REQUEST_SECONDS,
        help=(
            "Minimum seconds between NCBI requests. Defaults to the conservative "
            f"{DEFAULT_REQUEST_SECONDS}."
        ),
    )
    parser.add_argument(
        "--poll-seconds",
        type=int,
        default=DEFAULT_POLL_SECONDS,
        help=(
            "Minimum seconds between status checks for the same BLAST RID. "
            f"Defaults to the conservative {DEFAULT_POLL_SECONDS}."
        ),
    )
    parser.add_argument(
        "--filter-max-mismatches",
        type=int,
        default=DEFAULT_MAX_MISMATCHES,
        help=(
            "Max mismatches for the blast_hits_filtered workbook sheet. "
            f"Defaults to {DEFAULT_MAX_MISMATCHES}."
        ),
    )
    parser.add_argument(
        "--filter-max-gap-opens",
        type=int,
        default=0,
        help="Max gap opens for the blast_hits_filtered workbook sheet. Defaults to 0.",
    )
    parser.add_argument(
        "--filter-min-alignment-fraction",
        type=float,
        default=0.8,
        help=(
            "Minimum alignment_length/query_length for filtered BLAST hits. "
            "Defaults to 0.8."
        ),
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        help=(
            "Optional folder for one-record NCBI EFetch transcript FASTA files. "
            "Private panel mode defaults to '.ncbi_transcript_cache'."
        ),
    )
    parser.add_argument(
        "--rid-log",
        type=Path,
        help="Optional CSV log written as soon as each BLAST RID is submitted.",
    )
    parser.add_argument("--output", type=Path, help="Write local scan CSV to this path.")
    parser.add_argument("--blast-output", type=Path, help="Write BLAST CSV to this path.")
    parser.add_argument(
        "--terminal",
        action="store_true",
        help=(
            "Print a readable local transcript scan summary to the terminal. "
            "For one --as-sequence or --ss-sequence scan without --output, this is the default."
        ),
    )
    parser.add_argument(
        "--stdout-csv",
        action="store_true",
        help="Print local transcript scan CSV to the terminal instead of the quick summary.",
    )
    parser.add_argument(
        "--closest",
        type=int,
        help=(
            "Print the N closest local transcript windows in terminal output, "
            "without applying --max-mismatches. In quick terminal mode, the tool "
            f"automatically shows the top {DEFAULT_CLOSEST_MATCHES} closest windows "
            "when no matches pass --max-mismatches. Private panel mode writes "
            "these to closest_transcript_windows."
        ),
    )
    parser.add_argument(
        "--result-workbook",
        type=Path,
        help=(
            "Write an Excel workbook with input_queries, local scan, BLAST hits, "
            "transcript-panel status, batch metadata, and run metadata as "
            "applicable. For --as-file/--as-table or --ss-file/--ss-table, defaults "
            "to a workflow-specific '<input>_*_results.xlsx' name when no CSV "
            "output is requested."
        ),
    )
    parser.add_argument(
        "--gui",
        action="store_true",
        help=(
            "Choose an AS/SS Excel input and single-transcript or private-panel "
            "target source with dialogs."
        ),
    )
    return parser


def validate_runtime_args(args: argparse.Namespace) -> None:
    """Validate CLI/GUI settings before local work or network requests begin."""
    if args.max_mismatches < 0:
        raise ValueError("--max-mismatches must be 0 or greater.")
    if args.filter_max_mismatches < 0:
        raise ValueError("--filter-max-mismatches must be 0 or greater.")
    if args.filter_max_gap_opens < 0:
        raise ValueError("--filter-max-gap-opens must be 0 or greater.")
    if not math.isfinite(args.filter_min_alignment_fraction) or not (
        0 <= args.filter_min_alignment_fraction <= 1
    ):
        raise ValueError("--filter-min-alignment-fraction must be between 0 and 1.")
    if args.hitlist_size < 1:
        raise ValueError("--hitlist-size must be 1 or greater.")
    if args.max_batch_bases < 1:
        raise ValueError("--max-batch-bases must be 1 or greater.")
    if args.timeout_seconds < 1:
        raise ValueError("--timeout-seconds must be 1 or greater.")
    if args.request_seconds < 0:
        raise ValueError("--request-seconds must be 0 or greater.")
    if args.poll_seconds < 0:
        raise ValueError("--poll-seconds must be 0 or greater.")
    if args.closest is not None and args.closest < 1:
        raise ValueError("--closest must be 1 or greater.")
    try:
        expect = float(args.expect)
    except (TypeError, ValueError) as error:
        raise ValueError("--expect must be a number greater than 0.") from error
    if not math.isfinite(expect) or expect <= 0:
        raise ValueError("--expect must be a number greater than 0.")
    if (args.blast or args.blast_only) and not clean_text_for_id(args.database):
        raise ValueError("--database cannot be blank for remote BLAST.")
    if getattr(args, "target_column", None) and not getattr(args, "target_table", None):
        raise ValueError("--target-column requires --target-table.")
    if getattr(args, "target_sheet", None) and not getattr(args, "target_table", None):
        raise ValueError("--target-sheet requires --target-table.")
    if getattr(args, "offline", False) and getattr(args, "refresh_targets", False):
        raise ValueError("--offline cannot be combined with --refresh-targets.")
    if private_panel_requested(args):
        args.private_panel = True
        if args.blast or args.blast_only:
            raise ValueError(
                "Private panel mode cannot be combined with --blast or --blast-only; "
                "guide sequences must remain local."
            )
        if args.target_accession_column:
            raise ValueError(
                "Private panel mode cannot use the per-query --target-accession-column."
            )
        if args.target_file or args.target_sequence:
            raise ValueError(
                "Private panel mode uses versioned --target-accession values or "
                "--target-table, not --target-file/--target-sequence."
            )
    args.word_size = resolve_blast_word_size(args.word_size, args.megablast)


def args_antisense_queries(args: argparse.Namespace) -> list[AntisenseQuery]:
    return read_antisense_queries(
        as_sequence=args.as_sequence,
        as_name=args.as_name,
        as_file=args.as_file,
        as_table=args.as_table,
        as_column=args.as_column,
        as_name_column=args.as_name_column,
        ss_sequence=args.ss_sequence,
        ss_name=args.ss_name,
        ss_file=args.ss_file,
        ss_table=args.ss_table,
        ss_column=args.ss_column,
        ss_name_column=args.ss_name_column,
        target_accession_column=args.target_accession_column,
        as_sheet=args.as_sheet,
        ss_sheet=args.ss_sheet,
    )


def local_scan_config_from_args(args: argparse.Namespace) -> LocalScanConfig:
    """Adapt legacy arguments to the explicit local workflow configuration."""
    return LocalScanConfig(
        target_accessions=tuple(target_accession_values(args.target_accession)),
        use_query_target_accession=bool(args.target_accession_column),
        target_sequence=args.target_sequence,
        target_file=args.target_file,
        email=args.email,
        tool=args.tool,
        cache_dir=args.cache_dir,
        max_mismatches=args.max_mismatches,
    )


def remote_blast_config_from_args(args: argparse.Namespace) -> RemoteBlastConfig:
    """Adapt legacy arguments to an explicitly remote submission configuration."""
    return RemoteBlastConfig(
        email=args.email,
        tool=args.tool,
        database=args.database,
        expect=args.expect,
        word_size=args.word_size,
        hitlist_size=args.hitlist_size,
        megablast=args.megablast,
        timeout_seconds=args.timeout_seconds,
        max_batch_bases=args.max_batch_bases,
        request_seconds=args.request_seconds,
        poll_seconds=args.poll_seconds,
        rid_log=args.rid_log,
    )


def private_panel_config_from_args(
    args: argparse.Namespace,
) -> PrivatePanelWorkflowConfig:
    """Adapt legacy arguments to the explicit local panel configuration."""
    return PrivatePanelWorkflowConfig(
        accessions=tuple(panel_accessions_from_args(args)),
        cache_dir=private_panel_cache_dir(args),
        email=args.email,
        tool=args.tool,
        offline=args.offline,
        refresh=bool(getattr(args, "refresh_targets", False)),
        request_seconds=args.request_seconds,
        max_mismatches=args.max_mismatches,
        closest=args.closest,
        download_targets_only=args.download_targets_only,
    )


def run_local_scan(
    args: argparse.Namespace,
    queries: list[AntisenseQuery],
    scan_regions: list[AntisenseRegion],
    max_mismatches: int | None = DEFAULT_MAX_MISMATCHES,
) -> list[TranscriptMatch]:
    return run_local_scan_workflow(
        local_scan_config_from_args(args),
        queries,
        scan_regions,
        max_mismatches=max_mismatches,
    )


def run_local_scan_with_comparison(
    args: argparse.Namespace,
    queries: list[AntisenseQuery],
    scan_regions: list[AntisenseRegion],
) -> tuple[list[TranscriptMatch], list[ComparisonResult]]:
    return run_local_comparison_workflow(
        local_scan_config_from_args(args),
        queries,
        scan_regions,
    )


def run_blast_batches(
    args: argparse.Namespace,
    queries: list[AntisenseQuery],
    *,
    client_factory: Callable[..., NcbiBlastClient] | None = None,
    sleeper: Callable[[float], None] | None = None,
) -> list[BlastBatchResult]:
    return run_remote_blast_batches(
        remote_blast_config_from_args(args),
        queries,
        client_factory=client_factory,
        sleeper=sleeper,
    )


def run_private_panel_workflow(
    args: argparse.Namespace,
    started_at: str,
    *,
    progress_callback: Callable[[int, int, str, str], None] | None = None,
    cancel_check: Callable[[], bool] | None = None,
    include_comparison_results: bool = False,
) -> int:
    """Retrieve public references and scan private guides entirely locally."""
    print(
        "Private local panel mode: guide sequences remain on this computer; "
        "NCBI EFetch requests contain transcript accessions only.",
        file=sys.stderr,
    )
    queries = [] if args.download_targets_only else args_antisense_queries(args)
    scan_regions = parse_scan_regions(args.scan_region)
    panel_result = run_private_panel_domain_workflow(
        private_panel_config_from_args(args),
        queries,
        scan_regions,
        progress_callback=progress_callback,
        cancel_check=cancel_check,
    )
    local_matches = list(panel_result.matches)
    if args.output:
        write_text(args.output, transcript_matches_to_csv(local_matches))
        print(f"Wrote private local panel matches to: {args.output}")
    if args.stdout_csv:
        print(transcript_matches_to_csv(local_matches), end="")

    result_workbook = default_private_panel_workbook(args)
    write_result_workbook(
        result_workbook,
        args,
        queries,
        scan_regions,
        local_matches,
        [],
        started_at,
        datetime.now(timezone.utc).isoformat(),
        include_blast_sheets=False,
        comparison_results=(
            list(panel_result.comparison_results)
            if include_comparison_results and not args.download_targets_only
            else None
        ),
        transcript_targets=list(panel_result.targets),
        query_target_summaries=(
            None
            if args.download_targets_only or include_comparison_results
            else list(panel_result.summaries)
        ),
        closest_local_matches=(
            list(panel_result.closest_matches) if args.closest is not None else None
        ),
    )
    ready_count = sum(target.status == "ready" for target in panel_result.targets)
    error_count = len(panel_result.targets) - ready_count
    print(
        f"Private panel targets: {ready_count} ready, {error_count} error; "
        f"queries: {len(queries)}; local matches: {len(local_matches)}"
    )
    print(f"Wrote private transcript panel workbook to: {result_workbook}")
    if error_count:
        print(
            "One or more target references could not be used; see transcript_targets.",
            file=sys.stderr,
        )
    if not ready_count:
        return 1
    return 2 if error_count else 0


def run_cli_args(args: argparse.Namespace, *, started_at: str) -> int:
    """Dispatch one already-parsed and validated non-GUI command."""
    if private_panel_requested(args):
        return run_private_panel_workflow(args, started_at)
    queries = args_antisense_queries(args)
    if (
        not args.blast_only
        and not args.target_accession_column
        and not args.target_accession
        and not args.target_file
        and not args.target_sequence
        and any(query.target_accession for query in queries)
    ):
        args.target_accession_column = "target_accession"
    scan_regions = parse_scan_regions(args.scan_region)
    local_matches: list[TranscriptMatch] = []
    closest_matches: list[TranscriptMatch] = []
    blast_outputs: list[BlastBatchResult] = []
    if args.blast_only:
        args.blast = True
    if not args.blast_only:
        result_workbook = default_result_workbook(args)
        quick_terminal_default = (
            (args.as_sequence or args.ss_sequence)
            and len(queries) == 1
            and not args.output
            and not result_workbook
            and not args.stdout_csv
        )
        print_terminal = args.terminal or args.closest is not None or quick_terminal_default
        prepare_closest_matches = args.closest is not None or print_terminal
        if prepare_closest_matches:
            all_local_matches = run_local_scan(
                args,
                queries,
                scan_regions,
                max_mismatches=None,
            )
            local_matches = [
                match for match in all_local_matches if match.mismatches <= args.max_mismatches
            ]
            closest_limit = args.closest or DEFAULT_CLOSEST_MATCHES
            closest_matches = closest_transcript_matches(all_local_matches, closest_limit)
        else:
            local_matches = run_local_scan(args, queries, scan_regions)
        csv_text = transcript_matches_to_csv(local_matches)
        show_closest_matches = args.closest is not None or (
            print_terminal and not local_matches
        )
        closest_limit = args.closest or DEFAULT_CLOSEST_MATCHES
        if args.output:
            write_text(args.output, csv_text)
            print(f"Wrote local transcript scan to: {args.output}")
        if print_terminal:
            print(
                format_transcript_matches_for_terminal(
                    local_matches,
                    queries,
                    scan_regions,
                    args.max_mismatches,
                )
            )
            if show_closest_matches:
                print(
                    format_closest_transcript_matches_for_terminal(
                        closest_matches,
                        closest_limit,
                        args.max_mismatches,
                    )
                )
        elif not args.output and not result_workbook:
            print(csv_text, end="")
            if show_closest_matches:
                print(
                    format_closest_transcript_matches_for_terminal(
                        closest_matches,
                        closest_limit,
                        args.max_mismatches,
                    )
                )
        if not local_matches and not print_terminal:
            print(
                "No local transcript matches found within mismatch threshold.",
                file=sys.stderr,
            )

    if args.blast:
        blast_outputs = run_blast_batches(args, queries)
        blast_csv = combine_blast_csv(blast_outputs)
        if args.blast_output:
            write_text(args.blast_output, blast_csv)
            rids = ", ".join(result.submission.rid for result in blast_outputs)
            print(f"Wrote BLAST CSV for RID(s) {rids} to: {args.blast_output}")
        elif not default_result_workbook(args):
            print(blast_csv, end="")

    result_workbook = default_result_workbook(args)
    if result_workbook:
        write_result_workbook(
            result_workbook,
            args,
            queries,
            scan_regions,
            local_matches,
            blast_outputs,
            started_at,
            datetime.now(timezone.utc).isoformat(),
            include_blast_sheets=args.blast,
        )
        workflow_label = "NCBI BLAST" if args.blast else "local transcript scan"
        print(f"Wrote {workflow_label} result workbook to: {result_workbook}")
    return 0


def main(
    argv: list[str] | None = None,
    *,
    gui_runner: Callable[[], int] | None = None,
) -> int:
    """Parse and dispatch the CLI; the compatibility facade injects its GUI."""
    args = build_parser().parse_args(argv)
    if args.gui:
        if gui_runner is None:
            raise RuntimeError("A GUI runner is required for --gui dispatch.")
        return gui_runner()
    started_at = datetime.now(timezone.utc).isoformat()
    try:
        validate_runtime_args(args)
        return run_cli_args(args, started_at=started_at)
    except (ValueError, TimeoutError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1
