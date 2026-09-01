"""NCBI transcript fetch and BLAST helpers for AS/SS oligos.

This module has two related workflows:

1. Specific transcript check:
   Fetch an NM/XM/NR/XR accession with NCBI EFetch, then scan the transcript for
   the reverse-complement target of an antisense sequence or the direct target
   of a sense sequence.

2. BLAST database search:
   Submit the oligo sequence to the NCBI BLAST URL API and retrieve a CSV
   report. This is best for broad searches such as refseq_rna or core_nt.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field, replace
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
import logging
import math
import os
from pathlib import Path
import re
import sys
import time
from typing import Callable, Iterable

from tools_for_pharma.oligo.core import get_complementary_sequence
from tools_for_pharma.oligo.ncbi_transport import (
    BLAST_URL,
    BLASTN_WORD_SIZES,
    DEFAULT_DATABASE,
    DEFAULT_EMAIL,
    DEFAULT_EXPECT,
    DEFAULT_HITLIST_SIZE,
    DEFAULT_MEGABLAST_WORD_SIZE,
    DEFAULT_POLL_SECONDS,
    DEFAULT_PROGRAM,
    DEFAULT_REQUEST_SECONDS,
    DEFAULT_TOOL,
    DEFAULT_WORD_SIZE,
    EFETCH_URL,
    MEGABLAST_WORD_SIZES,
    BlastSubmission,
    NcbiBlastClient,
    NcbiHttpClient,
    efetch_fasta_params,
    parse_blast_field,
    require_email,
    resolve_blast_word_size,
)
from tools_for_pharma.oligo.transcript import fasta_or_plain_text_to_sequence, get_fasta_header
from tools_for_pharma.oligo.transcript_accessions import (
    VERSIONED_REFSEQ_GENOMIC_RE,
    VERSIONED_REFSEQ_TRANSCRIPT_RE,
    extract_refseq_accession_from_header,
    normalize_versioned_refseq_accession,
)
from tools_for_pharma.sequence.comparison import mismatch_positions_1based
from tools_for_pharma.sequence.fasta import (
    FastaRecord,
    format_fasta,
    parse_fasta,
)
from tools_for_pharma.sequence.nucleotides import (
    normalize_dna as normalize_dna_sequence,
    normalize_rna,
)
from tools_for_pharma.shared.excel_utils import list_excel_sheets


DEFAULT_MAX_MISMATCHES = 3
DEFAULT_CLOSEST_MATCHES = 10
DEFAULT_SINGLE_GUI_CLOSEST_MATCHES = 5
DEFAULT_BATCH_BASES = 1000
APP_DATA_DIR_NAME = "TranscriptScanData"
GUI_SETTINGS_FILE_NAME = "settings.json"
GUI_LOG_FILE_NAME = "transcript_scan.log"
CSV_COLUMNS = [
    "query_id",
    "subject_id",
    "percent_identity",
    "alignment_length",
    "mismatches",
    "gap_opens",
    "query_start",
    "query_end",
    "subject_start",
    "subject_end",
    "evalue",
    "bit_score",
]


@dataclass(frozen=True)
class AntisenseQuery:
    """One named AS or SS input sequence."""

    name: str
    sequence_5to3: str
    target_accession: str = ""
    target_gene: str = ""
    species: str = ""
    notes: str = ""
    sequence_type: str = "AS"
    blast_query_id: str = field(default="", compare=False)
    source_fields: dict[str, object] = field(default_factory=dict, compare=False, repr=False)


@dataclass(frozen=True)
class AntisenseRegion:
    """A 1-based inclusive AS subregion to scan."""

    name: str
    start: int | None = None
    end: int | None = None


@dataclass(frozen=True)
class TranscriptMatch:
    """One local AS/SS-vs-transcript match."""

    transcript_name: str
    antisense_name: str
    scan_region: str
    as_region_start: int
    as_region_end: int
    antisense_5to3: str
    antisense_region_5to3: str
    target_5to3: str
    transcript_start: int
    transcript_end: int
    mismatches: int
    transcript_window_5to3: str
    transcript_match_as_5to3: str
    mismatch_positions_1based: tuple[int, ...]
    as_mismatch_positions_1based: tuple[int, ...]
    sequence_type: str = "AS"


@dataclass(frozen=True)
class TranscriptTargetResult:
    """Retrieval and validation status for one private-panel transcript."""

    requested_accession: str
    retrieved_accession: str = ""
    transcript_name: str = ""
    sequence_5to3: str = field(default="", repr=False)
    sequence_length_nt: int = 0
    cache_path: str = ""
    cache_status: str = ""
    exact_version_match: bool = False
    sequence_sha256: str = ""
    retrieved_at_utc: str = ""
    status: str = "error"
    error: str = ""


@dataclass(frozen=True)
class QueryTargetSummary:
    """One status row for a guide-versus-transcript panel comparison."""

    query_name: str
    sequence_type: str
    requested_accession: str
    retrieved_accession: str
    target_status: str
    scan_status: str
    scan_regions: str
    match_count: int
    exact_match_count: int
    best_mismatches: int | None
    error: str = ""


@dataclass(frozen=True)
class ComparisonResult:
    """One compact best-result row for a query, transcript, and scan region."""

    input_order: int
    query_name: str
    target_accession: str
    scan_region: str
    region_start: int
    region_end: int
    result: str
    sites_within_threshold: int
    best_mismatches: int | None
    mismatch_positions_in_query_1based: tuple[int, ...] = ()
    best_transcript_start: int | None = None
    best_transcript_end: int | None = None
    query_region_5to3: str = ""
    best_match_in_query_orientation_5to3: str = ""
    differences: str = ""


@dataclass(frozen=True)
class PrivatePanelScanResult:
    """Complete result of a private local transcript-panel scan."""

    targets: tuple[TranscriptTargetResult, ...]
    matches: tuple[TranscriptMatch, ...]
    summaries: tuple[QueryTargetSummary, ...]
    closest_matches: tuple[TranscriptMatch, ...] = ()
    comparison_results: tuple[ComparisonResult, ...] = ()


@dataclass(frozen=True)
class BlastBatchResult:
    """One completed BLAST batch and its returned CSV text."""

    batch_index: int
    submission: BlastSubmission
    queries: tuple[AntisenseQuery, ...]
    csv_text: str


def normalize_dna(sequence: str) -> str:
    """Normalize a sequence to DNA letters for NCBI BLAST requests."""
    return normalize_dna_sequence(sequence)


def fasta_record(name: str, sequence: str, line_width: int = 80) -> str:
    """Return a simple FASTA record from a raw nucleotide sequence."""
    cleaned = normalize_dna(sequence)
    return format_fasta(
        FastaRecord(sanitize_fasta_name(name), "", cleaned),
        width=line_width,
        trailing_newline=False,
    )


def sanitize_fasta_name(name: str) -> str:
    """Return a FASTA-safe query identifier."""
    cleaned = re.sub(r"[^A-Za-z0-9_.:-]+", "_", clean_text_for_id(name)).strip("_")
    return cleaned or "oligo_query"


def assign_unique_blast_query_ids(records: Iterable[AntisenseQuery]) -> list[AntisenseQuery]:
    """Return records with stable, unique FASTA identifiers.

    Distinct input names can sanitize to the same FASTA identifier. Preserve the
    first identifier unchanged and suffix later collisions in input order.
    """
    assigned = []
    used: set[str] = set()
    next_suffix: dict[str, int] = {}
    for record in records:
        base = sanitize_fasta_name(record.blast_query_id or record.name)
        candidate = base
        suffix = next_suffix.get(base, 2)
        while candidate in used:
            candidate = f"{base}_{suffix}"
            suffix += 1
        next_suffix[base] = suffix
        used.add(candidate)
        assigned.append(replace(record, blast_query_id=candidate))
    return assigned


def normalize_sequence_type(value: str) -> str:
    cleaned = clean_text_for_id(value).upper()
    if cleaned not in {"AS", "SS"}:
        raise ValueError("Sequence type must be AS or SS.")
    return cleaned


def default_query_name(sequence_type: str, index: int | None = None) -> str:
    prefix = normalize_sequence_type(sequence_type)
    if index is None:
        return "antisense_query" if prefix == "AS" else "sense_query"
    return f"{prefix}_{index}"


def clean_text_for_id(value: object) -> str:
    return re.sub(r"\s+", " ", str(value)).strip()


def multi_fasta(records: Iterable[AntisenseQuery]) -> str:
    """Return a multi-FASTA string for one or more oligo queries."""
    prepared = assign_unique_blast_query_ids(records)
    return "\n".join(
        fasta_record(record.blast_query_id, record.sequence_5to3)
        for record in prepared
    )


def parse_fasta_records(text: str, sequence_type: str = "AS") -> list[AntisenseQuery]:
    """Parse FASTA text into AS or SS query records."""
    normalized_type = normalize_sequence_type(sequence_type)
    lines = str(text).splitlines()
    first_header_index = next(
        (
            index
            for index, raw_line in enumerate(lines)
            if raw_line.strip().startswith(">")
        ),
        None,
    )
    if first_header_index is None:
        return []

    compatibility_lines = lines[first_header_index:]
    query_names = []
    header_count = 0
    for index, raw_line in enumerate(compatibility_lines):
        line = raw_line.strip()
        if not line.startswith(">"):
            continue
        header_count += 1
        query_name = line[1:].strip() or default_query_name(
            normalized_type,
            header_count,
        )
        query_names.append(query_name)
        if not line[1:].strip():
            compatibility_lines[index] = f">{query_name}"

    fasta_records = parse_fasta(
        "\n".join(compatibility_lines),
        ignore_comments=False,
    )
    return [
        AntisenseQuery(
            query_name,
            normalize_rna(record.sequence),
            sequence_type=normalized_type,
        )
        for query_name, record in zip(query_names, fasta_records)
    ]


def parse_plain_antisense_lines(text: str, sequence_type: str = "AS") -> list[AntisenseQuery]:
    """Parse a plain text list of AS or SS sequences.

    Accepted line styles:
      AUGCUA...
      AS_001,AUGCUA...
      AS_001<TAB>AUGCUA...
      AS_001 AUGCUA...
    """
    records = []
    normalized_type = normalize_sequence_type(sequence_type)
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if "," in line:
            parts = [part.strip() for part in line.split(",", 1)]
        elif "\t" in line:
            parts = [part.strip() for part in line.split("\t", 1)]
        else:
            parts = line.split(maxsplit=1)

        if len(parts) == 2:
            name, sequence = parts
        else:
            name = default_query_name(normalized_type, len(records) + 1)
            sequence = parts[0]
        records.append(AntisenseQuery(name, normalize_rna(sequence), sequence_type=normalized_type))
    return records


def read_antisense_file(path: Path, sequence_type: str = "AS") -> list[AntisenseQuery]:
    """Read AS or SS queries from FASTA or plain text."""
    text = path.read_text(encoding="utf-8-sig")
    if any(line.lstrip().startswith(">") for line in text.splitlines()):
        records = parse_fasta_records(text, sequence_type=sequence_type)
    else:
        records = parse_plain_antisense_lines(text, sequence_type=sequence_type)
    if not records:
        raise ValueError(f"No {normalize_sequence_type(sequence_type)} sequences found in {path}.")
    return records


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
        raise ValueError(f"{normalized_type} table must be an Excel workbook or CSV/text file.")

    if table.empty:
        raise ValueError(f"{normalized_type} table is empty: {path}")

    columns_by_lower = {clean_text_for_id(column).lower(): str(column) for column in table.columns}
    if sequence_column is None:
        if normalized_type == "SS":
            candidates = ["sense", "ss", "ss_sequence", "ss sequence", "sequence"]
        else:
            candidates = ["antisense", "as", "as_sequence", "as sequence", "sequence"]
        for candidate in candidates:
            if candidate in columns_by_lower:
                sequence_column = columns_by_lower[candidate]
                break
        if sequence_column is None:
            sequence_column = str(table.columns[0])
    if sequence_column not in table.columns:
        raise ValueError(f"{normalized_type} table is missing sequence column: {sequence_column}")

    if name_column is None:
        for candidate in ["name", "id", "oligo", "oligo id", "oligo_id", "as name", "ss name"]:
            if candidate in columns_by_lower:
                name_column = columns_by_lower[candidate]
                break
    elif name_column not in table.columns:
        raise ValueError(f"{normalized_type} table is missing name column: {name_column}")

    if target_accession_column is None:
        for candidate in ["target_accession", "target accession", "refseq", "refseq accession"]:
            if candidate in columns_by_lower:
                target_accession_column = columns_by_lower[candidate]
                break
    elif target_accession_column not in table.columns:
        raise ValueError(f"{normalized_type} table is missing target accession column: {target_accession_column}")

    metadata_columns = {
        "target_gene": next(
            (columns_by_lower[candidate] for candidate in ["target_gene", "target gene", "gene"] if candidate in columns_by_lower),
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
        if name_column:
            raw_name = row[name_column]
            name = clean_text_for_id(raw_name) if raw_name is not None and not pd.isna(raw_name) else ""
        else:
            name = ""
        target_accession = ""
        if target_accession_column:
            raw_accession = row[target_accession_column]
            if raw_accession is not None and not pd.isna(raw_accession):
                target_accession = clean_text_for_id(raw_accession)

        metadata = {}
        for key, column in metadata_columns.items():
            value = ""
            if column:
                raw_value = row[column]
                if raw_value is not None and not pd.isna(raw_value):
                    value = clean_text_for_id(raw_value)
            metadata[key] = value

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


def duplicate_sequence_groups(records: list[AntisenseQuery]) -> dict[str, list[str]]:
    """Return normalized AS sequence to AS names for duplicate sequences only."""
    groups: dict[str, list[str]] = {}
    for record in records:
        groups.setdefault(normalize_rna(record.sequence_5to3), []).append(record.name)
    return {
        sequence: names
        for sequence, names in groups.items()
        if len(names) > 1
    }


def input_query_rows(records: list[AntisenseQuery]) -> list[dict[str, object]]:
    """Return input query rows with duplicate annotations."""
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


def batch_antisense_queries(
    records: list[AntisenseQuery],
    max_batch_bases: int = DEFAULT_BATCH_BASES,
) -> list[list[AntisenseQuery]]:
    """Group short oligo queries into multi-FASTA BLAST batches."""
    if max_batch_bases < 1:
        raise ValueError("--max-batch-bases must be 1 or greater.")

    batches: list[list[AntisenseQuery]] = []
    current: list[AntisenseQuery] = []
    current_bases = 0
    for record in records:
        sequence_bases = len(normalize_rna(record.sequence_5to3))
        if current and current_bases + sequence_bases > max_batch_bases:
            batches.append(current)
            current = []
            current_bases = 0
        current.append(record)
        current_bases += sequence_bases
    if current:
        batches.append(current)
    return batches


def parse_scan_region(value: str) -> AntisenseRegion:
    """Parse scan region specs such as full, 2-18, or seed:2-8."""
    text = clean_text_for_id(value)
    if not text:
        raise ValueError("Scan region cannot be blank.")
    if text.lower() == "full":
        return AntisenseRegion("full")

    if ":" in text:
        name, range_text = [part.strip() for part in text.split(":", 1)]
    else:
        name = text
        range_text = text
    match = re.fullmatch(r"(\d+)\s*-\s*(\d+)", range_text)
    if not match:
        raise ValueError(
            f"Invalid scan region '{value}'. Use 'full', '2-18', or 'seed:2-8'."
        )

    start = int(match.group(1))
    end = int(match.group(2))
    if start < 1 or end < start:
        raise ValueError(f"Invalid scan region coordinates: {value}")
    return AntisenseRegion(name or f"{start}-{end}", start, end)


def parse_scan_regions(values: list[str] | None) -> list[AntisenseRegion]:
    if not values:
        return [AntisenseRegion("full")]
    return [parse_scan_region(value) for value in values]


def antisense_region_sequence(sequence: str, region: AntisenseRegion) -> tuple[str, int, int]:
    antisense = normalize_rna(sequence)
    if region.start is None or region.end is None:
        return antisense, 1, len(antisense)
    if region.end > len(antisense):
        raise ValueError(
            f"Scan region {region.name} ends at {region.end}, but AS sequence "
            f"is only {len(antisense)} nt."
        )
    return antisense[region.start - 1 : region.end], region.start, region.end


def transcript_cache_path(cache_dir: Path, accession: str) -> Path:
    return cache_dir / f"{sanitize_fasta_name(accession)}.fasta"


def target_accession_values(value: object) -> list[str]:
    """Normalize a scalar or repeatable argparse accession value."""
    if value is None:
        return []
    if isinstance(value, str):
        values = [value]
    else:
        values = list(value)
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

    columns_by_lower = {clean_text_for_id(column).lower(): str(column) for column in table.columns}
    if accession_column is None:
        for candidate in ["target_accession", "accession", "refseq", "refseq accession"]:
            if candidate in columns_by_lower:
                accession_column = columns_by_lower[candidate]
                break
        if accession_column is None:
            accession_column = str(table.columns[0])
    if accession_column not in table.columns:
        raise ValueError(f"Target table is missing accession column: {accession_column}")

    accessions = []
    for value in table[accession_column]:
        if value is not None and not pd.isna(value):
            cleaned = clean_text_for_id(value)
            if cleaned:
                accessions.append(cleaned)
    if not accessions:
        raise ValueError(f"No target accessions found in: {path}")
    return accessions


def panel_accessions_from_args(args: argparse.Namespace) -> list[str]:
    """Return de-duplicated, exact-version target accessions for panel mode."""
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


def format_cached_transcript_fasta(header: str, sequence: str, width: int = 80) -> str:
    """Format a verified transcript while preserving its descriptive header."""
    dna = normalize_dna(sequence)
    return format_fasta(
        FastaRecord.from_header(clean_text_for_id(header), dna),
        width=width,
    )


def transcript_target_from_fasta(
    requested_accession: str,
    fasta_text: str,
    cache_path: Path,
    cache_status: str,
    retrieved_at_utc: str,
) -> TranscriptTargetResult:
    """Validate one fetched/cached transcript and build its target record."""
    validate_single_transcript_record(fasta_text, f"Transcript {requested_accession}")
    header = get_fasta_header(fasta_text)
    if not header:
        raise ValueError(f"Transcript {requested_accession} FASTA header is missing.")
    retrieved_accession = extract_refseq_accession_from_header(header)
    if retrieved_accession != requested_accession:
        raise ValueError(
            f"Requested exact RefSeq version {requested_accession}, but retrieved "
            f"{retrieved_accession}."
        )
    sequence = fasta_or_plain_text_to_sequence(fasta_text)
    sequence_dna = normalize_dna(sequence)
    return TranscriptTargetResult(
        requested_accession=requested_accession,
        retrieved_accession=retrieved_accession,
        transcript_name=header,
        sequence_5to3=sequence,
        sequence_length_nt=len(sequence),
        cache_path=str(cache_path),
        cache_status=cache_status,
        exact_version_match=True,
        sequence_sha256=hashlib.sha256(sequence_dna.encode("ascii")).hexdigest(),
        retrieved_at_utc=retrieved_at_utc,
        status="ready",
    )


def retrieve_transcript_targets(
    accessions: list[str],
    *,
    email: str,
    tool: str = DEFAULT_TOOL,
    cache_dir: Path,
    offline: bool = False,
    refresh: bool = False,
    request_seconds: int = DEFAULT_REQUEST_SECONDS,
    client: NcbiHttpClient | None = None,
    progress_callback: Callable[[int, int, str, str], None] | None = None,
    cancel_check: Callable[[], bool] | None = None,
) -> list[TranscriptTargetResult]:
    """Retrieve public transcript references without transmitting guide sequences."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    http_client = client

    results = []
    total = len(accessions)
    for index, accession in enumerate(accessions, start=1):
        if progress_callback:
            progress_callback(index - 1, total, accession, "starting")
        if cancel_check and cancel_check():
            for cancelled_index, cancelled_accession in enumerate(
                accessions[index - 1 :],
                start=index,
            ):
                cancelled_path = transcript_cache_path(cache_dir, cancelled_accession)
                results.append(
                    TranscriptTargetResult(
                        requested_accession=cancelled_accession,
                        cache_path=str(cancelled_path),
                        cache_status=(
                            "cache" if cancelled_path.exists() else "missing"
                        ),
                        status="error",
                        error="Transcript retrieval cancelled by user.",
                    )
                )
                if progress_callback:
                    progress_callback(
                        cancelled_index,
                        total,
                        cancelled_accession,
                        "cancelled",
                    )
            break
        cache_path = transcript_cache_path(cache_dir, accession)
        cache_existed = cache_path.exists()
        try:
            if cache_existed and not refresh:
                fasta_text = cache_path.read_text(encoding="utf-8-sig")
                retrieved_at = datetime.fromtimestamp(
                    cache_path.stat().st_mtime,
                    tz=timezone.utc,
                ).isoformat()
                target = transcript_target_from_fasta(
                    accession,
                    fasta_text,
                    cache_path,
                    "cache",
                    retrieved_at,
                )
            elif offline:
                raise ValueError(
                    f"Offline mode requires cached transcript {accession} at {cache_path}."
                )
            else:
                if http_client is None:
                    http_client = NcbiHttpClient(
                        email=require_email(email),
                        tool=tool,
                        request_seconds=max(request_seconds, DEFAULT_REQUEST_SECONDS),
                    )
                request_email = require_email(
                    getattr(http_client, "email", None) or email
                )
                fasta_text = http_client.get_text(
                    EFETCH_URL,
                    efetch_fasta_params(
                        accession,
                        email=request_email,
                        tool=tool,
                    ),
                )
                retrieved_at = datetime.now(timezone.utc).isoformat()
                target = transcript_target_from_fasta(
                    accession,
                    fasta_text,
                    cache_path,
                    "refreshed" if cache_existed else "downloaded",
                    retrieved_at,
                )
                cache_path.write_text(
                    format_cached_transcript_fasta(target.transcript_name, target.sequence_5to3),
                    encoding="utf-8",
                )
            results.append(target)
        except Exception as error:
            results.append(
                TranscriptTargetResult(
                    requested_accession=accession,
                    cache_path=str(cache_path),
                    cache_status=(
                        "refresh_failed"
                        if refresh and cache_path.exists()
                        else "missing" if not cache_path.exists() else "invalid"
                    ),
                    status="error",
                    error=str(error),
                )
            )
        if progress_callback:
            progress_status = (
                results[-1].cache_status
                if results[-1].status == "ready"
                else results[-1].status
            )
            progress_callback(index, total, accession, progress_status)
    return results


def validate_single_transcript_record(text: str, source_label: str) -> None:
    """Reject multi-record FASTA instead of concatenating transcript records."""
    record_count = sum(
        1
        for raw_line in str(text).splitlines()
        if raw_line.lstrip().startswith(">")
    )
    if record_count > 1:
        raise ValueError(
            f"{source_label} contains {record_count} FASTA records. "
            "The current local transcript scanner accepts exactly one transcript "
            "record per target; use separate one-record FASTA files."
        )


def fetch_transcript_fasta(
    accession: str,
    email: str | None,
    tool: str = DEFAULT_TOOL,
    cache_dir: Path | None = None,
    *,
    client: NcbiHttpClient | None = None,
) -> str:
    """Fetch transcript FASTA from NCBI EFetch by accession or UID."""
    if cache_dir is not None:
        cache_path = transcript_cache_path(cache_dir, accession)
        if cache_path.exists():
            return cache_path.read_text(encoding="utf-8-sig")

    request_email = require_email(getattr(client, "email", None) or email)
    http_client = client or NcbiHttpClient(email=request_email, tool=tool)
    text = http_client.get_text(
        EFETCH_URL,
        efetch_fasta_params(
            accession,
            email=getattr(http_client, "email", None) or request_email,
            tool=tool,
        ),
    )
    if not text.lstrip().startswith(">"):
        raise ValueError(f"NCBI EFetch did not return FASTA for {accession}:\n{text[:500]}")
    if cache_dir is not None:
        cache_dir.mkdir(parents=True, exist_ok=True)
        transcript_cache_path(cache_dir, accession).write_text(text, encoding="utf-8")
    return text


def read_transcript_input(
    transcript_sequence: str | None = None,
    transcript_file: Path | None = None,
    accession: str | None = None,
    email: str | None = None,
    tool: str = DEFAULT_TOOL,
    cache_dir: Path | None = None,
) -> tuple[str, str]:
    """Return transcript name and normalized RNA sequence from one input source."""
    provided = [
        transcript_sequence is not None,
        transcript_file is not None,
        accession is not None,
    ]
    if sum(provided) != 1:
        raise ValueError("Provide exactly one of --target-sequence, --target-file, or --target-accession.")

    if accession:
        fasta_text = fetch_transcript_fasta(
            accession,
            email=email,
            tool=tool,
            cache_dir=cache_dir,
        )
        validate_single_transcript_record(fasta_text, f"NCBI record {accession}")
        return get_fasta_header(fasta_text) or accession, fasta_or_plain_text_to_sequence(fasta_text)

    if transcript_file:
        text = transcript_file.read_text(encoding="utf-8-sig")
        validate_single_transcript_record(text, str(transcript_file))
        return get_fasta_header(text) or transcript_file.name, fasta_or_plain_text_to_sequence(text)

    assert transcript_sequence is not None
    validate_single_transcript_record(transcript_sequence, "--target-sequence")
    return get_fasta_header(transcript_sequence) or "target_transcript", fasta_or_plain_text_to_sequence(transcript_sequence)


def prepare_pasted_transcript_sequence(text: str, target_name: str | None = None) -> str:
    """Return one canonical FASTA record for a pasted local transcript target."""
    validate_single_transcript_record(text, "Pasted transcript sequence")
    sequence = fasta_or_plain_text_to_sequence(text)
    header = clean_text_for_id(target_name or "") or get_fasta_header(text) or "pasted_transcript"
    return format_cached_transcript_fasta(header, sequence)


def local_transcript_target_from_args(args: argparse.Namespace) -> TranscriptTargetResult:
    """Build a ready target record from a pasted sequence or one local file."""
    transcript_name, sequence = read_transcript_input(
        transcript_sequence=args.target_sequence,
        transcript_file=args.target_file,
    )
    sequence_dna = normalize_dna(sequence)
    if args.target_file:
        source = "local file"
        source_path = str(args.target_file)
    else:
        source = "pasted sequence"
        source_path = ""
    return TranscriptTargetResult(
        requested_accession=transcript_name,
        transcript_name=transcript_name,
        sequence_5to3=sequence,
        sequence_length_nt=len(sequence),
        cache_path=source_path,
        cache_status=source,
        exact_version_match=False,
        sequence_sha256=hashlib.sha256(sequence_dna.encode("ascii")).hexdigest(),
        status="ready",
    )


def mismatch_positions(query: str, target: str) -> tuple[int, ...]:
    """Return 1-based mismatch positions between equal-length RNA strings."""
    return mismatch_positions_1based(query, target)


def scan_antisense_against_transcript(
    antisense_5to3: str,
    transcript_sequence: str,
    transcript_name: str = "target_transcript",
    antisense_name: str = "antisense_query",
    scan_region: AntisenseRegion | None = None,
    max_mismatches: int | None = DEFAULT_MAX_MISMATCHES,
    sequence_type: str = "AS",
) -> list[TranscriptMatch]:
    """Find AS or SS target windows in a transcript.

    The transcript window is reported in transcript 5'->3' orientation. AS
    oligos are reverse-complemented to the expected transcript target. SS
    oligos are compared directly to the transcript.
    """
    normalized_type = normalize_sequence_type(sequence_type)
    antisense = normalize_rna(antisense_5to3)
    transcript = normalize_rna(transcript_sequence)
    region = scan_region or AntisenseRegion("full")
    region_sequence, region_start, region_end = antisense_region_sequence(antisense, region)
    target = (
        get_complementary_sequence(region_sequence, reverse=True)
        if normalized_type == "AS"
        else region_sequence
    )
    if len(transcript) < len(target):
        return []

    matches = []
    for start_index in range(0, len(transcript) - len(target) + 1):
        window = transcript[start_index : start_index + len(target)]
        mismatches = mismatch_positions(target, window)
        if max_mismatches is None or len(mismatches) <= max_mismatches:
            transcript_match_as = (
                get_complementary_sequence(window, reverse=True)
                if normalized_type == "AS"
                else window
            )
            matches.append(
                TranscriptMatch(
                    transcript_name=transcript_name,
                    antisense_name=antisense_name,
                    scan_region=region.name,
                    as_region_start=region_start,
                    as_region_end=region_end,
                    antisense_5to3=antisense,
                    antisense_region_5to3=region_sequence,
                    target_5to3=target,
                    transcript_start=start_index + 1,
                    transcript_end=start_index + len(target),
                    mismatches=len(mismatches),
                    transcript_window_5to3=window,
                    transcript_match_as_5to3=transcript_match_as,
                    mismatch_positions_1based=mismatches,
                    as_mismatch_positions_1based=mismatch_positions(region_sequence, transcript_match_as),
                    sequence_type=normalized_type,
                )
            )
    return sorted(matches, key=lambda item: (item.mismatches, item.transcript_start))


def scan_sense_against_transcript(
    sense_5to3: str,
    transcript_sequence: str,
    transcript_name: str = "target_transcript",
    sense_name: str = "ss_query",
    scan_region: AntisenseRegion | None = None,
    max_mismatches: int | None = DEFAULT_MAX_MISMATCHES,
) -> list[TranscriptMatch]:
    """Find direct sense target windows in a transcript."""
    return scan_antisense_against_transcript(
        antisense_5to3=sense_5to3,
        transcript_sequence=transcript_sequence,
        transcript_name=transcript_name,
        antisense_name=sense_name,
        scan_region=scan_region,
        max_mismatches=max_mismatches,
        sequence_type="SS",
    )


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
                ";".join(str(position) for position in match.mismatch_positions_1based),
                ";".join(str(position) for position in match.as_mismatch_positions_1based),
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
    output.write("  ".join(header.ljust(widths[index]) for index, header in enumerate(headers)))
    output.write("\n")
    output.write("  ".join("-" * width for width in widths))
    output.write("\n")
    for row in text_rows:
        output.write("  ".join(value.ljust(widths[index]) for index, value in enumerate(row)))
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
            ";".join(str(position) for position in match.mismatch_positions_1based) or "-",
            ";".join(str(position) for position in match.as_mismatch_positions_1based) or "-",
        ]
        for match in match_list
    ]
    return terminal_table(headers, rows)


def closest_transcript_matches(matches: Iterable[TranscriptMatch], limit: int) -> list[TranscriptMatch]:
    if limit < 1:
        raise ValueError("--closest must be 1 or greater.")
    return sorted(
        matches,
        key=lambda match: (
            match.mismatches,
            match.transcript_start,
            match.antisense_name,
            match.scan_region,
        ),
    )[:limit]


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


def parse_blast_csv(text: str) -> list[dict[str, str]]:
    """Parse NCBI tabular CSV BLAST output into dictionaries."""
    rows = []
    for row in csv.reader(io.StringIO(text)):
        if not row or row[0].startswith("#"):
            continue
        if len(row) == len(CSV_COLUMNS):
            rows.append(dict(zip(CSV_COLUMNS, row)))
    return rows


def transcript_match_rows(
    matches: Iterable[TranscriptMatch],
    *,
    include_as_oriented_match: bool = True,
) -> list[dict[str, object]]:
    """Return local-match rows, optionally including the AS-oriented match.

    The default preserves the existing row schema for callers outside the Excel
    workbook.  Result workbooks omit the redundant AS-oriented sequence while
    CSV and terminal output continue to expose it.
    """
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


def comparison_result_for_region(
    *,
    input_order: int,
    query: AntisenseQuery,
    target_accession: str,
    scan_region: AntisenseRegion,
    passing_matches: Iterable[TranscriptMatch] = (),
    all_matches: Iterable[TranscriptMatch] = (),
    target_error: str = "",
) -> ComparisonResult:
    """Build one user-facing best comparison for a query/target/region."""
    query_region, region_start, region_end = antisense_region_sequence(
        query.sequence_5to3,
        scan_region,
    )
    passing = sorted(
        passing_matches,
        key=lambda match: (match.mismatches, match.transcript_start),
    )
    candidates = sorted(
        all_matches,
        key=lambda match: (match.mismatches, match.transcript_start),
    )
    best = passing[0] if passing else candidates[0] if candidates else None

    if target_error:
        result = "target_error"
    elif passing:
        result = "exact_match" if best and best.mismatches == 0 else "match"
    else:
        result = "no_match"

    mismatch_positions: tuple[int, ...] = ()
    differences = ""
    if best is not None:
        mismatch_positions = tuple(
            region_start + position - 1
            for position in best.as_mismatch_positions_1based
        )
        difference_items = []
        for relative_position, (expected_base, observed_base) in enumerate(
            zip(query_region, best.transcript_match_as_5to3),
            start=1,
        ):
            if expected_base != observed_base:
                full_position = region_start + relative_position - 1
                difference_items.append(
                    f"{full_position}:{expected_base}>{observed_base}"
                )
        differences = "; ".join(difference_items) or "None"

    return ComparisonResult(
        input_order=input_order,
        query_name=query.name,
        target_accession=target_accession,
        scan_region=scan_region.name,
        region_start=region_start,
        region_end=region_end,
        result=result,
        sites_within_threshold=len(passing),
        best_mismatches=best.mismatches if best is not None else None,
        mismatch_positions_in_query_1based=mismatch_positions,
        best_transcript_start=best.transcript_start if best is not None else None,
        best_transcript_end=best.transcript_end if best is not None else None,
        query_region_5to3=query_region,
        best_match_in_query_orientation_5to3=(
            best.transcript_match_as_5to3 if best is not None else ""
        ),
        differences=differences,
    )


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


def query_length_by_blast_id(queries: Iterable[AntisenseQuery]) -> dict[str, int]:
    return {
        query.blast_query_id: len(normalize_rna(query.sequence_5to3))
        for query in assign_unique_blast_query_ids(queries)
    }


def blast_raw_rows(
    batch_results: Iterable[BlastBatchResult],
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
                    "rid": result.submission.rid,
                    "batch_index": result.batch_index,
                    **row,
                    "query_length": query_length,
                    "alignment_fraction": (
                        alignment_length / query_length
                        if query_length
                        else None
                    ),
                }
            )
    return rows


def filter_blast_rows(
    rows: Iterable[dict[str, object]],
    max_mismatches: int,
    max_gap_opens: int,
    min_alignment_fraction: float,
) -> list[dict[str, object]]:
    filtered = []
    for row in rows:
        try:
            mismatches = int(float(row["mismatches"]))
            gap_opens = int(float(row["gap_opens"]))
            alignment_fraction = float(row["alignment_fraction"])
        except (TypeError, ValueError):
            continue
        if (
            mismatches <= max_mismatches
            and gap_opens <= max_gap_opens
            and alignment_fraction >= min_alignment_fraction
        ):
            filtered.append(row)
    return filtered


def blast_batch_rows(batch_results: Iterable[BlastBatchResult]) -> list[dict[str, object]]:
    rows = []
    for result in batch_results:
        sequences = [normalize_rna(query.sequence_5to3) for query in result.queries]
        rows.append(
            {
                "batch_index": result.batch_index,
                "rid": result.submission.rid,
                "rtoe_seconds": result.submission.rtoe_seconds,
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
        "sequence_types": ";".join(sorted({normalize_sequence_type(query.sequence_type) for query in queries})),
        "total_query_bases": sum(len(normalize_rna(query.sequence_5to3)) for query in queries),
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
        "panel_targets_ready": sum(target.status == "ready" for target in target_list),
        "panel_targets_error": sum(target.status == "error" for target in target_list),
    }
    return [{"key": key, "value": value} for key, value in metadata.items()]


def write_excel_workbook(path: Path, sheets: dict[str, list[dict[str, object]]]) -> None:
    import pandas as pd

    path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(path) as writer:
        for sheet_name, rows in sheets.items():
            pd.DataFrame(rows).to_excel(writer, sheet_name=sheet_name[:31], index=False)


def default_result_workbook(args: argparse.Namespace) -> Path | None:
    if args.result_workbook:
        return args.result_workbook
    source = args.as_table or args.as_file or getattr(args, "ss_table", None) or getattr(args, "ss_file", None)
    if source and not args.output and not args.blast_output:
        workflow = "ncbi_blast" if args.blast or args.blast_only else "ncbi_transcript_scan"
        return source.with_name(f"{source.stem}_{workflow}_results.xlsx")
    return None


def private_panel_requested(args: argparse.Namespace) -> bool:
    return bool(
        getattr(args, "private_panel", False)
        or getattr(args, "target_table", None)
        or getattr(args, "offline", False)
        or getattr(args, "refresh_targets", False)
        or getattr(args, "download_targets_only", False)
        or len(target_accession_values(args.target_accession)) > 1
    )


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


def private_panel_cache_dir(args: argparse.Namespace) -> Path:
    return args.cache_dir or Path(".ncbi_transcript_cache")


def application_base_dir() -> Path:
    """Return the portable app folder or the repository root during development."""
    if getattr(sys, "frozen", False):
        return Path(sys.executable).resolve().parent
    return Path(__file__).resolve().parents[2]


def application_data_dir() -> Path:
    """Return the writable data subfolder kept beside the portable app."""
    return application_base_dir() / APP_DATA_DIR_NAME


def gui_settings_path() -> Path:
    return application_data_dir() / GUI_SETTINGS_FILE_NAME


def gui_log_path() -> Path:
    return application_data_dir() / "logs" / GUI_LOG_FILE_NAME


def load_gui_settings() -> dict[str, object]:
    """Load portable per-user settings, returning an empty mapping if unavailable."""
    path = gui_settings_path()
    if not path.exists():
        return {}
    try:
        loaded = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return loaded if isinstance(loaded, dict) else {}


def save_gui_settings(settings: dict[str, object]) -> Path:
    """Atomically save portable settings beside the packaged application."""
    path = gui_settings_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = path.with_suffix(path.suffix + ".tmp")
    temporary_path.write_text(
        json.dumps(settings, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary_path.replace(path)
    return path


def shared_gui_transcript_cache_dir() -> Path:
    """Return the persistent transcript cache inside the portable data folder."""
    return application_data_dir() / "transcript_cache"


def default_gui_result_workbook(input_file: Path) -> Path:
    return input_file.with_name(f"{input_file.stem}_ncbi_transcript_scan_results.xlsx")


def write_result_workbook(
    path: Path,
    args: argparse.Namespace,
    queries: list[AntisenseQuery],
    scan_regions: list[AntisenseRegion],
    local_matches: list[TranscriptMatch],
    blast_results: list[BlastBatchResult],
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
    sheets = {
        "input_queries": input_query_rows(queries),
    }
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
        sheets["query_target_summary"] = query_target_summary_rows(query_target_summaries)
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


def excel_headers(input_file: Path, sheet_name: str | None = None) -> list[str]:
    import pandas as pd

    table = pd.read_excel(input_file, sheet_name=sheet_name or 0, nrows=0)
    return [str(column) for column in table.columns]


def choose_sheet_gui(root, input_file: Path) -> str | None:
    import tkinter as tk
    from tkinter import ttk

    sheets = list_excel_sheets(input_file)
    if len(sheets) <= 1:
        return None

    selected = {"value": sheets[0]}
    window = tk.Toplevel(root)
    window.title("Select transcript scan sheet")
    window.resizable(False, False)
    window.columnconfigure(1, weight=1)

    ttk.Label(window, text="Worksheet").grid(
        row=0, column=0, padx=16, pady=(16, 8), sticky="w"
    )
    sheet_var = tk.StringVar(value=sheets[0])
    sheet_box = ttk.Combobox(
        window,
        textvariable=sheet_var,
        values=sheets,
        state="readonly",
        width=max(30, min(60, max(len(sheet) for sheet in sheets) + 2)),
    )
    sheet_box.grid(row=0, column=1, padx=16, pady=(16, 8), sticky="ew")

    buttons = ttk.Frame(window)
    buttons.grid(row=1, column=0, columnspan=2, padx=16, pady=(8, 16), sticky="e")

    def use_sheet() -> None:
        selected["value"] = sheet_var.get()
        window.destroy()

    def cancel() -> None:
        selected["value"] = None
        window.destroy()

    ttk.Button(buttons, text="Cancel", command=cancel).grid(row=0, column=0, padx=(0, 8))
    ttk.Button(buttons, text="Continue", command=use_sheet).grid(row=0, column=1)
    window.protocol("WM_DELETE_WINDOW", cancel)
    window.bind("<Return>", lambda _event: use_sheet())
    window.bind("<Escape>", lambda _event: cancel())
    window.grab_set()
    sheet_box.focus_set()
    window.wait_window()
    return selected["value"]


def default_header(headers: list[str], candidates: list[str], fallback: str | None = None) -> str:
    by_lower = {header.lower(): header for header in headers}
    for candidate in candidates:
        if candidate.lower() in by_lower:
            return by_lower[candidate.lower()]
    return fallback if fallback is not None else headers[0]


def prompt_and_save_ncbi_email(root, current_email: str = "") -> str | None:
    """Ask for a valid NCBI contact email and save it in portable settings."""
    from tkinter import messagebox, simpledialog

    initial_value = clean_text_for_id(current_email)
    while True:
        entered = simpledialog.askstring(
            "NCBI contact email",
            "Enter your contact email for NCBI transcript requests.\n\n"
            "It is saved locally in TranscriptScanData\\settings.json.",
            initialvalue=initial_value,
            parent=root,
        )
        if entered is None:
            return None
        try:
            email = require_email(entered)
        except ValueError as error:
            messagebox.showerror("Invalid email", str(error), parent=root)
            initial_value = clean_text_for_id(entered)
            continue

        settings = load_gui_settings()
        settings["ncbi_email"] = email
        try:
            save_gui_settings(settings)
        except OSError as error:
            messagebox.showerror(
                "Cannot save settings",
                "The app folder must be writable so the email and transcript "
                f"cache can be saved.\n\n{error}",
                parent=root,
            )
            return None
        return email


def saved_or_prompted_ncbi_email(root) -> str | None:
    """Return the saved email, prompting on first use or invalid settings."""
    saved = clean_text_for_id(load_gui_settings().get("ncbi_email", ""))
    try:
        return require_email(saved)
    except ValueError:
        return prompt_and_save_ncbi_email(root, saved)


def choose_ncbi_gui_mode(root, ncbi_email: str) -> tuple[str | None, str]:
    """Choose the simple single-sequence workflow or the existing Excel workflow."""
    import tkinter as tk
    from tkinter import ttk

    selected = {"mode": None, "email": ncbi_email}
    window = tk.Toplevel(root)
    window.title("Private local transcript scan")
    window.resizable(False, False)

    ttk.Label(
        window,
        text="Choose an input workflow",
        font=("TkDefaultFont", 11, "bold"),
    ).grid(row=0, column=0, columnspan=2, padx=20, pady=(20, 12), sticky="w")
    ttk.Label(
        window,
        text=(
            "Both workflows compare locally. Transcript accessions may be sent "
            "to NCBI, but oligo sequences are not."
        ),
        wraplength=480,
    ).grid(row=1, column=0, columnspan=2, padx=20, pady=(0, 16), sticky="w")

    email_var = tk.StringVar(value=ncbi_email)
    email_frame = ttk.Frame(window)
    email_frame.grid(row=2, column=0, columnspan=2, padx=20, pady=(0, 16), sticky="ew")
    ttk.Label(email_frame, text="NCBI contact email:").grid(row=0, column=0, sticky="w")
    ttk.Label(email_frame, textvariable=email_var).grid(
        row=0, column=1, padx=(6, 12), sticky="w"
    )

    def change_email() -> None:
        changed = prompt_and_save_ncbi_email(window, str(selected["email"]))
        if changed:
            selected["email"] = changed
            email_var.set(changed)

    ttk.Button(email_frame, text="Change", command=change_email).grid(
        row=0, column=2, sticky="e"
    )

    def open_app_data_folder() -> None:
        try:
            data_dir = application_data_dir()
            data_dir.mkdir(parents=True, exist_ok=True)
            os.startfile(str(data_dir))
        except Exception as error:
            from tkinter import messagebox

            messagebox.showerror(
                "Cannot open app data folder",
                str(error),
                parent=window,
            )

    ttk.Button(
        email_frame,
        text="Open app data",
        command=open_app_data_folder,
    ).grid(row=0, column=3, padx=(8, 0), sticky="e")

    def choose(mode: str) -> None:
        selected["mode"] = mode
        window.destroy()

    ttk.Button(
        window,
        text="Single sequence and one transcript",
        command=lambda: choose("single"),
        width=34,
    ).grid(row=3, column=0, padx=(20, 8), pady=(0, 20), sticky="ew")
    ttk.Button(
        window,
        text="Excel sequence table",
        command=lambda: choose("excel"),
        width=26,
    ).grid(row=3, column=1, padx=(8, 20), pady=(0, 20), sticky="ew")

    window.protocol("WM_DELETE_WINDOW", window.destroy)
    window.bind("<Escape>", lambda _event: window.destroy())
    window.grab_set()
    window.wait_window()
    return selected["mode"], str(selected["email"])


def single_sequence_gui_draft(
    previous: dict[str, object] | None = None,
) -> dict[str, object]:
    """Return session-only single-scan values, retaining a previous form draft."""
    draft: dict[str, object] = {
        "sequence_type": "AS",
        "sequence_name": "",
        "sequence": "",
        "target_mode": "accession",
        "target_accession": "",
        "target_name": "",
        "target_sequence": "",
        "target_file": "",
        "scan_regions": ["full"],
        "max_mismatches": DEFAULT_MAX_MISMATCHES,
        "closest": DEFAULT_SINGLE_GUI_CLOSEST_MATCHES,
        "refresh_targets": False,
        "cache_dir": shared_gui_transcript_cache_dir(),
    }
    if previous:
        draft.update(previous)
        draft["scan_regions"] = list(previous.get("scan_regions", ["full"]))
    if draft.get("target_mode") not in {"accession", "paste", "file"}:
        draft["target_mode"] = "accession"
    return draft


def choose_single_sequence_gui_settings(
    root,
    previous: dict[str, object] | None = None,
) -> dict[str, object] | None:
    """Collect one local AS/SS-versus-transcript scan from the user."""
    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk

    selected: dict[str, object] = {}
    draft = single_sequence_gui_draft(previous)
    cache_dir = Path(draft.get("cache_dir") or shared_gui_transcript_cache_dir())
    scan_region_values = {str(value) for value in draft["scan_regions"]}

    window = tk.Toplevel(root)
    window.title("Single sequence transcript scan")
    window.minsize(760, 760)
    window.resizable(True, True)
    window.columnconfigure(1, weight=1)

    sequence_type_var = tk.StringVar(value=str(draft.get("sequence_type") or "AS"))
    sequence_name_var = tk.StringVar(value=str(draft.get("sequence_name") or ""))
    target_mode_var = tk.StringVar(value=str(draft.get("target_mode") or "accession"))
    accession_var = tk.StringVar(value=str(draft.get("target_accession") or ""))
    target_name_var = tk.StringVar(value=str(draft.get("target_name") or ""))
    target_file_var = tk.StringVar(value=str(draft.get("target_file") or ""))
    full_region_var = tk.BooleanVar(value="full" in scan_region_values)
    seed_region_var = tk.BooleanVar(value="seed:2-8" in scan_region_values)
    core_region_var = tk.BooleanVar(value="core:2-18" in scan_region_values)
    max_mismatches_var = tk.StringVar(value=str(draft["max_mismatches"]))
    closest_var = tk.StringVar(value=str(draft["closest"]))
    refresh_var = tk.BooleanVar(value=bool(draft.get("refresh_targets", False)))

    ttk.Label(window, text="Sequence type").grid(
        row=0, column=0, padx=16, pady=(16, 8), sticky="w"
    )
    ttk.Combobox(
        window,
        textvariable=sequence_type_var,
        values=["AS", "SS"],
        state="readonly",
        width=10,
    ).grid(row=0, column=1, padx=16, pady=(16, 8), sticky="w")

    ttk.Label(window, text="Sequence name (optional)").grid(
        row=1, column=0, padx=16, pady=8, sticky="w"
    )
    ttk.Entry(window, textvariable=sequence_name_var, width=52).grid(
        row=1, column=1, padx=16, pady=8, sticky="ew"
    )

    ttk.Label(window, text="Sequence, 5' to 3'").grid(
        row=2, column=0, padx=16, pady=8, sticky="nw"
    )
    sequence_text = tk.Text(window, width=56, height=4, wrap="word")
    sequence_text.grid(row=2, column=1, padx=16, pady=8, sticky="ew")
    sequence_text.insert("1.0", str(draft.get("sequence") or ""))

    target_frame = ttk.LabelFrame(window, text="Transcript target source")
    target_frame.grid(row=3, column=0, columnspan=2, padx=16, pady=8, sticky="nsew")
    target_frame.columnconfigure(1, weight=1)

    ttk.Radiobutton(
        target_frame,
        text="Download/reuse a RefSeq transcript accession",
        variable=target_mode_var,
        value="accession",
    ).grid(row=0, column=0, columnspan=3, padx=12, pady=(10, 4), sticky="w")
    ttk.Label(target_frame, text="Exact NM/XM/NR/XR version").grid(
        row=1, column=0, padx=(28, 8), pady=4, sticky="w"
    )
    accession_entry = ttk.Entry(target_frame, textvariable=accession_var, width=32)
    accession_entry.grid(row=1, column=1, padx=8, pady=4, sticky="ew")
    ttk.Label(target_frame, text="Example: NM_002439.5").grid(
        row=2, column=1, padx=8, pady=(0, 8), sticky="w"
    )

    ttk.Radiobutton(
        target_frame,
        text="Paste one transcript sequence (plain sequence or one FASTA record)",
        variable=target_mode_var,
        value="paste",
    ).grid(row=3, column=0, columnspan=3, padx=12, pady=(8, 4), sticky="w")
    ttk.Label(target_frame, text="Target name (optional)").grid(
        row=4, column=0, padx=(28, 8), pady=4, sticky="w"
    )
    target_name_entry = ttk.Entry(target_frame, textvariable=target_name_var, width=32)
    target_name_entry.grid(row=4, column=1, padx=8, pady=4, sticky="ew")
    target_sequence_text = tk.Text(target_frame, width=56, height=6, wrap="word")
    target_sequence_text.grid(
        row=5,
        column=0,
        columnspan=3,
        padx=(28, 12),
        pady=(4, 8),
        sticky="ew",
    )
    target_sequence_text.insert("1.0", str(draft.get("target_sequence") or ""))

    ttk.Radiobutton(
        target_frame,
        text="Load one transcript FASTA/text file",
        variable=target_mode_var,
        value="file",
    ).grid(row=6, column=0, columnspan=3, padx=12, pady=(8, 4), sticky="w")
    target_file_entry = ttk.Entry(target_frame, textvariable=target_file_var, width=42)
    target_file_entry.grid(row=7, column=0, columnspan=2, padx=(28, 8), pady=(4, 10), sticky="ew")

    def choose_target_file() -> None:
        path = filedialog.askopenfilename(
            parent=window,
            title="Select one transcript FASTA or text file",
            filetypes=[
                ("Sequence text files", "*.txt *.fa *.fasta *.fna *.ffn"),
                ("All files", "*.*"),
            ],
        )
        if path:
            target_file_var.set(path)
            target_mode_var.set("file")
            update_target_controls()

    browse_button = ttk.Button(target_frame, text="Browse", command=choose_target_file)
    browse_button.grid(row=7, column=2, padx=(0, 12), pady=(4, 10), sticky="e")

    ttk.Label(window, text="Scan regions").grid(
        row=4, column=0, padx=16, pady=8, sticky="nw"
    )
    region_frame = ttk.Frame(window)
    region_frame.grid(row=4, column=1, padx=16, pady=8, sticky="w")
    ttk.Checkbutton(region_frame, text="Full sequence", variable=full_region_var).grid(
        row=0, column=0, padx=(0, 16), sticky="w"
    )
    ttk.Checkbutton(
        region_frame,
        text="Seed (positions 2-8)",
        variable=seed_region_var,
    ).grid(row=0, column=1, padx=(0, 16), sticky="w")
    ttk.Checkbutton(
        region_frame,
        text="Core (positions 2-18)",
        variable=core_region_var,
    ).grid(row=0, column=2, sticky="w")

    ttk.Label(window, text="Maximum mismatches").grid(
        row=5, column=0, padx=16, pady=8, sticky="w"
    )
    ttk.Entry(window, textvariable=max_mismatches_var, width=8).grid(
        row=5, column=1, padx=16, pady=8, sticky="w"
    )

    ttk.Label(window, text="Closest windows per region").grid(
        row=6, column=0, padx=16, pady=8, sticky="w"
    )
    ttk.Entry(window, textvariable=closest_var, width=8).grid(
        row=6, column=1, padx=16, pady=8, sticky="w"
    )

    option_frame = ttk.Frame(window)
    option_frame.grid(row=7, column=0, columnspan=2, padx=16, pady=8, sticky="w")
    ttk.Label(
        option_frame,
        text=(
            "Accession mode uses the saved transcript when available and downloads "
            "it automatically when missing. Pasted and file targets stay local."
        ),
        wraplength=680,
    ).grid(row=0, column=0, sticky="w")
    refresh_check = ttk.Checkbutton(
        option_frame,
        text="Refresh this transcript from NCBI once",
        variable=refresh_var,
    )
    refresh_check.grid(row=1, column=0, pady=(4, 0), sticky="w")

    cache_frame = ttk.Frame(window)
    cache_frame.grid(row=8, column=0, columnspan=2, padx=16, pady=8, sticky="ew")
    ttk.Label(
        cache_frame,
        text=f"Shared transcript cache: {cache_dir}",
        wraplength=560,
    ).grid(row=0, column=0, padx=(0, 8), sticky="w")

    def open_cache_folder() -> None:
        try:
            cache_dir.mkdir(parents=True, exist_ok=True)
            os.startfile(str(cache_dir))
        except Exception as error:
            messagebox.showerror("Cannot open cache folder", str(error), parent=window)

    ttk.Button(cache_frame, text="Open cache folder", command=open_cache_folder).grid(
        row=0, column=1, sticky="e"
    )
    ttk.Label(
        window,
        text=(
            "Privacy: the oligo sequence and any pasted/local transcript stay on "
            "this computer. In accession mode, only accession/contact metadata are "
            "sent to NCBI. Form values are discarded when the app closes."
        ),
        wraplength=680,
    ).grid(row=9, column=0, columnspan=2, padx=16, pady=(4, 8), sticky="w")

    buttons = ttk.Frame(window)
    buttons.grid(row=10, column=0, columnspan=2, padx=16, pady=(8, 16), sticky="e")

    def update_target_controls() -> None:
        mode = target_mode_var.get()
        accession_entry.configure(state="normal" if mode == "accession" else "disabled")
        target_name_entry.configure(state="normal" if mode == "paste" else "disabled")
        target_sequence_text.configure(state="normal" if mode == "paste" else "disabled")
        target_file_entry.configure(state="normal" if mode == "file" else "disabled")
        browse_button.configure(state="normal" if mode == "file" else "disabled")
        refresh_check.configure(state="normal" if mode == "accession" else "disabled")

    target_mode_var.trace_add("write", lambda *_args: update_target_controls())
    update_target_controls()

    def use_settings() -> None:
        try:
            sequence = normalize_rna(sequence_text.get("1.0", "end").strip())
            mode = target_mode_var.get()
            if mode == "accession":
                accession = normalize_versioned_refseq_accession(accession_var.get())
            elif mode == "paste":
                prepare_pasted_transcript_sequence(
                    target_sequence_text.get("1.0", "end").strip(),
                    target_name_var.get(),
                )
                accession = accession_var.get().strip()
            elif mode == "file":
                target_file = Path(target_file_var.get().strip())
                if not target_file.exists() or not target_file.is_file():
                    raise ValueError(f"Transcript file does not exist: {target_file}")
                file_text = target_file.read_text(encoding="utf-8-sig")
                validate_single_transcript_record(file_text, str(target_file))
                fasta_or_plain_text_to_sequence(file_text)
                accession = accession_var.get().strip()
            else:
                raise ValueError("Choose a transcript target source.")

            scan_regions = []
            if full_region_var.get():
                scan_regions.append("full")
            if seed_region_var.get():
                scan_regions.append("seed:2-8")
            if core_region_var.get():
                scan_regions.append("core:2-18")
            if not scan_regions:
                raise ValueError("Select at least one scan region.")
            for region in parse_scan_regions(scan_regions):
                antisense_region_sequence(sequence, region)
            max_mismatches = int(max_mismatches_var.get())
            if max_mismatches < 0:
                raise ValueError("Maximum mismatches must be 0 or greater.")
            closest = int(closest_var.get())
            if closest < 1:
                raise ValueError("Closest windows must be 1 or greater.")
        except (OSError, UnicodeError, ValueError) as error:
            messagebox.showerror("Invalid settings", str(error), parent=window)
            return

        selected.update(
            {
                "sequence_type": sequence_type_var.get(),
                "sequence_name": sequence_name_var.get().strip(),
                "sequence": sequence,
                "target_mode": mode,
                "target_accession": accession,
                "target_name": target_name_var.get().strip(),
                "target_sequence": target_sequence_text.get("1.0", "end").strip(),
                "target_file": target_file_var.get().strip(),
                "scan_regions": scan_regions,
                "max_mismatches": max_mismatches,
                "closest": closest,
                "refresh_targets": refresh_var.get() if mode == "accession" else False,
                "cache_dir": cache_dir,
            }
        )
        window.destroy()

    ttk.Button(buttons, text="Cancel", command=window.destroy).grid(
        row=0, column=0, padx=(0, 8)
    )
    ttk.Button(buttons, text="Run local scan", command=use_settings).grid(
        row=0, column=1
    )
    window.protocol("WM_DELETE_WINDOW", window.destroy)
    window.bind("<Escape>", lambda _event: window.destroy())
    window.grab_set()
    window.wait_window()
    return selected or None


def single_sequence_gui_args(settings: dict[str, object]) -> argparse.Namespace:
    """Translate validated single-sequence GUI settings into normal CLI arguments."""
    sequence_type = normalize_sequence_type(str(settings["sequence_type"]))
    sequence_flag = "--as-sequence" if sequence_type == "AS" else "--ss-sequence"
    name_flag = "--as-name" if sequence_type == "AS" else "--ss-name"
    target_mode = str(settings.get("target_mode") or "accession")
    argv = [
        sequence_flag,
        str(settings["sequence"]),
        "--cache-dir",
        str(settings.get("cache_dir") or shared_gui_transcript_cache_dir()),
        "--email",
        str(settings.get("email", "")),
        "--max-mismatches",
        str(settings["max_mismatches"]),
        "--closest",
        str(settings.get("closest", DEFAULT_SINGLE_GUI_CLOSEST_MATCHES)),
    ]
    if target_mode == "accession":
        argv.extend(
            [
                "--private-panel",
                "--target-accession",
                normalize_versioned_refseq_accession(settings.get("target_accession", "")),
            ]
        )
    elif target_mode == "paste":
        argv.extend(
            [
                "--target-sequence",
                prepare_pasted_transcript_sequence(
                    str(settings.get("target_sequence") or ""),
                    str(settings.get("target_name") or ""),
                ),
            ]
        )
    elif target_mode == "file":
        target_file = Path(str(settings.get("target_file") or ""))
        if not str(settings.get("target_file") or "").strip():
            raise ValueError("Choose one transcript FASTA/text file.")
        argv.extend(["--target-file", str(target_file)])
    else:
        raise ValueError(f"Unsupported single-sequence target mode: {target_mode}")
    if settings.get("sequence_name"):
        argv.extend([name_flag, str(settings["sequence_name"])])
    for region in settings.get("scan_regions", ["full"]):
        argv.extend(["--scan-region", str(region)])
    if target_mode == "accession" and settings.get("refresh_targets"):
        argv.append("--refresh-targets")
    return build_parser().parse_args(argv)


def run_single_sequence_scan(
    args: argparse.Namespace,
    *,
    progress_callback: Callable[[int, int, str, str], None] | None = None,
    client: NcbiHttpClient | None = None,
) -> tuple[list[AntisenseQuery], list[AntisenseRegion], PrivatePanelScanResult]:
    """Run one private local query against one accession, pasted, or file target."""
    accessions = target_accession_values(args.target_accession)
    if accessions:
        if len(accessions) != 1:
            raise ValueError("Single-sequence mode requires exactly one transcript accession.")
        normalized_accessions = [normalize_versioned_refseq_accession(accessions[0])]
        targets = retrieve_transcript_targets(
            normalized_accessions,
            email=args.email,
            tool=args.tool,
            cache_dir=private_panel_cache_dir(args),
            offline=False,
            refresh=bool(getattr(args, "refresh_targets", False)),
            request_seconds=args.request_seconds,
            client=client,
            progress_callback=progress_callback,
        )
    else:
        targets = [local_transcript_target_from_args(args)]
        if progress_callback:
            progress_callback(
                1,
                1,
                targets[0].transcript_name,
                targets[0].cache_status,
            )
    target = targets[0]
    if target.status != "ready":
        target_label = accessions[0] if accessions else target.transcript_name
        raise ValueError(target.error or f"Transcript {target_label} could not be prepared.")

    queries = args_antisense_queries(args)
    if len(queries) != 1:
        raise ValueError("Single-sequence mode requires exactly one AS or SS sequence.")
    scan_regions = parse_scan_regions(args.scan_region)
    result = run_private_panel_scan(
        queries,
        targets,
        scan_regions,
        args.max_mismatches,
        closest=args.closest,
    )
    return queries, scan_regions, result


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
        + ("Used accession/cache workflow" if target.retrieved_accession else "Not used")
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
        "Closest windows are not filtered by the mismatch threshold and do not model gaps.\n"
    )
    for region in scan_regions:
        region_matches = [
            match for match in result.closest_matches if match.scan_region == region.name
        ]
        output.write(f"\nRegion: {region.name}\n")
        if region_matches:
            output.write(transcript_match_terminal_table(region_matches))
            output.write("\n")
        else:
            output.write("No windows available for this region.\n")
    return output.getvalue().rstrip() + "\n"


def show_single_sequence_result_gui(root, result_text: str) -> bool:
    """Show a scrollable result and return True when the user requests another scan."""
    import tkinter as tk
    from tkinter import ttk

    new_scan = {"value": False}
    window = tk.Toplevel(root)
    window.title("Single sequence transcript scan result")
    window.geometry("1100x700")
    window.minsize(760, 480)
    window.rowconfigure(0, weight=1)
    window.columnconfigure(0, weight=1)

    result_frame = ttk.Frame(window)
    result_frame.grid(row=0, column=0, padx=12, pady=12, sticky="nsew")
    result_frame.rowconfigure(0, weight=1)
    result_frame.columnconfigure(0, weight=1)
    output_text = tk.Text(
        result_frame,
        wrap="none",
        font="TkFixedFont",
        padx=8,
        pady=8,
    )
    vertical = ttk.Scrollbar(result_frame, orient="vertical", command=output_text.yview)
    horizontal = ttk.Scrollbar(result_frame, orient="horizontal", command=output_text.xview)
    output_text.configure(yscrollcommand=vertical.set, xscrollcommand=horizontal.set)
    output_text.grid(row=0, column=0, sticky="nsew")
    vertical.grid(row=0, column=1, sticky="ns")
    horizontal.grid(row=1, column=0, sticky="ew")
    output_text.insert("1.0", result_text)
    output_text.configure(state="disabled")

    buttons = ttk.Frame(window)
    buttons.grid(row=1, column=0, padx=12, pady=(0, 12), sticky="e")

    def copy_all() -> None:
        root.clipboard_clear()
        root.clipboard_append(result_text)
        root.update()

    def request_new_scan() -> None:
        new_scan["value"] = True
        window.destroy()

    ttk.Button(buttons, text="Copy all", command=copy_all).grid(
        row=0, column=0, padx=(0, 8)
    )
    ttk.Button(buttons, text="Edit and run again", command=request_new_scan).grid(
        row=0, column=1, padx=(0, 8)
    )
    ttk.Button(buttons, text="Close", command=window.destroy).grid(row=0, column=2)

    window.protocol("WM_DELETE_WINDOW", window.destroy)
    window.bind("<Escape>", lambda _event: window.destroy())
    window.grab_set()
    window.wait_window()
    return new_scan["value"]


def run_single_sequence_gui(root, ncbi_email: str) -> int:
    """Run the one-sequence/one-transcript private local GUI workflow."""
    import tkinter as tk
    from tkinter import messagebox, ttk

    draft = single_sequence_gui_draft()
    while True:
        settings = choose_single_sequence_gui_settings(root, draft)
        if not settings:
            return 0
        settings["email"] = ncbi_email
        draft = {key: value for key, value in settings.items() if key != "email"}
        try:
            args = single_sequence_gui_args(settings)
            validate_runtime_args(args)
        except Exception as error:
            logging.exception("Invalid single-sequence scan settings")
            messagebox.showerror(
                "Single sequence transcript scan failed",
                str(error),
                parent=root,
            )
            continue

        progress_window = tk.Toplevel(root)
        progress_window.title("Single sequence transcript scan")
        progress_window.resizable(False, False)
        status_var = tk.StringVar(value="Preparing transcript reference...")
        ttk.Label(progress_window, textvariable=status_var, width=64).grid(
            row=0, column=0, padx=16, pady=(16, 8), sticky="w"
        )
        progress_bar = ttk.Progressbar(
            progress_window,
            mode="determinate",
            length=440,
            maximum=1,
        )
        progress_bar.grid(row=1, column=0, padx=16, pady=(8, 16), sticky="ew")

        def update_progress(
            completed: int,
            total: int,
            accession: str,
            status: str,
        ) -> None:
            progress_bar.configure(maximum=max(total, 1), value=completed)
            status_var.set(f"{accession}: {status}")
            progress_window.update()

        try:
            try:
                queries, scan_regions, result = run_single_sequence_scan(
                    args,
                    progress_callback=update_progress,
                )
            finally:
                progress_window.destroy()
        except Exception as error:
            logging.exception("Single-sequence transcript scan failed")
            messagebox.showerror(
                "Single sequence transcript scan failed",
                f"{error}\n\nYour form entries have been retained.",
                parent=root,
            )
            continue

        result_text = format_single_sequence_scan_result(
            args,
            queries,
            scan_regions,
            result,
        )
        draft["refresh_targets"] = False
        if not show_single_sequence_result_gui(root, result_text):
            return 0


def choose_ncbi_gui_settings(root, headers: list[str]) -> dict[str, object] | None:
    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk

    selected: dict[str, object] = {}
    target_accession_default = default_header(
        headers,
        ["target_accession", "target accession", "refseq", "refseq accession"],
        "",
    )

    window = tk.Toplevel(root)
    window.title("Private local transcript scan settings")
    window.resizable(False, False)
    window.columnconfigure(1, weight=1)

    sequence_type_var = tk.StringVar(value="AS")
    as_column_var = tk.StringVar(
        value=default_header(headers, ["AS_5to3", "antisense", "as", "sequence"])
    )
    name_column_var = tk.StringVar(
        value=default_header(headers, ["oligo_id", "oligo id", "id", "name"], "")
    )
    target_mode_var = tk.StringVar(
        value="column" if target_accession_default else "refseq"
    )
    target_column_var = tk.StringVar(value=target_accession_default or headers[0])
    refseq_var = tk.StringVar()
    transcript_file_var = tk.StringVar()
    panel_accessions_var = tk.StringVar()
    target_table_var = tk.StringVar()
    scan_regions_var = tk.StringVar(value="full")
    max_mismatches_var = tk.StringVar(value=str(DEFAULT_MAX_MISMATCHES))
    offline_var = tk.BooleanVar(value=False)

    ttk.Label(window, text="Sequence type").grid(
        row=0, column=0, padx=16, pady=(16, 8), sticky="w"
    )
    ttk.Combobox(
        window,
        textvariable=sequence_type_var,
        values=["AS", "SS"],
        state="readonly",
        width=12,
    ).grid(row=0, column=1, padx=16, pady=(16, 8), sticky="w")

    ttk.Label(window, text="Sequence column").grid(
        row=1, column=0, padx=16, pady=8, sticky="w"
    )
    ttk.Combobox(
        window,
        textvariable=as_column_var,
        values=headers,
        state="readonly",
        width=max(30, min(60, max(len(header) for header in headers) + 2)),
    ).grid(row=1, column=1, padx=16, pady=8, sticky="ew")

    ttk.Label(window, text="Sequence name column").grid(
        row=2, column=0, padx=16, pady=8, sticky="w"
    )
    ttk.Combobox(
        window,
        textvariable=name_column_var,
        values=["", *headers],
        state="readonly",
        width=36,
    ).grid(row=2, column=1, padx=16, pady=8, sticky="ew")

    ttk.Label(window, text="Transcript source").grid(
        row=3, column=0, padx=16, pady=8, sticky="nw"
    )
    source_frame = ttk.Frame(window)
    source_frame.grid(row=3, column=1, padx=16, pady=8, sticky="ew")
    source_frame.columnconfigure(1, weight=1)

    ttk.Radiobutton(
        source_frame,
        text="Use target_accession column",
        variable=target_mode_var,
        value="column",
    ).grid(row=0, column=0, columnspan=2, sticky="w")
    ttk.Combobox(
        source_frame,
        textvariable=target_column_var,
        values=headers,
        state="readonly",
        width=36,
    ).grid(row=1, column=0, columnspan=2, pady=(4, 8), sticky="ew")

    ttk.Radiobutton(
        source_frame,
        text="Use this RefSeq accession for all rows",
        variable=target_mode_var,
        value="refseq",
    ).grid(row=2, column=0, columnspan=2, sticky="w")
    ttk.Entry(source_frame, textvariable=refseq_var, width=38).grid(
        row=3, column=0, columnspan=2, pady=(4, 8), sticky="ew"
    )

    ttk.Radiobutton(
        source_frame,
        text="Use transcript FASTA/text file for all rows",
        variable=target_mode_var,
        value="file",
    ).grid(row=4, column=0, columnspan=2, sticky="w")

    def choose_transcript_file() -> None:
        path = filedialog.askopenfilename(
            parent=window,
            title="Select transcript FASTA or text file",
            filetypes=[
                ("Sequence text files", "*.txt *.fa *.fasta *.fna *.ffn"),
                ("FASTA files", "*.fa *.fasta *.fna *.ffn"),
                ("Text files", "*.txt"),
                ("All files", "*.*"),
            ],
        )
        if path:
            transcript_file_var.set(path)
            target_mode_var.set("file")

    ttk.Entry(source_frame, textvariable=transcript_file_var, width=38).grid(
        row=5, column=0, pady=(4, 0), sticky="ew"
    )
    ttk.Button(source_frame, text="Browse", command=choose_transcript_file).grid(
        row=5, column=1, padx=(8, 0), pady=(4, 0), sticky="e"
    )

    ttk.Radiobutton(
        source_frame,
        text="Private panel: versioned RefSeq accessions",
        variable=target_mode_var,
        value="panel",
    ).grid(row=6, column=0, columnspan=2, pady=(8, 0), sticky="w")
    ttk.Entry(source_frame, textvariable=panel_accessions_var, width=38).grid(
        row=7, column=0, columnspan=2, pady=(4, 8), sticky="ew"
    )

    ttk.Radiobutton(
        source_frame,
        text="Private panel: accession table",
        variable=target_mode_var,
        value="panel_table",
    ).grid(row=8, column=0, columnspan=2, sticky="w")

    def choose_target_table() -> None:
        path = filedialog.askopenfilename(
            parent=window,
            title="Select target accession table",
            filetypes=[
                ("Target lists", "*.txt *.list *.csv *.tsv *.xlsx *.xlsm *.xls"),
                ("All files", "*.*"),
            ],
        )
        if path:
            target_table_var.set(path)
            target_mode_var.set("panel_table")

    ttk.Entry(source_frame, textvariable=target_table_var, width=38).grid(
        row=9, column=0, pady=(4, 0), sticky="ew"
    )
    ttk.Button(source_frame, text="Browse", command=choose_target_table).grid(
        row=9, column=1, padx=(8, 0), pady=(4, 0), sticky="e"
    )

    ttk.Label(window, text="Scan regions").grid(
        row=4, column=0, padx=16, pady=8, sticky="w"
    )
    ttk.Entry(window, textvariable=scan_regions_var, width=38).grid(
        row=4, column=1, padx=16, pady=8, sticky="ew"
    )

    ttk.Label(window, text="Max mismatches").grid(
        row=5, column=0, padx=16, pady=8, sticky="w"
    )
    ttk.Entry(window, textvariable=max_mismatches_var, width=12).grid(
        row=5, column=1, padx=16, pady=8, sticky="w"
    )

    ttk.Checkbutton(
        window,
        text="Offline: require every panel transcript in the local cache",
        variable=offline_var,
    ).grid(row=6, column=0, columnspan=2, padx=16, pady=8, sticky="w")
    ttk.Label(
        window,
        text="Private panel mode sends transcript accessions to NCBI; guide sequences remain local.",
    ).grid(row=7, column=0, columnspan=2, padx=16, pady=(0, 8), sticky="w")

    buttons = ttk.Frame(window)
    buttons.grid(row=8, column=0, columnspan=2, padx=16, pady=(8, 16), sticky="e")

    def use_settings() -> None:
        try:
            max_mismatches = int(max_mismatches_var.get())
            if max_mismatches < 0:
                raise ValueError
            scan_regions = [
                value.strip()
                for value in re.split(r"[;,]", scan_regions_var.get())
                if value.strip()
            ]
            parse_scan_regions(scan_regions)
        except ValueError as error:
            messagebox.showerror("Invalid settings", str(error), parent=window)
            return

        mode = target_mode_var.get()
        if mode == "column" and not target_column_var.get():
            messagebox.showerror("Missing target accession column", "Choose a target accession column.", parent=window)
            return
        if mode == "refseq" and not refseq_var.get().strip():
            messagebox.showerror("Missing RefSeq", "Enter a RefSeq accession.", parent=window)
            return
        if mode == "file" and not transcript_file_var.get().strip():
            messagebox.showerror("Missing transcript file", "Choose a transcript FASTA/text file.", parent=window)
            return
        panel_accessions = [
            value
            for value in re.split(r"[,;\s]+", panel_accessions_var.get().strip())
            if value
        ]
        if mode == "panel":
            if not panel_accessions:
                messagebox.showerror(
                    "Missing panel accessions",
                    "Enter one or more versioned RefSeq accessions.",
                    parent=window,
                )
                return
            try:
                panel_accessions = [
                    normalize_versioned_refseq_accession(value)
                    for value in panel_accessions
                ]
            except ValueError as error:
                messagebox.showerror("Invalid panel accession", str(error), parent=window)
                return
        if mode == "panel_table" and not target_table_var.get().strip():
            messagebox.showerror(
                "Missing target table",
                "Choose a transcript accession table.",
                parent=window,
            )
            return
        if offline_var.get() and mode not in {"panel", "panel_table"}:
            messagebox.showerror(
                "Offline panel mode",
                "Offline mode is available only for private transcript panels.",
                parent=window,
            )
            return

        selected.update(
            {
                "sequence_type": sequence_type_var.get(),
                "as_column": as_column_var.get(),
                "as_name_column": name_column_var.get() or None,
                "target_mode": mode,
                "target_accession_column": target_column_var.get() if mode == "column" else None,
                "target_accession": (
                    refseq_var.get().strip()
                    if mode == "refseq"
                    else panel_accessions if mode == "panel" else None
                ),
                "target_file": Path(transcript_file_var.get()) if mode == "file" else None,
                "target_table": (
                    Path(target_table_var.get()) if mode == "panel_table" else None
                ),
                "private_panel": mode in {"panel", "panel_table"},
                "offline": offline_var.get(),
                "scan_regions": scan_regions,
                "max_mismatches": max_mismatches,
            }
        )
        window.destroy()

    def cancel() -> None:
        selected.clear()
        window.destroy()

    ttk.Button(buttons, text="Cancel", command=cancel).grid(row=0, column=0, padx=(0, 8))
    ttk.Button(buttons, text="Run", command=use_settings).grid(row=0, column=1)

    window.protocol("WM_DELETE_WINDOW", cancel)
    window.bind("<Return>", lambda _event: use_settings())
    window.bind("<Escape>", lambda _event: cancel())
    window.grab_set()
    window.wait_window()
    return selected or None


def gui_args(input_file: Path, sheet_name: str | None, settings: dict[str, object]) -> argparse.Namespace:
    sequence_type = normalize_sequence_type(str(settings.get("sequence_type", "AS")))
    private_panel = bool(settings.get("private_panel", False))
    output_path = (
        input_file.with_name(f"{input_file.stem}_private_transcript_panel_results.xlsx")
        if private_panel
        else default_gui_result_workbook(input_file)
    )
    return argparse.Namespace(
        as_sequence=None,
        as_name=None,
        as_file=None,
        as_table=input_file if sequence_type == "AS" else None,
        as_column=settings["as_column"] if sequence_type == "AS" else None,
        as_name_column=settings["as_name_column"] if sequence_type == "AS" else None,
        as_sheet=sheet_name if sequence_type == "AS" else None,
        ss_sequence=None,
        ss_name=None,
        ss_file=None,
        ss_table=input_file if sequence_type == "SS" else None,
        ss_column=settings["as_column"] if sequence_type == "SS" else None,
        ss_name_column=settings["as_name_column"] if sequence_type == "SS" else None,
        ss_sheet=sheet_name if sequence_type == "SS" else None,
        target_accession=settings["target_accession"],
        target_accession_column=settings["target_accession_column"],
        target_file=settings["target_file"],
        target_sequence=None,
        target_table=settings.get("target_table"),
        target_column=None,
        target_sheet=None,
        private_panel=private_panel,
        offline=bool(settings.get("offline", False)),
        refresh_targets=False,
        download_targets_only=False,
        scan_region=settings["scan_regions"],
        max_mismatches=settings["max_mismatches"],
        email=str(settings.get("email", "")),
        tool=DEFAULT_TOOL,
        blast=False,
        blast_only=False,
        database=DEFAULT_DATABASE,
        expect=DEFAULT_EXPECT,
        word_size=DEFAULT_WORD_SIZE,
        hitlist_size=DEFAULT_HITLIST_SIZE,
        megablast=False,
        timeout_seconds=1800,
        max_batch_bases=DEFAULT_BATCH_BASES,
        request_seconds=DEFAULT_REQUEST_SECONDS,
        poll_seconds=DEFAULT_POLL_SECONDS,
        filter_max_mismatches=DEFAULT_MAX_MISMATCHES,
        filter_max_gap_opens=0,
        filter_min_alignment_fraction=0.8,
        cache_dir=shared_gui_transcript_cache_dir(),
        rid_log=None,
        output=None,
        blast_output=None,
        terminal=False,
        stdout_csv=False,
        closest=None,
        result_workbook=output_path,
        gui=False,
    )


def run_gui() -> int:
    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk

    root = tk.Tk()
    root.withdraw()

    try:
        ncbi_email = saved_or_prompted_ncbi_email(root)
        if ncbi_email is None:
            return 0
        gui_mode, ncbi_email = choose_ncbi_gui_mode(root, ncbi_email)
        if gui_mode is None:
            return 0
        if gui_mode == "single":
            return run_single_sequence_gui(root, ncbi_email)

        input_path = filedialog.askopenfilename(
            title="Select local transcript scan Excel input file",
            filetypes=[
                ("Excel files", "*.xlsx *.xlsm *.xls"),
                ("All files", "*.*"),
            ],
        )
        if not input_path:
            return 0

        input_file = Path(input_path)
        sheet_name = choose_sheet_gui(root, input_file)
        if sheet_name is None and len(list_excel_sheets(input_file)) > 1:
            return 0

        headers = excel_headers(input_file, sheet_name)
        if not headers:
            messagebox.showerror("No columns", "The selected sheet has no headers.")
            return 1

        settings = choose_ncbi_gui_settings(root, headers)
        if not settings:
            return 0
        settings["email"] = ncbi_email

        args = gui_args(input_file, sheet_name, settings)
        validate_runtime_args(args)
        started_at = datetime.now(timezone.utc).isoformat()
        if private_panel_requested(args):
            progress_window = tk.Toplevel(root)
            progress_window.title("Private transcript panel")
            progress_window.resizable(False, False)
            status_var = tk.StringVar(value="Preparing transcript references...")
            ttk.Label(progress_window, textvariable=status_var, width=64).grid(
                row=0,
                column=0,
                padx=16,
                pady=(16, 8),
                sticky="w",
            )
            progress_bar = ttk.Progressbar(
                progress_window,
                mode="determinate",
                length=440,
            )
            progress_bar.grid(row=1, column=0, padx=16, pady=8, sticky="ew")
            cancelled = {"value": False}

            def request_cancel() -> None:
                cancelled["value"] = True
                status_var.set("Cancelling after the current retrieval request...")

            ttk.Button(
                progress_window,
                text="Cancel",
                command=request_cancel,
            ).grid(row=2, column=0, padx=16, pady=(8, 16), sticky="e")
            progress_window.protocol("WM_DELETE_WINDOW", request_cancel)

            def update_progress(
                completed: int,
                total: int,
                accession: str,
                status: str,
            ) -> None:
                progress_bar.configure(maximum=max(total, 1), value=completed)
                status_var.set(
                    f"{completed}/{total} targets - {accession}: {status}"
                )
                progress_window.update()

            try:
                exit_code = run_private_panel_workflow(
                    args,
                    started_at,
                    progress_callback=update_progress,
                    cancel_check=lambda: cancelled["value"],
                    include_comparison_results=True,
                )
            finally:
                progress_window.destroy()
            if cancelled["value"]:
                messagebox.showwarning(
                    "Private transcript panel cancelled",
                    "The partial status workbook was saved. Guide sequences "
                    f"remained local.\n\n{args.result_workbook}",
                )
            elif exit_code == 0:
                messagebox.showinfo(
                    "Private transcript panel complete",
                    "Guide sequences remained local.\n\n"
                    f"Wrote panel workbook to:\n{args.result_workbook}",
                )
            elif exit_code == 2:
                messagebox.showwarning(
                    "Private transcript panel completed with target errors",
                    "Some targets could not be scanned. Review comparison_results "
                    f"and transcript_targets in:\n{args.result_workbook}",
                )
            else:
                messagebox.showwarning(
                    "Private transcript panel incomplete",
                    "No target transcript was ready. Review the transcript_targets "
                    f"sheet in:\n{args.result_workbook}",
                )
            return exit_code
        queries = args_antisense_queries(args)
        scan_regions = parse_scan_regions(args.scan_region)
        local_matches, comparison_results = run_local_scan_with_comparison(
            args,
            queries,
            scan_regions,
        )
        completed_at = datetime.now(timezone.utc).isoformat()
        write_result_workbook(
            args.result_workbook,
            args,
            queries,
            scan_regions,
            local_matches,
            [],
            started_at,
            completed_at,
            include_blast_sheets=False,
            comparison_results=comparison_results,
        )
        messagebox.showinfo(
            "Local transcript scan complete",
            f"Wrote local transcript scan workbook to:\n{args.result_workbook}",
        )
        return 0
    except Exception as error:
        logging.exception("Local transcript scan failed")
        messagebox.showerror(
            "Local transcript scan failed",
            f"{error}\n\nDiagnostic log:\n{gui_log_path()}",
        )
        return 1
    finally:
        root.destroy()


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
    parser.add_argument(
        "--as-name-column",
        help="Optional AS name/id column for --as-table.",
    )
    parser.add_argument("--as-sheet", help="Excel worksheet for --as-table. Defaults to first sheet.")
    parser.add_argument("--ss-sequence", help="One SS/sense oligo sequence in 5'->3' transcript orientation.")
    parser.add_argument("--ss-name", help="Name for the single --ss-sequence input.")
    parser.add_argument("--ss-file", type=Path, help="Text or FASTA file containing SS/sense sequences.")
    parser.add_argument("--ss-table", type=Path, help="Excel/CSV table containing SS/sense sequences.")
    parser.add_argument(
        "--ss-column",
        help="SS/sense sequence column for --ss-table. Defaults to sense/ss/sequence or first column.",
    )
    parser.add_argument(
        "--ss-name-column",
        help="Optional SS/sense name/id column for --ss-table.",
    )
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
            "Column in --as-table or --ss-table containing per-row NM/XM/NR/XR accessions. "
            "Defaults to target_accession if present."
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


def run_local_scan(
    args: argparse.Namespace,
    queries: list[AntisenseQuery],
    scan_regions: list[AntisenseRegion],
    max_mismatches: int | None = DEFAULT_MAX_MISMATCHES,
) -> list[TranscriptMatch]:
    if max_mismatches == DEFAULT_MAX_MISMATCHES:
        max_mismatches = args.max_mismatches
    matches = []
    shared_transcript: tuple[str, str] | None = None
    if args.target_accession_column:
        missing_accessions = [query.name for query in queries if not query.target_accession]
        if missing_accessions:
            names = ", ".join(missing_accessions)
            raise ValueError(
                "Missing target accession for query row(s): "
                f"{names}. Fill --target-accession-column for every query."
            )
    if not args.target_accession_column:
        accessions = target_accession_values(args.target_accession)
        if len(accessions) > 1:
            raise ValueError(
                "Multiple --target-accession values require private panel mode."
            )
        shared_transcript = read_transcript_input(
            transcript_sequence=args.target_sequence,
            transcript_file=args.target_file,
            accession=accessions[0] if accessions else None,
            email=args.email,
            tool=args.tool,
            cache_dir=args.cache_dir,
        )

    for query in queries:
        if args.target_accession_column:
            transcript_name, transcript = read_transcript_input(
                accession=query.target_accession,
                email=args.email,
                tool=args.tool,
                cache_dir=args.cache_dir,
            )
        else:
            assert shared_transcript is not None
            transcript_name, transcript = shared_transcript

        for scan_region in scan_regions:
            matches.extend(
                scan_antisense_against_transcript(
                    query.sequence_5to3,
                    transcript,
                    transcript_name=transcript_name,
                    antisense_name=query.name,
                    scan_region=scan_region,
                    max_mismatches=max_mismatches,
                    sequence_type=query.sequence_type,
                )
            )
    return matches


def run_local_scan_with_comparison(
    args: argparse.Namespace,
    queries: list[AntisenseQuery],
    scan_regions: list[AntisenseRegion],
) -> tuple[list[TranscriptMatch], list[ComparisonResult]]:
    """Run the GUI table scan and build one compact result per selected region."""
    matches: list[TranscriptMatch] = []
    comparison_results: list[ComparisonResult] = []
    shared_transcript: tuple[str, str] | None = None
    shared_target_label = ""

    if args.target_accession_column:
        missing_accessions = [query.name for query in queries if not query.target_accession]
        if missing_accessions:
            names = ", ".join(missing_accessions)
            raise ValueError(
                "Missing target accession for query row(s): "
                f"{names}. Fill --target-accession-column for every query."
            )
    else:
        accessions = target_accession_values(args.target_accession)
        if len(accessions) > 1:
            raise ValueError(
                "Multiple --target-accession values require private panel mode."
            )
        shared_transcript = read_transcript_input(
            transcript_sequence=args.target_sequence,
            transcript_file=args.target_file,
            accession=accessions[0] if accessions else None,
            email=args.email,
            tool=args.tool,
            cache_dir=args.cache_dir,
        )
        shared_target_label = accessions[0] if accessions else shared_transcript[0]

    for input_order, query in enumerate(queries, start=1):
        if args.target_accession_column:
            transcript_name, transcript = read_transcript_input(
                accession=query.target_accession,
                email=args.email,
                tool=args.tool,
                cache_dir=args.cache_dir,
            )
            target_label = query.target_accession
        else:
            assert shared_transcript is not None
            transcript_name, transcript = shared_transcript
            target_label = shared_target_label

        for scan_region in scan_regions:
            all_region_matches = scan_antisense_against_transcript(
                query.sequence_5to3,
                transcript,
                transcript_name=transcript_name,
                antisense_name=query.name,
                scan_region=scan_region,
                max_mismatches=None,
                sequence_type=query.sequence_type,
            )
            passing_matches = [
                match
                for match in all_region_matches
                if match.mismatches <= args.max_mismatches
            ]
            matches.extend(passing_matches)
            comparison_results.append(
                comparison_result_for_region(
                    input_order=input_order,
                    query=query,
                    target_accession=target_label,
                    scan_region=scan_region,
                    passing_matches=passing_matches,
                    all_matches=all_region_matches,
                )
            )

    return matches, comparison_results


def run_private_panel_scan(
    queries: list[AntisenseQuery],
    targets: list[TranscriptTargetResult],
    scan_regions: list[AntisenseRegion],
    max_mismatches: int,
    closest: int | None = None,
) -> PrivatePanelScanResult:
    """Scan every private guide against every transcript target locally."""
    matches: list[TranscriptMatch] = []
    panel_closest_matches: list[TranscriptMatch] = []
    summaries: list[QueryTargetSummary] = []
    comparison_results: list[ComparisonResult] = []
    region_names = ";".join(region.name for region in scan_regions)

    for input_order, query in enumerate(queries, start=1):
        sequence_type = normalize_sequence_type(query.sequence_type)
        for target in targets:
            if target.status != "ready":
                for scan_region in scan_regions:
                    comparison_results.append(
                        comparison_result_for_region(
                            input_order=input_order,
                            query=query,
                            target_accession=target.requested_accession,
                            scan_region=scan_region,
                            target_error=target.error,
                        )
                    )
                summaries.append(
                    QueryTargetSummary(
                        query_name=query.name,
                        sequence_type=sequence_type,
                        requested_accession=target.requested_accession,
                        retrieved_accession=target.retrieved_accession,
                        target_status=target.status,
                        scan_status="target_error",
                        scan_regions=region_names,
                        match_count=0,
                        exact_match_count=0,
                        best_mismatches=None,
                        error=target.error,
                    )
                )
                continue

            pair_matches = []
            pair_all_matches = []
            for scan_region in scan_regions:
                region_matches = scan_antisense_against_transcript(
                    query.sequence_5to3,
                    target.sequence_5to3,
                    transcript_name=target.transcript_name or target.retrieved_accession,
                    antisense_name=query.name,
                    scan_region=scan_region,
                    max_mismatches=None,
                    sequence_type=sequence_type,
                )
                pair_all_matches.extend(region_matches)
                passing_region_matches = [
                    match
                    for match in region_matches
                    if match.mismatches <= max_mismatches
                ]
                pair_matches.extend(passing_region_matches)
                comparison_results.append(
                    comparison_result_for_region(
                        input_order=input_order,
                        query=query,
                        target_accession=target.requested_accession,
                        scan_region=scan_region,
                        passing_matches=passing_region_matches,
                        all_matches=region_matches,
                    )
                )
                if closest is not None:
                    panel_closest_matches.extend(
                        closest_transcript_matches(region_matches, closest)
                    )
            matches.extend(pair_matches)
            summaries.append(
                QueryTargetSummary(
                    query_name=query.name,
                    sequence_type=sequence_type,
                    requested_accession=target.requested_accession,
                    retrieved_accession=target.retrieved_accession,
                    target_status=target.status,
                    scan_status="matched" if pair_matches else "no_match",
                    scan_regions=region_names,
                    match_count=len(pair_matches),
                    exact_match_count=sum(match.mismatches == 0 for match in pair_matches),
                    best_mismatches=(
                        min(
                            match.mismatches
                            for match in (
                                pair_all_matches if closest is not None else pair_matches
                            )
                        )
                        if (pair_all_matches if closest is not None else pair_matches)
                        else None
                    ),
                )
            )

    return PrivatePanelScanResult(
        targets=tuple(targets),
        matches=tuple(matches),
        summaries=tuple(summaries),
        closest_matches=tuple(panel_closest_matches),
        comparison_results=tuple(comparison_results),
    )


def combine_blast_csv(outputs: Iterable[BlastBatchResult]) -> str:
    combined = io.StringIO()
    writer = csv.writer(combined, lineterminator="\n")
    writer.writerow(["rid", *CSV_COLUMNS])
    for result in outputs:
        for row in parse_blast_csv(result.csv_text):
            writer.writerow([result.submission.rid, *[row[column] for column in CSV_COLUMNS]])
    return combined.getvalue()


def append_rid_log(path: Path, result: BlastBatchResult) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    write_header = not path.exists()
    with path.open("a", encoding="utf-8", newline="") as output_file:
        writer = csv.writer(output_file)
        if write_header:
            writer.writerow(
                [
                    "batch_index",
                    "rid",
                    "rtoe_seconds",
                    "query_count",
                    "total_query_bases",
                    "query_names",
                    "submitted_at_utc",
                ]
            )
        writer.writerow(
            [
                result.batch_index,
                result.submission.rid,
                result.submission.rtoe_seconds,
                len(result.queries),
                sum(len(normalize_rna(query.sequence_5to3)) for query in result.queries),
                ";".join(query.name for query in result.queries),
                datetime.now(timezone.utc).isoformat(),
            ]
        )


def run_blast_batches(
    args: argparse.Namespace,
    queries: list[AntisenseQuery],
    *,
    client_factory: Callable[..., NcbiBlastClient] | None = None,
    sleeper: Callable[[float], None] | None = None,
) -> list[BlastBatchResult]:
    queries = assign_unique_blast_query_ids(queries)
    print(
        "Privacy warning: remote NCBI BLAST will transmit "
        f"{len(queries)} oligo sequence(s) outside this computer.",
        file=sys.stderr,
    )
    make_client = client_factory or NcbiBlastClient
    sleep = sleeper or time.sleep
    client = make_client(
        email=require_email(args.email),
        tool=args.tool,
        request_seconds=max(args.request_seconds, DEFAULT_REQUEST_SECONDS),
    )
    outputs = []
    batches = batch_antisense_queries(queries, args.max_batch_bases)
    for batch_index, batch in enumerate(batches, start=1):
        print(
            f"Submitting BLAST batch {batch_index}/{len(batches)} "
            f"({len(batch)} oligo sequences)...",
            file=sys.stderr,
        )
        submission = client.submit_blastn(
            query_fasta=multi_fasta(batch),
            database=args.database,
            expect=args.expect,
            word_size=args.word_size,
            hitlist_size=args.hitlist_size,
            megablast=args.megablast,
        )
        submitted_result = BlastBatchResult(
            batch_index=batch_index,
            submission=submission,
            queries=tuple(batch),
            csv_text="",
        )
        if args.rid_log:
            append_rid_log(args.rid_log, submitted_result)
        if submission.rtoe_seconds:
            sleep(max(submission.rtoe_seconds, DEFAULT_REQUEST_SECONDS))
        client.wait_for_result(
            submission.rid,
            poll_seconds=max(args.poll_seconds, DEFAULT_POLL_SECONDS),
            timeout_seconds=args.timeout_seconds,
        )
        result = BlastBatchResult(
            batch_index=batch_index,
            submission=submission,
            queries=tuple(batch),
            csv_text=client.fetch_csv(submission.rid, alignments=args.hitlist_size),
        )
        outputs.append(result)
    return outputs


def run_private_panel_workflow(
    args: argparse.Namespace,
    started_at: str,
    *,
    progress_callback: Callable[[int, int, str, str], None] | None = None,
    cancel_check: Callable[[], bool] | None = None,
    include_comparison_results: bool = False,
) -> int:
    """Retrieve public references and scan private guides entirely locally."""
    accessions = panel_accessions_from_args(args)
    cache_dir = private_panel_cache_dir(args)
    print(
        "Private local panel mode: guide sequences remain on this computer; "
        "NCBI EFetch requests contain transcript accessions only.",
        file=sys.stderr,
    )
    targets = retrieve_transcript_targets(
        accessions,
        email=args.email,
        tool=args.tool,
        cache_dir=cache_dir,
        offline=args.offline,
        refresh=bool(getattr(args, "refresh_targets", False)),
        request_seconds=args.request_seconds,
        progress_callback=progress_callback,
        cancel_check=cancel_check,
    )

    queries: list[AntisenseQuery] = []
    scan_regions = parse_scan_regions(args.scan_region)
    if args.download_targets_only:
        panel_result = PrivatePanelScanResult(
            targets=tuple(targets),
            matches=(),
            summaries=(),
        )
    else:
        queries = args_antisense_queries(args)
        panel_result = run_private_panel_scan(
            queries,
            targets,
            scan_regions,
            args.max_mismatches,
            closest=args.closest,
        )

    local_matches = list(panel_result.matches)
    if args.output:
        write_text(args.output, transcript_matches_to_csv(local_matches))
        print(f"Wrote private local panel matches to: {args.output}")
    if args.stdout_csv:
        print(transcript_matches_to_csv(local_matches), end="")

    result_workbook = default_private_panel_workbook(args)
    completed_at = datetime.now(timezone.utc).isoformat()
    write_result_workbook(
        result_workbook,
        args,
        queries,
        scan_regions,
        local_matches,
        [],
        started_at,
        completed_at,
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


def main() -> int:
    args = build_parser().parse_args()
    if args.gui:
        return run_gui()

    started_at = datetime.now(timezone.utc).isoformat()
    try:
        validate_runtime_args(args)
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
                all_local_matches = run_local_scan(args, queries, scan_regions, max_mismatches=None)
                local_matches = [
                    match
                    for match in all_local_matches
                    if match.mismatches <= args.max_mismatches
                ]
                closest_limit = args.closest or DEFAULT_CLOSEST_MATCHES
                closest_matches = closest_transcript_matches(all_local_matches, closest_limit)
            else:
                local_matches = run_local_scan(args, queries, scan_regions)
            csv_text = transcript_matches_to_csv(local_matches)
            show_closest_matches = args.closest is not None or (print_terminal and not local_matches)
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
                print("No local transcript matches found within mismatch threshold.", file=sys.stderr)

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
            completed_at = datetime.now(timezone.utc).isoformat()
            write_result_workbook(
                result_workbook,
                args,
                queries,
                scan_regions,
                local_matches,
                blast_outputs,
                started_at,
                completed_at,
                include_blast_sheets=args.blast,
            )
            workflow_label = "NCBI BLAST" if args.blast else "local transcript scan"
            print(f"Wrote {workflow_label} result workbook to: {result_workbook}")

        return 0
    except (ValueError, TimeoutError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
