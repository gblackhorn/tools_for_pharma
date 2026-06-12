"""NCBI transcript fetch and BLAST helpers for antisense oligos.

This module has two related workflows:

1. Specific transcript check:
   Fetch an NM/XM/NR/XR accession with NCBI EFetch, then scan the transcript for
   the reverse-complement target of an antisense sequence.

2. BLAST database search:
   Submit the antisense sequence to the NCBI BLAST URL API and retrieve a CSV
   report. This is best for broad searches such as refseq_rna or core_nt.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import csv
from datetime import datetime, timezone
import io
from pathlib import Path
import re
import sys
import time
from typing import Iterable
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from tools_for_pharma.oligo.core import get_complementary_sequence, normalize_rna
from tools_for_pharma.oligo.transcript import fasta_or_plain_text_to_sequence, get_fasta_header
from tools_for_pharma.shared.excel_utils import list_excel_sheets


BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
DEFAULT_TOOL = "tools_for_pharma_oligo"
DEFAULT_EMAIL = "da.guo@argobiopharma.com"
DEFAULT_DATABASE = "refseq_rna"
DEFAULT_PROGRAM = "blastn"
DEFAULT_EXPECT = "1000"
DEFAULT_WORD_SIZE = 7
DEFAULT_HITLIST_SIZE = 50
DEFAULT_MAX_MISMATCHES = 3
DEFAULT_BATCH_BASES = 1000
DEFAULT_POLL_SECONDS = 75
DEFAULT_REQUEST_SECONDS = 15
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
    """One named antisense input sequence."""

    name: str
    sequence_5to3: str
    target_accession: str = ""
    target_gene: str = ""
    species: str = ""
    notes: str = ""


@dataclass(frozen=True)
class AntisenseRegion:
    """A 1-based inclusive AS subregion to scan."""

    name: str
    start: int | None = None
    end: int | None = None


@dataclass(frozen=True)
class TranscriptMatch:
    """One local antisense-vs-transcript match."""

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


@dataclass(frozen=True)
class BlastSubmission:
    """NCBI BLAST request metadata returned after CMD=Put."""

    rid: str
    rtoe_seconds: int | None


@dataclass(frozen=True)
class BlastBatchResult:
    """One completed BLAST batch and its returned CSV text."""

    batch_index: int
    submission: BlastSubmission
    queries: tuple[AntisenseQuery, ...]
    csv_text: str


def normalize_dna(sequence: str) -> str:
    """Normalize a sequence to DNA letters for NCBI BLAST requests."""
    return normalize_rna(sequence).replace("U", "T")


def fasta_record(name: str, sequence: str, line_width: int = 80) -> str:
    """Return a simple FASTA record from a raw nucleotide sequence."""
    cleaned = normalize_dna(sequence)
    lines = [f">{sanitize_fasta_name(name)}"]
    lines.extend(cleaned[index : index + line_width] for index in range(0, len(cleaned), line_width))
    return "\n".join(lines)


def sanitize_fasta_name(name: str) -> str:
    """Return a FASTA-safe query identifier."""
    cleaned = re.sub(r"[^A-Za-z0-9_.:-]+", "_", clean_text_for_id(name)).strip("_")
    return cleaned or "antisense_query"


def clean_text_for_id(value: object) -> str:
    return re.sub(r"\s+", " ", str(value)).strip()


def multi_fasta(records: Iterable[AntisenseQuery]) -> str:
    """Return a multi-FASTA string for one or more AS queries."""
    return "\n".join(fasta_record(record.name, record.sequence_5to3) for record in records)


def require_email(email: str | None) -> str:
    """Return a usable NCBI contact email or raise a clear error."""
    return email or DEFAULT_EMAIL


class NcbiHttpClient:
    """Small HTTP client with NCBI-friendly request spacing."""

    def __init__(
        self,
        email: str,
        tool: str = DEFAULT_TOOL,
        request_seconds: int = DEFAULT_REQUEST_SECONDS,
    ) -> None:
        self.email = require_email(email)
        self.tool = tool
        self.request_seconds = request_seconds
        self._last_request_time = 0.0

    def _wait_if_needed(self) -> None:
        elapsed = time.monotonic() - self._last_request_time
        wait_seconds = self.request_seconds - elapsed
        if wait_seconds > 0:
            time.sleep(wait_seconds)

    def get_text(self, url: str, params: dict[str, object]) -> str:
        self._wait_if_needed()
        clean_params = {
            key: value
            for key, value in params.items()
            if value is not None and value != ""
        }
        query = urlencode(clean_params)
        request = Request(
            f"{url}?{query}",
            headers={"User-Agent": f"{self.tool}/1.0 ({self.email})"},
        )
        try:
            with urlopen(request, timeout=120) as response:
                text = response.read().decode("utf-8", errors="replace")
        except HTTPError as error:
            raise ValueError(f"NCBI HTTP error {error.code}: {error.reason}") from error
        except URLError as error:
            raise ValueError(f"NCBI request failed: {error.reason}") from error
        finally:
            self._last_request_time = time.monotonic()
        return text


class NcbiBlastClient(NcbiHttpClient):
    """Client for the NCBI BLAST Common URL API."""

    def submit_blastn(
        self,
        query_sequence: str | None = None,
        query_fasta: str | None = None,
        database: str = DEFAULT_DATABASE,
        expect: str = DEFAULT_EXPECT,
        word_size: int = DEFAULT_WORD_SIZE,
        hitlist_size: int = DEFAULT_HITLIST_SIZE,
        megablast: bool = False,
        short_query_adjust: bool = True,
    ) -> BlastSubmission:
        if query_fasta is None:
            if query_sequence is None:
                raise ValueError("Provide query_sequence or query_fasta for BLAST submission.")
            query_fasta = fasta_record("antisense_query", query_sequence)
        params = {
            "CMD": "Put",
            "PROGRAM": DEFAULT_PROGRAM,
            "DATABASE": database,
            "QUERY": query_fasta,
            "EXPECT": expect,
            "WORD_SIZE": word_size,
            "HITLIST_SIZE": hitlist_size,
            "SHORT_QUERY_ADJUST": str(short_query_adjust).lower(),
            "FILTER": "F",
            "MEGABLAST": "on" if megablast else None,
            "tool": self.tool,
            "email": self.email,
        }
        text = self.get_text(BLAST_URL, params)
        rid = parse_blast_field(text, "RID")
        if not rid:
            raise ValueError(f"NCBI BLAST submission did not return an RID:\n{text[:500]}")
        rtoe = parse_blast_field(text, "RTOE")
        return BlastSubmission(rid=rid, rtoe_seconds=int(rtoe) if rtoe and rtoe.isdigit() else None)

    def blast_status(self, rid: str) -> str:
        text = self.get_text(
            BLAST_URL,
            {
                "CMD": "Get",
                "RID": rid,
                "FORMAT_OBJECT": "SearchInfo",
                "tool": self.tool,
                "email": self.email,
            },
        )
        status = parse_blast_field(text, "Status")
        return status or "UNKNOWN"

    def wait_for_result(
        self,
        rid: str,
        poll_seconds: int = DEFAULT_POLL_SECONDS,
        timeout_seconds: int = 1800,
    ) -> None:
        deadline = time.monotonic() + timeout_seconds
        while time.monotonic() < deadline:
            status = self.blast_status(rid)
            if status == "READY":
                return
            if status in {"FAILED", "UNKNOWN", "EXPIRED"}:
                raise ValueError(f"NCBI BLAST RID {rid} returned status {status}.")
            time.sleep(max(poll_seconds, DEFAULT_POLL_SECONDS))
        raise TimeoutError(f"Timed out waiting for NCBI BLAST RID {rid}.")

    def fetch_csv(self, rid: str, alignments: int = DEFAULT_HITLIST_SIZE) -> str:
        return self.get_text(
            BLAST_URL,
            {
                "CMD": "Get",
                "RID": rid,
                "FORMAT_TYPE": "CSV",
                "ALIGNMENT_VIEW": "Tabular",
                "ALIGNMENTS": alignments,
                "DESCRIPTIONS": alignments,
                "tool": self.tool,
                "email": self.email,
            },
        )

    def run_blastn(
        self,
        query_sequence: str | None = None,
        query_fasta: str | None = None,
        database: str = DEFAULT_DATABASE,
        expect: str = DEFAULT_EXPECT,
        word_size: int = DEFAULT_WORD_SIZE,
        hitlist_size: int = DEFAULT_HITLIST_SIZE,
        megablast: bool = False,
        timeout_seconds: int = 1800,
    ) -> tuple[BlastSubmission, str]:
        submission = self.submit_blastn(
            query_sequence=query_sequence,
            query_fasta=query_fasta,
            database=database,
            expect=expect,
            word_size=word_size,
            hitlist_size=hitlist_size,
            megablast=megablast,
        )
        if submission.rtoe_seconds:
            time.sleep(max(submission.rtoe_seconds, DEFAULT_REQUEST_SECONDS))
        self.wait_for_result(submission.rid, timeout_seconds=timeout_seconds)
        return submission, self.fetch_csv(submission.rid, alignments=hitlist_size)


def parse_blast_field(text: str, field_name: str) -> str | None:
    """Parse BLAST API fields such as RID, RTOE, or Status."""
    pattern = re.compile(rf"^\s*{re.escape(field_name)}\s*=\s*(\S+)\s*$", re.MULTILINE)
    match = pattern.search(text)
    return match.group(1) if match else None


def parse_fasta_records(text: str) -> list[AntisenseQuery]:
    """Parse FASTA text into AS query records."""
    records = []
    name: str | None = None
    sequence_lines: list[str] = []
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append(AntisenseQuery(name, normalize_rna("".join(sequence_lines))))
            name = line[1:].strip() or f"AS_{len(records) + 1}"
            sequence_lines = []
        else:
            sequence_lines.append(line)

    if name is not None:
        records.append(AntisenseQuery(name, normalize_rna("".join(sequence_lines))))
    return records


def parse_plain_antisense_lines(text: str) -> list[AntisenseQuery]:
    """Parse a plain text list of AS sequences.

    Accepted line styles:
      AUGCUA...
      AS_001,AUGCUA...
      AS_001<TAB>AUGCUA...
      AS_001 AUGCUA...
    """
    records = []
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
            name = f"AS_{len(records) + 1}"
            sequence = parts[0]
        records.append(AntisenseQuery(name, normalize_rna(sequence)))
    return records


def read_antisense_file(path: Path) -> list[AntisenseQuery]:
    """Read AS queries from FASTA or plain text."""
    text = path.read_text(encoding="utf-8-sig")
    if any(line.lstrip().startswith(">") for line in text.splitlines()):
        records = parse_fasta_records(text)
    else:
        records = parse_plain_antisense_lines(text)
    if not records:
        raise ValueError(f"No AS sequences found in {path}.")
    return records


def read_antisense_table(
    path: Path,
    sequence_column: str | None = None,
    name_column: str | None = None,
    target_accession_column: str | None = None,
    sheet_name: str | int | None = None,
) -> list[AntisenseQuery]:
    """Read AS queries from an Excel or CSV table."""
    import pandas as pd

    suffix = path.suffix.lower()
    if suffix in {".xlsx", ".xlsm", ".xls"}:
        table = pd.read_excel(path, sheet_name=sheet_name or 0)
    elif suffix in {".csv", ".txt"}:
        table = pd.read_csv(path)
    else:
        raise ValueError("AS table must be an Excel workbook or CSV/text file.")

    if table.empty:
        raise ValueError(f"AS table is empty: {path}")

    columns_by_lower = {clean_text_for_id(column).lower(): str(column) for column in table.columns}
    if sequence_column is None:
        for candidate in ["antisense", "as", "as_sequence", "as sequence", "sequence"]:
            if candidate in columns_by_lower:
                sequence_column = columns_by_lower[candidate]
                break
        if sequence_column is None:
            sequence_column = str(table.columns[0])
    if sequence_column not in table.columns:
        raise ValueError(f"AS table is missing sequence column: {sequence_column}")

    if name_column is None:
        for candidate in ["name", "id", "oligo", "oligo id", "oligo_id", "as name"]:
            if candidate in columns_by_lower:
                name_column = columns_by_lower[candidate]
                break
    elif name_column not in table.columns:
        raise ValueError(f"AS table is missing name column: {name_column}")

    if target_accession_column is None:
        for candidate in ["target_accession", "target accession", "refseq", "refseq accession"]:
            if candidate in columns_by_lower:
                target_accession_column = columns_by_lower[candidate]
                break
    elif target_accession_column not in table.columns:
        raise ValueError(f"AS table is missing target accession column: {target_accession_column}")

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

        records.append(
            AntisenseQuery(
                name=name or f"AS_{row_index + 1}",
                sequence_5to3=sequence,
                target_accession=target_accession,
                target_gene=metadata["target_gene"],
                species=metadata["species"],
                notes=metadata["notes"],
            )
        )

    if not records:
        raise ValueError(f"No AS sequences found in table: {path}")
    return records


def read_antisense_queries(
    as_sequence: str | None = None,
    as_name: str | None = None,
    as_file: Path | None = None,
    as_table: Path | None = None,
    as_column: str | None = None,
    as_name_column: str | None = None,
    target_accession_column: str | None = None,
    as_sheet: str | int | None = None,
) -> list[AntisenseQuery]:
    """Read AS queries from exactly one input source."""
    provided = [as_sequence is not None, as_file is not None, as_table is not None]
    if sum(provided) != 1:
        raise ValueError("Provide exactly one of --as-sequence, --as-file, or --as-table.")

    if as_sequence is not None:
        return [AntisenseQuery(as_name or "antisense_query", normalize_rna(as_sequence))]
    if as_file is not None:
        return read_antisense_file(as_file)
    assert as_table is not None
    return read_antisense_table(
        as_table,
        as_column,
        as_name_column,
        target_accession_column,
        as_sheet,
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
    duplicate_groups = duplicate_sequence_groups(records)
    rows = []
    for index, record in enumerate(records, start=1):
        sequence = normalize_rna(record.sequence_5to3)
        duplicate_names = duplicate_groups.get(sequence, [])
        rows.append(
            {
                "input_order": index,
                "antisense_name": record.name,
                "blast_query_id": sanitize_fasta_name(record.name),
                "antisense_5to3": sequence,
                "length_nt": len(sequence),
                "target_accession": record.target_accession,
                "target_gene": record.target_gene,
                "species": record.species,
                "notes": record.notes,
                "is_duplicate_sequence": bool(duplicate_names),
                "duplicate_group_names": ";".join(duplicate_names),
            }
        )
    return rows


def batch_antisense_queries(
    records: list[AntisenseQuery],
    max_batch_bases: int = DEFAULT_BATCH_BASES,
) -> list[list[AntisenseQuery]]:
    """Group short AS queries into multi-FASTA BLAST batches."""
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


def fetch_transcript_fasta(
    accession: str,
    email: str,
    tool: str = DEFAULT_TOOL,
    cache_dir: Path | None = None,
) -> str:
    """Fetch transcript FASTA from NCBI EFetch by accession or UID."""
    if cache_dir is not None:
        cache_path = transcript_cache_path(cache_dir, accession)
        if cache_path.exists():
            return cache_path.read_text(encoding="utf-8-sig")

    client = NcbiHttpClient(email=email, tool=tool)
    text = client.get_text(
        EFETCH_URL,
        {
            "db": "nuccore",
            "id": accession,
            "rettype": "fasta",
            "retmode": "text",
            "tool": tool,
            "email": email,
        },
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
            email=require_email(email),
            tool=tool,
            cache_dir=cache_dir,
        )
        return get_fasta_header(fasta_text) or accession, fasta_or_plain_text_to_sequence(fasta_text)

    if transcript_file:
        text = transcript_file.read_text(encoding="utf-8-sig")
        return get_fasta_header(text) or transcript_file.name, fasta_or_plain_text_to_sequence(text)

    assert transcript_sequence is not None
    return get_fasta_header(transcript_sequence) or "target_transcript", fasta_or_plain_text_to_sequence(transcript_sequence)


def mismatch_positions(query: str, target: str) -> tuple[int, ...]:
    """Return 1-based mismatch positions between equal-length RNA strings."""
    return tuple(index + 1 for index, (left, right) in enumerate(zip(query, target)) if left != right)


def scan_antisense_against_transcript(
    antisense_5to3: str,
    transcript_sequence: str,
    transcript_name: str = "target_transcript",
    antisense_name: str = "antisense_query",
    scan_region: AntisenseRegion | None = None,
    max_mismatches: int = DEFAULT_MAX_MISMATCHES,
) -> list[TranscriptMatch]:
    """Find reverse-complement antisense target windows in a transcript.

    The transcript window is reported in transcript 5'->3' orientation. The AS
    oligo is first reverse-complemented to the expected transcript target.
    """
    antisense = normalize_rna(antisense_5to3)
    transcript = normalize_rna(transcript_sequence)
    region = scan_region or AntisenseRegion("full")
    region_sequence, region_start, region_end = antisense_region_sequence(antisense, region)
    target = get_complementary_sequence(region_sequence, reverse=True)
    if len(transcript) < len(target):
        return []

    matches = []
    for start_index in range(0, len(transcript) - len(target) + 1):
        window = transcript[start_index : start_index + len(target)]
        mismatches = mismatch_positions(target, window)
        if len(mismatches) <= max_mismatches:
            transcript_match_as = get_complementary_sequence(window, reverse=True)
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
                )
            )
    return sorted(matches, key=lambda item: (item.mismatches, item.transcript_start))


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
            ]
        )
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


def transcript_match_rows(matches: Iterable[TranscriptMatch]) -> list[dict[str, object]]:
    rows = []
    for match in matches:
        rows.append(
            {
                "transcript_name": match.transcript_name,
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
                "transcript_match_as_5to3": match.transcript_match_as_5to3,
                "mismatch_positions_1based": ";".join(
                    str(position) for position in match.mismatch_positions_1based
                ),
                "as_mismatch_positions_1based": ";".join(
                    str(position) for position in match.as_mismatch_positions_1based
                ),
            }
        )
    return rows


def query_length_by_blast_id(queries: Iterable[AntisenseQuery]) -> dict[str, int]:
    return {
        sanitize_fasta_name(query.name): len(normalize_rna(query.sequence_5to3))
        for query in queries
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
) -> list[dict[str, object]]:
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
        "total_query_bases": sum(len(normalize_rna(query.sequence_5to3)) for query in queries),
        "scan_regions": ";".join(region.name for region in scan_regions),
        "max_mismatches_local_scan": args.max_mismatches,
        "blast_filter_max_mismatches": args.filter_max_mismatches,
        "blast_filter_max_gap_opens": args.filter_max_gap_opens,
        "blast_filter_min_alignment_fraction": args.filter_min_alignment_fraction,
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
    source = args.as_table or args.as_file
    if source and not args.output and not args.blast_output:
        return source.with_name(f"{source.stem}_ncbi_blast_results.xlsx")
    return None


def default_gui_result_workbook(input_file: Path) -> Path:
    return input_file.with_name(f"{input_file.stem}_ncbi_blast_results.xlsx")


def write_result_workbook(
    path: Path,
    args: argparse.Namespace,
    queries: list[AntisenseQuery],
    scan_regions: list[AntisenseRegion],
    local_matches: list[TranscriptMatch],
    blast_results: list[BlastBatchResult],
    started_at: str,
    completed_at: str,
) -> None:
    raw_blast_rows = blast_raw_rows(blast_results, queries)
    sheets = {
        "input_queries": input_query_rows(queries),
        "local_transcript_scan": transcript_match_rows(local_matches),
        "blast_hits_raw": raw_blast_rows,
        "blast_hits_filtered": filter_blast_rows(
            raw_blast_rows,
            args.filter_max_mismatches,
            args.filter_max_gap_opens,
            args.filter_min_alignment_fraction,
        ),
        "blast_batches": blast_batch_rows(blast_results),
        "run_metadata": metadata_rows(args, queries, scan_regions, started_at, completed_at),
    }
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
    window.title("Select AS table sheet")
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
    window.title("NCBI transcript scan settings")
    window.resizable(False, False)
    window.columnconfigure(1, weight=1)

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
    scan_regions_var = tk.StringVar(value="full")
    max_mismatches_var = tk.StringVar(value=str(DEFAULT_MAX_MISMATCHES))

    ttk.Label(window, text="AS sequence column").grid(
        row=0, column=0, padx=16, pady=(16, 8), sticky="w"
    )
    ttk.Combobox(
        window,
        textvariable=as_column_var,
        values=headers,
        state="readonly",
        width=max(30, min(60, max(len(header) for header in headers) + 2)),
    ).grid(row=0, column=1, padx=16, pady=(16, 8), sticky="ew")

    ttk.Label(window, text="AS name column").grid(
        row=1, column=0, padx=16, pady=8, sticky="w"
    )
    ttk.Combobox(
        window,
        textvariable=name_column_var,
        values=["", *headers],
        state="readonly",
        width=36,
    ).grid(row=1, column=1, padx=16, pady=8, sticky="ew")

    ttk.Label(window, text="Transcript source").grid(
        row=2, column=0, padx=16, pady=8, sticky="nw"
    )
    source_frame = ttk.Frame(window)
    source_frame.grid(row=2, column=1, padx=16, pady=8, sticky="ew")
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

    ttk.Label(window, text="Scan regions").grid(
        row=3, column=0, padx=16, pady=8, sticky="w"
    )
    ttk.Entry(window, textvariable=scan_regions_var, width=38).grid(
        row=3, column=1, padx=16, pady=8, sticky="ew"
    )

    ttk.Label(window, text="Max mismatches").grid(
        row=4, column=0, padx=16, pady=8, sticky="w"
    )
    ttk.Entry(window, textvariable=max_mismatches_var, width=12).grid(
        row=4, column=1, padx=16, pady=8, sticky="w"
    )

    buttons = ttk.Frame(window)
    buttons.grid(row=5, column=0, columnspan=2, padx=16, pady=(8, 16), sticky="e")

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

        selected.update(
            {
                "as_column": as_column_var.get(),
                "as_name_column": name_column_var.get() or None,
                "target_mode": mode,
                "target_accession_column": target_column_var.get() if mode == "column" else None,
                "target_accession": refseq_var.get().strip() if mode == "refseq" else None,
                "target_file": Path(transcript_file_var.get()) if mode == "file" else None,
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
    output_path = default_gui_result_workbook(input_file)
    return argparse.Namespace(
        as_sequence=None,
        as_name=None,
        as_file=None,
        as_table=input_file,
        as_column=settings["as_column"],
        as_name_column=settings["as_name_column"],
        as_sheet=sheet_name,
        target_accession=settings["target_accession"],
        target_accession_column=settings["target_accession_column"],
        target_file=settings["target_file"],
        target_sequence=None,
        scan_region=settings["scan_regions"],
        max_mismatches=settings["max_mismatches"],
        email=DEFAULT_EMAIL,
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
        cache_dir=input_file.with_name(f"{input_file.stem}_ncbi_cache"),
        rid_log=None,
        output=None,
        blast_output=None,
        result_workbook=output_path,
        gui=False,
    )


def run_gui() -> int:
    import tkinter as tk
    from tkinter import filedialog, messagebox

    root = tk.Tk()
    root.withdraw()

    try:
        input_path = filedialog.askopenfilename(
            title="Select AS Excel input file",
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

        args = gui_args(input_file, sheet_name, settings)
        started_at = datetime.now(timezone.utc).isoformat()
        queries = args_antisense_queries(args)
        scan_regions = parse_scan_regions(args.scan_region)
        local_matches = run_local_scan(args, queries, scan_regions)
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
        )
        messagebox.showinfo(
            "Done",
            f"Wrote NCBI transcript scan workbook to:\n{args.result_workbook}",
        )
        return 0
    except Exception as error:
        messagebox.showerror("NCBI transcript scan failed", str(error))
        return 1
    finally:
        root.destroy()


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compare AS oligos to transcripts and optionally run NCBI BLAST."
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
    parser.add_argument("--target-accession", help="NM/XM/NR/XR accession to fetch from NCBI.")
    parser.add_argument(
        "--target-accession-column",
        help=(
            "Column in --as-table containing per-row NM/XM/NR/XR accessions. "
            "Defaults to target_accession if present."
        ),
    )
    parser.add_argument("--target-file", type=Path, help="FASTA/plain transcript file.")
    parser.add_argument("--target-sequence", help="Pasted FASTA/plain transcript sequence.")
    parser.add_argument(
        "--scan-region",
        action="append",
        help=(
            "AS region to scan locally: full, 2-18, or seed:2-8. "
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
            "Contact email for NCBI API usage guidelines. "
            f"Defaults to {DEFAULT_EMAIL}."
        ),
    )
    parser.add_argument("--tool", default=DEFAULT_TOOL, help=f"NCBI tool name. Defaults to {DEFAULT_TOOL}.")
    parser.add_argument(
        "--blast",
        action="store_true",
        help="Also run NCBI BLAST URL API against a nucleotide database.",
    )
    parser.add_argument(
        "--blast-only",
        action="store_true",
        help="Run BLAST without requiring a specific target transcript.",
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
    parser.add_argument("--megablast", action="store_true", help="Enable megablast for near-identical hits.")
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
            "Maximum total AS bases per BLAST multi-FASTA submission. "
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
        help="Optional folder for cached NCBI EFetch transcript FASTA files.",
    )
    parser.add_argument(
        "--rid-log",
        type=Path,
        help="Optional CSV log written as soon as each BLAST RID is submitted.",
    )
    parser.add_argument("--output", type=Path, help="Write local scan CSV to this path.")
    parser.add_argument("--blast-output", type=Path, help="Write BLAST CSV to this path.")
    parser.add_argument(
        "--result-workbook",
        type=Path,
        help=(
            "Write an Excel workbook with input_queries, local scan, BLAST hits, "
            "batch metadata, and run metadata. For --as-file/--as-table, defaults "
            "to '<input>_ncbi_blast_results.xlsx' when no CSV output is requested."
        ),
    )
    parser.add_argument(
        "--gui",
        action="store_true",
        help="Choose AS Excel input and transcript target source with dialogs.",
    )
    return parser


def args_antisense_queries(args: argparse.Namespace) -> list[AntisenseQuery]:
    return read_antisense_queries(
        as_sequence=args.as_sequence,
        as_name=args.as_name,
        as_file=args.as_file,
        as_table=args.as_table,
        as_column=args.as_column,
        as_name_column=args.as_name_column,
        target_accession_column=args.target_accession_column,
        as_sheet=args.as_sheet,
    )


def run_local_scan(
    args: argparse.Namespace,
    queries: list[AntisenseQuery],
    scan_regions: list[AntisenseRegion],
) -> list[TranscriptMatch]:
    matches = []
    shared_transcript: tuple[str, str] | None = None
    if not args.target_accession_column:
        shared_transcript = read_transcript_input(
            transcript_sequence=args.target_sequence,
            transcript_file=args.target_file,
            accession=args.target_accession,
            email=args.email,
            tool=args.tool,
            cache_dir=args.cache_dir,
        )

    for query in queries:
        if args.target_accession_column:
            if not query.target_accession:
                continue
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
                    max_mismatches=args.max_mismatches,
                )
            )
    return matches


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
) -> list[BlastBatchResult]:
    client = NcbiBlastClient(
        email=require_email(args.email),
        tool=args.tool,
        request_seconds=max(args.request_seconds, DEFAULT_REQUEST_SECONDS),
    )
    outputs = []
    batches = batch_antisense_queries(queries, args.max_batch_bases)
    for batch_index, batch in enumerate(batches, start=1):
        print(
            f"Submitting BLAST batch {batch_index}/{len(batches)} "
            f"({len(batch)} AS sequences)...",
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
            time.sleep(max(submission.rtoe_seconds, DEFAULT_REQUEST_SECONDS))
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


def main() -> int:
    args = build_parser().parse_args()
    if args.gui:
        return run_gui()

    started_at = datetime.now(timezone.utc).isoformat()
    try:
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
        blast_outputs: list[BlastBatchResult] = []
        if args.blast_only:
            args.blast = True
        if not args.blast_only:
            local_matches = run_local_scan(args, queries, scan_regions)
            csv_text = transcript_matches_to_csv(local_matches)
            if args.output:
                write_text(args.output, csv_text)
                print(f"Wrote local transcript scan to: {args.output}")
            elif not default_result_workbook(args):
                print(csv_text, end="")
            if not local_matches:
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
            )
            print(f"Wrote NCBI oligo result workbook to: {result_workbook}")

        return 0
    except (ValueError, TimeoutError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
