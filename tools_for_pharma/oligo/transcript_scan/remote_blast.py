"""Explicit remote NCBI BLAST submission workflow.

Calling :func:`run_blast_batches` transmits the supplied guide sequences to
NCBI.  Local transcript scanning does not import this module.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from datetime import datetime, timezone
import io
from pathlib import Path
import sys
import time
from typing import Callable, Iterable

from tools_for_pharma.oligo.ncbi_transport import (
    DEFAULT_DATABASE,
    DEFAULT_EXPECT,
    DEFAULT_HITLIST_SIZE,
    DEFAULT_REQUEST_SECONDS,
    DEFAULT_POLL_SECONDS,
    DEFAULT_TOOL,
    DEFAULT_WORD_SIZE,
    BlastSubmission,
    NcbiBlastClient,
    require_email,
)
from tools_for_pharma.oligo.transcript_scan.models import AntisenseQuery
from tools_for_pharma.oligo.transcript_scan.queries import (
    DEFAULT_BATCH_BASES,
    assign_unique_blast_query_ids,
    batch_antisense_queries,
    sanitize_fasta_name,
)
from tools_for_pharma.sequence.fasta import FastaRecord, format_fasta
from tools_for_pharma.sequence.nucleotides import normalize_dna, normalize_rna


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
class BlastBatchResult:
    """One completed BLAST batch and its returned CSV text."""

    batch_index: int
    submission: BlastSubmission
    queries: tuple[AntisenseQuery, ...]
    csv_text: str


@dataclass(frozen=True)
class RemoteBlastConfig:
    """Configuration that makes remote guide submission explicit."""

    email: str | None
    tool: str = DEFAULT_TOOL
    database: str = DEFAULT_DATABASE
    expect: float | str = DEFAULT_EXPECT
    word_size: int = DEFAULT_WORD_SIZE
    hitlist_size: int = DEFAULT_HITLIST_SIZE
    megablast: bool = False
    timeout_seconds: float = 1800
    max_batch_bases: int = DEFAULT_BATCH_BASES
    request_seconds: float = DEFAULT_REQUEST_SECONDS
    poll_seconds: float = DEFAULT_POLL_SECONDS
    rid_log: Path | None = None


MessageCallback = Callable[[str], None]


def parse_blast_csv(text: str) -> list[dict[str, str]]:
    """Parse NCBI tabular CSV BLAST output into dictionaries."""
    rows = []
    for row in csv.reader(io.StringIO(text)):
        if not row or row[0].startswith("#"):
            continue
        if len(row) == len(CSV_COLUMNS):
            rows.append(dict(zip(CSV_COLUMNS, row)))
    return rows


def filter_blast_rows(
    rows: Iterable[dict[str, object]],
    max_mismatches: int,
    max_gap_opens: int,
    min_alignment_fraction: float,
) -> list[dict[str, object]]:
    """Apply the configured technical filters to parsed remote BLAST rows."""
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


def fasta_record(name: str, sequence: str, line_width: int = 80) -> str:
    """Return one DNA FASTA record suitable for a remote BLAST request."""
    cleaned = normalize_dna(sequence)
    return format_fasta(
        FastaRecord(sanitize_fasta_name(name), "", cleaned),
        width=line_width,
        trailing_newline=False,
    )


def multi_fasta(records: Iterable[AntisenseQuery]) -> str:
    """Return multi-FASTA for a batch of remote BLAST queries."""
    prepared = assign_unique_blast_query_ids(records)
    return "\n".join(
        fasta_record(record.blast_query_id, record.sequence_5to3)
        for record in prepared
    )


def combine_blast_csv(outputs: Iterable[BlastBatchResult]) -> str:
    combined = io.StringIO()
    writer = csv.writer(combined, lineterminator="\n")
    writer.writerow(["rid", *CSV_COLUMNS])
    for result in outputs:
        for row in parse_blast_csv(result.csv_text):
            writer.writerow(
                [result.submission.rid, *[row[column] for column in CSV_COLUMNS]]
            )
    return combined.getvalue()


def append_rid_log(
    path: Path,
    result: BlastBatchResult,
    *,
    now: Callable[[], datetime] | None = None,
) -> None:
    """Persist a BLAST RID immediately after submission for recovery/audit."""
    path.parent.mkdir(parents=True, exist_ok=True)
    write_header = not path.exists()
    current_time = now or (lambda: datetime.now(timezone.utc))
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
                current_time().isoformat(),
            ]
        )


def _stderr_message(message: str) -> None:
    print(message, file=sys.stderr)


def run_blast_batches(
    config: RemoteBlastConfig,
    queries: list[AntisenseQuery],
    *,
    client_factory: Callable[..., NcbiBlastClient] | None = None,
    sleeper: Callable[[float], None] | None = None,
    message_callback: MessageCallback | None = None,
) -> list[BlastBatchResult]:
    """Submit, poll, and retrieve explicitly remote BLAST query batches."""
    prepared_queries = assign_unique_blast_query_ids(queries)
    emit = message_callback or _stderr_message
    emit(
        "Privacy warning: remote NCBI BLAST will transmit "
        f"{len(prepared_queries)} oligo sequence(s) outside this computer."
    )
    make_client = client_factory or NcbiBlastClient
    sleep = sleeper or time.sleep
    client = make_client(
        email=require_email(config.email),
        tool=config.tool,
        request_seconds=max(config.request_seconds, DEFAULT_REQUEST_SECONDS),
    )
    outputs: list[BlastBatchResult] = []
    batches = batch_antisense_queries(prepared_queries, config.max_batch_bases)
    for batch_index, batch in enumerate(batches, start=1):
        emit(
            f"Submitting BLAST batch {batch_index}/{len(batches)} "
            f"({len(batch)} oligo sequences)..."
        )
        submission = client.submit_blastn(
            query_fasta=multi_fasta(batch),
            database=config.database,
            expect=config.expect,
            word_size=config.word_size,
            hitlist_size=config.hitlist_size,
            megablast=config.megablast,
        )
        submitted_result = BlastBatchResult(
            batch_index=batch_index,
            submission=submission,
            queries=tuple(batch),
            csv_text="",
        )
        if config.rid_log:
            append_rid_log(config.rid_log, submitted_result)
        if submission.rtoe_seconds:
            sleep(max(submission.rtoe_seconds, DEFAULT_REQUEST_SECONDS))
        client.wait_for_result(
            submission.rid,
            poll_seconds=max(config.poll_seconds, DEFAULT_POLL_SECONDS),
            timeout_seconds=config.timeout_seconds,
        )
        outputs.append(
            BlastBatchResult(
                batch_index=batch_index,
                submission=submission,
                queries=tuple(batch),
                csv_text=client.fetch_csv(
                    submission.rid,
                    alignments=config.hitlist_size,
                ),
            )
        )
    return outputs
