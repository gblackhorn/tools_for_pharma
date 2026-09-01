"""Pure AS/SS query preparation and scan-region parsing."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import re
from typing import Iterable

from tools_for_pharma.oligo.transcript_scan.models import (
    AntisenseQuery,
    AntisenseRegion,
)
from tools_for_pharma.sequence.fasta import parse_fasta
from tools_for_pharma.sequence.nucleotides import normalize_rna


DEFAULT_BATCH_BASES = 1000


def clean_text_for_id(value: object) -> str:
    return re.sub(r"\s+", " ", str(value)).strip()


def sanitize_fasta_name(name: str) -> str:
    """Return a FASTA-safe query identifier."""
    cleaned = re.sub(
        r"[^A-Za-z0-9_.:-]+",
        "_",
        clean_text_for_id(name),
    ).strip("_")
    return cleaned or "oligo_query"


def assign_unique_blast_query_ids(
    records: Iterable[AntisenseQuery],
) -> list[AntisenseQuery]:
    """Return records with stable, unique FASTA identifiers."""
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


def parse_fasta_records(
    text: str,
    sequence_type: str = "AS",
) -> list[AntisenseQuery]:
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


def parse_plain_antisense_lines(
    text: str,
    sequence_type: str = "AS",
) -> list[AntisenseQuery]:
    """Parse named or unnamed AS/SS sequences from plain text lines."""
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
        records.append(
            AntisenseQuery(
                name,
                normalize_rna(sequence),
                sequence_type=normalized_type,
            )
        )
    return records


def read_antisense_file(
    path: Path,
    sequence_type: str = "AS",
) -> list[AntisenseQuery]:
    """Read AS or SS queries from FASTA or plain text."""
    text = path.read_text(encoding="utf-8-sig")
    if any(line.lstrip().startswith(">") for line in text.splitlines()):
        records = parse_fasta_records(text, sequence_type=sequence_type)
    else:
        records = parse_plain_antisense_lines(text, sequence_type=sequence_type)
    if not records:
        normalized_type = normalize_sequence_type(sequence_type)
        raise ValueError(f"No {normalized_type} sequences found in {path}.")
    return records


def duplicate_sequence_groups(
    records: list[AntisenseQuery],
) -> dict[str, list[str]]:
    """Return normalized sequences and names for duplicate groups only."""
    groups: dict[str, list[str]] = {}
    for record in records:
        groups.setdefault(normalize_rna(record.sequence_5to3), []).append(record.name)
    return {
        sequence: names
        for sequence, names in groups.items()
        if len(names) > 1
    }


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


def antisense_region_sequence(
    sequence: str,
    region: AntisenseRegion,
) -> tuple[str, int, int]:
    antisense = normalize_rna(sequence)
    if region.start is None or region.end is None:
        return antisense, 1, len(antisense)
    if region.end > len(antisense):
        raise ValueError(
            f"Scan region {region.name} ends at {region.end}, but AS sequence "
            f"is only {len(antisense)} nt."
        )
    return antisense[region.start - 1 : region.end], region.start, region.end
