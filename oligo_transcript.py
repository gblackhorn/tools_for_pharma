"""Transcript FASTA/plain-text helpers for oligo design."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from oligo_core import get_complementary_sequence, get_subsequence, normalize_rna


@dataclass(frozen=True)
class TranscriptOligoResult:
    """A transcript range formatted as sense and antisense oligo strands."""

    source_name: str
    transcript_length: int
    start: int
    end: int
    selected_length: int
    ss_5to3: str
    as_5to3: str


def read_text_sequence_file(input_path: str) -> str:
    """Read a text-like sequence file."""
    source = Path(input_path)
    if not source.exists():
        raise ValueError(f"Input file does not exist: {source}")
    if not source.is_file():
        raise ValueError(f"Input path is not a file: {source}")

    return source.read_text(encoding="utf-8-sig")


def fasta_or_plain_text_to_sequence(text: str) -> str:
    """Extract sequence letters from FASTA or plain sequence text as RNA."""
    lines = []
    for raw_line in str(text).splitlines():
        line = raw_line.strip()
        if not line or line.startswith(">") or line.startswith(";"):
            continue
        lines.append(line)

    return normalize_rna("".join(lines))


def get_fasta_header(text: str) -> Optional[str]:
    """Return the first FASTA header without the leading '>', if present."""
    for raw_line in str(text).splitlines():
        line = raw_line.strip()
        if line.startswith(">"):
            return line[1:].strip() or None
    return None


def transcript_region_to_oligo(
    transcript_sequence: str,
    start: int,
    end: int,
    source_name: str = "transcript",
) -> TranscriptOligoResult:
    """Return transcript sense strand and antisense reverse-complement.

    Coordinates are 1-based and inclusive, matching common NCBI transcript
    position numbering.
    """
    normalized_transcript = normalize_rna(transcript_sequence)
    ss_5to3 = get_subsequence(normalized_transcript, start=start, end=end)
    as_5to3 = get_complementary_sequence(ss_5to3, reverse=True)

    return TranscriptOligoResult(
        source_name=source_name,
        transcript_length=len(normalized_transcript),
        start=start,
        end=end,
        selected_length=len(ss_5to3),
        ss_5to3=ss_5to3,
        as_5to3=as_5to3,
    )


def transcript_text_to_oligo(
    text: str,
    start: int,
    end: int,
    source_name: Optional[str] = None,
) -> TranscriptOligoResult:
    """Parse FASTA/plain text and return SS/AS oligo strands."""
    header = get_fasta_header(text)
    name = source_name or header or "transcript"
    sequence = fasta_or_plain_text_to_sequence(text)
    return transcript_region_to_oligo(sequence, start=start, end=end, source_name=name)


def transcript_file_to_oligo(
    input_path: str,
    start: int,
    end: int,
) -> TranscriptOligoResult:
    """Read a FASTA/plain-text transcript file and return SS/AS oligo strands."""
    source = Path(input_path)
    text = read_text_sequence_file(str(source))
    return transcript_text_to_oligo(
        text,
        start=start,
        end=end,
        source_name=get_fasta_header(text) or source.name,
    )


def format_transcript_oligo_result(result: TranscriptOligoResult) -> str:
    """Format transcript SS/AS output for CLI, GUI, or text files."""
    return "\n".join(
        [
            f"source: {result.source_name}",
            f"transcript_length_nt: {result.transcript_length}",
            f"range_1based_inclusive: {result.start}-{result.end}",
            f"selected_length_nt: {result.selected_length}",
            f"SS_5to3: {result.ss_5to3}",
            f"AS_5to3: {result.as_5to3}",
        ]
    )
