"""Domain-independent FASTA structure parsing and formatting."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable


@dataclass(frozen=True)
class FastaRecord:
    """One FASTA record without an imposed nucleotide or protein alphabet."""

    identifier: str
    description: str
    sequence: str

    @property
    def header(self) -> str:
        """Return the complete header without the leading ``>``."""
        if self.description:
            return f"{self.identifier} {self.description}"
        return self.identifier

    @classmethod
    def from_header(cls, header: str, sequence: str) -> "FastaRecord":
        """Build a record by splitting a complete header at its first space."""
        cleaned_header = str(header).strip()
        if not cleaned_header:
            raise ValueError("FASTA header cannot be blank.")
        identifier, _, description = cleaned_header.partition(" ")
        return cls(
            identifier=identifier,
            description=description.strip(),
            sequence=str(sequence),
        )


def parse_fasta(text: object, *, ignore_comments: bool = True) -> list[FastaRecord]:
    """Parse FASTA structure while preserving every record boundary.

    Alphabet normalization is deliberately left to the calling workflow.
    Blank lines are ignored. Lines beginning with ``;`` are treated as FASTA
    comments by default.
    """
    records: list[FastaRecord] = []
    header: str | None = None
    sequence_lines: list[str] = []

    def append_record() -> None:
        if header is None:
            return
        records.append(FastaRecord.from_header(header, "".join(sequence_lines)))

    for raw_line in str(text).splitlines():
        line = raw_line.strip()
        if not line or (ignore_comments and line.startswith(";")):
            continue
        if line.startswith(">"):
            append_record()
            header = line[1:].strip()
            if not header:
                raise ValueError("FASTA header cannot be blank.")
            sequence_lines = []
        else:
            if header is None:
                raise ValueError(
                    "Sequence content appeared before the first FASTA header."
                )
            sequence_lines.append(line)

    append_record()
    if not records:
        raise ValueError("No FASTA records found.")
    return records


def require_single_fasta_record(
    records: Iterable[FastaRecord],
    *,
    source_label: str = "FASTA input",
) -> FastaRecord:
    """Return the only FASTA record or raise with the observed record count."""
    materialized = list(records)
    if len(materialized) != 1:
        raise ValueError(
            f"{source_label} must contain exactly one FASTA record; "
            f"found {len(materialized)}."
        )
    return materialized[0]


def format_fasta(
    record: FastaRecord,
    width: int = 70,
    *,
    trailing_newline: bool = True,
) -> str:
    """Format one FASTA record without changing its sequence alphabet."""
    if isinstance(width, bool) or not isinstance(width, int) or width < 1:
        raise ValueError("FASTA line width must be a positive integer.")
    if not record.identifier.strip():
        raise ValueError("FASTA identifier cannot be blank.")

    sequence_lines = [
        record.sequence[index : index + width]
        for index in range(0, len(record.sequence), width)
    ]
    formatted = f">{record.header}\n" + "\n".join(sequence_lines)
    if trailing_newline:
        formatted += "\n"
    return formatted
