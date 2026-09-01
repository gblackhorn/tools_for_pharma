"""Validation and parsing helpers for versioned RefSeq transcripts."""

from __future__ import annotations

import re


VERSIONED_REFSEQ_TRANSCRIPT_RE = re.compile(
    r"^(?:NM|XM|NR|XR)_\d+\.\d+$",
    re.IGNORECASE,
)
VERSIONED_REFSEQ_GENOMIC_RE = re.compile(
    r"^(?:NC|NG|NT|NW)_\d+\.\d+$",
    re.IGNORECASE,
)


def _clean_identifier_text(value: object) -> str:
    return re.sub(r"\s+", " ", str(value)).strip()


def normalize_versioned_refseq_accession(accession: object) -> str:
    """Return an uppercase, exactly versioned RefSeq transcript accession."""
    normalized = _clean_identifier_text(accession).upper()
    if not VERSIONED_REFSEQ_TRANSCRIPT_RE.fullmatch(normalized):
        if VERSIONED_REFSEQ_GENOMIC_RE.fullmatch(normalized):
            raise ValueError(
                f"{normalized} is a genomic RefSeq accession, not a transcript "
                "accession. Transcript accession mode accepts exact NM/XM/NR/XR "
                "versions. To compare locally, paste one transcript sequence or "
                "choose a one-record FASTA/text file. Whole-chromosome NC records "
                "are not supported in transcript mode."
            )
        raise ValueError(
            "Private panel accessions must include an exact RefSeq transcript "
            f"version such as NM_000041.4; received: {accession}"
        )
    return normalized


def extract_refseq_accession_from_header(header: str) -> str:
    """Extract a versioned RefSeq transcript accession from a FASTA header."""
    cleaned_header = _clean_identifier_text(header)
    for token in re.split(r"[|\s]+", cleaned_header):
        if VERSIONED_REFSEQ_TRANSCRIPT_RE.fullmatch(token):
            return token.upper()
    return cleaned_header.split(maxsplit=1)[0].strip("|").upper()
