from __future__ import annotations

import pytest

from tools_for_pharma.oligo.transcript_accessions import (
    extract_refseq_accession_from_header,
    normalize_versioned_refseq_accession,
)


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        (" nm_000041.4 ", "NM_000041.4"),
        ("xm_017001.2", "XM_017001.2"),
        ("NR_001234.1", "NR_001234.1"),
        ("xr_123456.7", "XR_123456.7"),
    ],
)
def test_normalize_versioned_refseq_transcript_accessions(
    raw: str,
    expected: str,
) -> None:
    assert normalize_versioned_refseq_accession(raw) == expected


@pytest.mark.parametrize("prefix", ["NC", "NG", "NT", "NW"])
def test_genomic_refseq_accessions_keep_distinct_error(prefix: str) -> None:
    with pytest.raises(ValueError) as error_info:
        normalize_versioned_refseq_accession(f"{prefix}_000005.10")

    message = str(error_info.value)
    assert "genomic RefSeq accession" in message
    assert "not a transcript accession" in message
    assert "paste one transcript sequence" in message
    assert "one-record FASTA/text file" in message


def test_unversioned_or_non_refseq_accession_uses_transcript_error() -> None:
    with pytest.raises(ValueError, match="must include an exact RefSeq transcript"):
        normalize_versioned_refseq_accession("NM_000041")


def test_header_parser_finds_refseq_token_and_preserves_fallback() -> None:
    assert (
        extract_refseq_accession_from_header(
            "gi|123|ref|nm_000041.4| transcript description"
        )
        == "NM_000041.4"
    )
    assert extract_refseq_accession_from_header("custom_id description") == "CUSTOM_ID"
