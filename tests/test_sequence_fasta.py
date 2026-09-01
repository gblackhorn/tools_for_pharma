"""Tests for reusable FASTA structure primitives."""

from __future__ import annotations

import pytest

from tools_for_pharma.oligo.transcript import fasta_or_plain_text_to_sequence
from tools_for_pharma.sequence.fasta import (
    FastaRecord,
    format_fasta,
    parse_fasta,
    require_single_fasta_record,
)


def test_parse_fasta_preserves_records_and_leaves_alphabet_to_caller() -> None:
    records = parse_fasta(
        "; file comment\n>ref description\nacgu\n; record comment\n>alt\nACGN\n"
    )

    assert records == [
        FastaRecord("ref", "description", "acgu"),
        FastaRecord("alt", "", "ACGN"),
    ]
    assert records[0].header == "ref description"


def test_parse_fasta_rejects_invalid_structure() -> None:
    with pytest.raises(ValueError, match="before the first FASTA header"):
        parse_fasta("ACGT\n>record\nACGT\n")
    with pytest.raises(ValueError, match="header cannot be blank"):
        parse_fasta(">\nACGT\n")
    with pytest.raises(ValueError, match="No FASTA records found"):
        parse_fasta("\n; comment only\n")


def test_require_single_fasta_record_reports_cardinality() -> None:
    one = FastaRecord("one", "", "ACGT")
    assert require_single_fasta_record([one]) is one

    with pytest.raises(ValueError, match="exactly one FASTA record; found 2"):
        require_single_fasta_record(
            [one, FastaRecord("two", "", "TGCA")],
            source_label="Transcript input",
        )


def test_format_fasta_controls_width_and_trailing_newline() -> None:
    record = FastaRecord("ref", "description", "AACCGGTT")

    assert format_fasta(record, width=4) == ">ref description\nAACC\nGGTT\n"
    assert format_fasta(record, width=4, trailing_newline=False) == (
        ">ref description\nAACC\nGGTT"
    )

    with pytest.raises(ValueError, match="positive integer"):
        format_fasta(record, width=0)


def test_single_transcript_adapter_accepts_plain_or_one_fasta_record() -> None:
    assert fasta_or_plain_text_to_sequence("ACGT\n; comment\n") == "ACGU"
    assert fasta_or_plain_text_to_sequence(">target\nACGT\n") == "ACGU"

    with pytest.raises(ValueError, match="exactly one FASTA record; found 2"):
        fasta_or_plain_text_to_sequence(">first\nACGT\n>second\nTGCA\n")
