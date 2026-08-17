from __future__ import annotations

import pytest

from tools_for_pharma.oligo.transcript_panel import (
    FastaRecord,
    classify_reference_region,
    compare_transcript_to_reference,
    parse_fasta_records,
)


def test_parse_fasta_records_keeps_records_separate() -> None:
    records = parse_fasta_records(">ref description\nACGT\n>alt\nACGA\n")

    assert [record.identifier for record in records] == ["ref", "alt"]
    assert [record.sequence for record in records] == ["ACGT", "ACGA"]


def test_parse_fasta_records_rejects_sequence_before_header() -> None:
    with pytest.raises(ValueError, match="before the first FASTA header"):
        parse_fasta_records("ACGT")


def test_classify_reference_region() -> None:
    assert classify_reference_region(1, 4, 5, 10) == "5' UTR"
    assert classify_reference_region(5, 10, 5, 10) == "CDS"
    assert classify_reference_region(11, 12, 5, 10) == "3' UTR"
    assert classify_reference_region(4, 5, 5, 10) == "CDS/UTR boundary"


def test_compare_transcript_exact_match() -> None:
    reference = FastaRecord("NM_ref.1", "", "AACCGGTT")
    transcript = FastaRecord("ENST_alt.1", "", "AACCGGTT")

    summary, differences, conserved = compare_transcript_to_reference(
        reference,
        transcript,
        reference_cds_start_1based=3,
        reference_cds_end_1based=6,
    )

    assert summary.exact_full_length_match is True
    assert summary.conserved_bases_nt == 8
    assert summary.reference_exactly_conserved_pct == 1.0
    assert differences == []
    assert len(conserved) == 1


def test_compare_transcript_reports_substitution_and_deletion() -> None:
    reference = FastaRecord("NM_ref.1", "", "AAAACCCCGGGG")
    transcript = FastaRecord("ENST_alt.1", "", "AAAATCCCGG")

    summary, differences, _ = compare_transcript_to_reference(
        reference,
        transcript,
        reference_cds_start_1based=5,
        reference_cds_end_1based=10,
    )

    assert summary.exact_full_length_match is False
    assert summary.difference_block_count >= 2
    assert any(block.difference_type == "substitution" for block in differences)
    assert any(block.difference_type == "deletion from transcript" for block in differences)
    assert all(block.transcript_accession == "ENST_alt.1" for block in differences)


def test_compare_transcript_reports_insertion_boundary() -> None:
    reference = FastaRecord("NM_ref.1", "", "AAAACCCC")
    transcript = FastaRecord("ENST_alt.1", "", "AAAAGGCCCC")

    _, differences, _ = compare_transcript_to_reference(reference, transcript)

    insertion = next(
        block
        for block in differences
        if block.difference_type == "insertion in transcript"
    )
    assert insertion.reference_start_1based is None
    assert insertion.reference_boundary_after_1based == 4
