"""Tests for reusable nucleotide sequence primitives."""

from __future__ import annotations

import pytest

from tools_for_pharma.sequence.nucleotides import (
    SequenceNormalizationError,
    complement_dna,
    complement_rna,
    normalize_dna,
    normalize_nucleotides,
    normalize_rna,
    reverse_complement_dna,
    reverse_complement_rna,
    subsequence_1based,
)
from tools_for_pharma.oligo.transcript_panel import normalize_dna as normalize_panel_dna


def test_normalize_rna_and_dna_convert_alphabets() -> None:
    assert normalize_rna(" a-c g.t 123\n") == "ACGU"
    assert normalize_dna(" a-c g.u 123\n") == "ACGT"


def test_normalization_requires_explicit_ambiguity_codes() -> None:
    with pytest.raises(SequenceNormalizationError) as error_info:
        normalize_dna("ACGN")

    assert error_info.value.reason == "invalid"
    assert error_info.value.invalid_bases == ("N",)
    assert normalize_dna("ACGN", allowed_ambiguity_codes="N") == "ACGN"


def test_whitespace_cleanup_preserves_other_invalid_characters() -> None:
    with pytest.raises(SequenceNormalizationError) as error_info:
        normalize_nucleotides(
            "AC-G N",
            alphabet="DNA",
            allowed_ambiguity_codes="N",
            cleanup="whitespace",
        )

    assert error_info.value.invalid_bases == ("-",)


def test_transcript_panel_adapter_preserves_existing_validation_messages() -> None:
    assert normalize_panel_dna("A C G U N") == "ACGTN"

    with pytest.raises(ValueError, match=r"^Sequence is empty\.$"):
        normalize_panel_dna(" \n ")
    with pytest.raises(ValueError, match=r"^Unsupported sequence characters: -$"):
        normalize_panel_dna("AC-G")


def test_normalization_reports_empty_input_after_cleanup() -> None:
    with pytest.raises(SequenceNormalizationError) as error_info:
        normalize_rna(" - 123 ")

    assert error_info.value.reason == "empty"
    assert str(error_info.value) == "Sequence is empty after cleanup."


def test_explicit_complement_operations_preserve_direction_intent() -> None:
    assert complement_rna("AUGC") == "UACG"
    assert reverse_complement_rna("AUGC") == "GCAU"
    assert complement_dna("ATGC") == "TACG"
    assert reverse_complement_dna("ATGC") == "GCAT"
    assert reverse_complement_dna("ACGN", allowed_ambiguity_codes="N") == "NCGT"


def test_subsequence_1based_uses_inclusive_coordinates() -> None:
    assert subsequence_1based("AUGCU", start=2, end=4) == "UGC"

    with pytest.raises(ValueError, match="Start position must be 1 or greater"):
        subsequence_1based("AUGCU", start=0, end=4)
    with pytest.raises(ValueError, match="greater than or equal to start"):
        subsequence_1based("AUGCU", start=4, end=3)
    with pytest.raises(ValueError, match="at least 6 nt"):
        subsequence_1based("AUGCU", start=2, end=6)
