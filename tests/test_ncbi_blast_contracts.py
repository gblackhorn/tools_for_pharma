"""Compatibility contracts to preserve while refactoring transcript scanning."""

from __future__ import annotations

from pathlib import Path

import pytest

import transcript_scan_app
from tools_for_pharma.oligo import ncbi_blast
from tools_for_pharma.oligo.core import (
    get_complementary_sequence,
    get_subsequence,
    normalize_rna,
)
from tools_for_pharma.oligo.transcript_panel import (
    normalize_dna as normalize_panel_dna,
)


REQUIRED_NCBI_BLAST_EXPORTS = {
    "AntisenseQuery",
    "AntisenseRegion",
    "BlastSubmission",
    "NcbiBlastClient",
    "NcbiHttpClient",
    "TranscriptTargetResult",
    "application_data_dir",
    "assign_unique_blast_query_ids",
    "batch_antisense_queries",
    "build_parser",
    "closest_transcript_matches",
    "comparison_result_for_region",
    "comparison_result_rows",
    "default_gui_result_workbook",
    "default_result_workbook",
    "format_closest_transcript_matches_for_terminal",
    "format_single_sequence_scan_result",
    "format_transcript_matches_for_terminal",
    "gui_args",
    "gui_log_path",
    "input_query_rows",
    "multi_fasta",
    "panel_accessions_from_args",
    "parse_blast_field",
    "parse_plain_antisense_lines",
    "parse_scan_region",
    "read_antisense_queries",
    "read_target_accession_table",
    "read_transcript_input",
    "retrieve_transcript_targets",
    "run_blast_batches",
    "run_gui",
    "run_local_scan",
    "run_local_scan_with_comparison",
    "run_private_panel_scan",
    "run_private_panel_workflow",
    "run_single_sequence_scan",
    "scan_antisense_against_transcript",
    "scan_sense_against_transcript",
    "shared_gui_transcript_cache_dir",
    "single_sequence_gui_args",
    "single_sequence_gui_draft",
    "transcript_match_rows",
    "transcript_matches_to_csv",
    "validate_runtime_args",
    "write_result_workbook",
}


def test_ncbi_blast_keeps_repository_compatibility_exports() -> None:
    missing = sorted(
        name for name in REQUIRED_NCBI_BLAST_EXPORTS if not hasattr(ncbi_blast, name)
    )

    assert missing == []


def test_portable_entry_point_uses_ncbi_blast_compatibility_exports() -> None:
    assert transcript_scan_app.application_data_dir is ncbi_blast.application_data_dir
    assert transcript_scan_app.gui_log_path is ncbi_blast.gui_log_path
    assert transcript_scan_app.run_gui is ncbi_blast.run_gui
    assert (
        transcript_scan_app.shared_gui_transcript_cache_dir
        is ncbi_blast.shared_gui_transcript_cache_dir
    )


def test_application_base_dir_uses_repository_root_during_source_run(
    monkeypatch,
) -> None:
    monkeypatch.delattr(ncbi_blast.sys, "frozen", raising=False)

    assert ncbi_blast.application_base_dir() == Path(ncbi_blast.__file__).resolve().parents[2]


def test_application_base_dir_uses_executable_folder_when_frozen(
    tmp_path,
    monkeypatch,
) -> None:
    executable = tmp_path / "portable" / "TranscriptScan.exe"
    monkeypatch.setattr(ncbi_blast.sys, "frozen", True, raising=False)
    monkeypatch.setattr(ncbi_blast.sys, "executable", str(executable))

    assert ncbi_blast.application_base_dir() == executable.parent
    assert ncbi_blast.application_data_dir() == executable.parent / "TranscriptScanData"


def test_current_rna_normalization_and_coordinate_contract() -> None:
    assert normalize_rna(" a-c g.t 123\n") == "ACGU"
    assert get_subsequence("AUGCU", start=2, end=4) == "UGC"
    assert get_complementary_sequence("AUGC", reverse=False) == "UACG"
    assert get_complementary_sequence("AUGC", reverse=True) == "GCAU"

    with pytest.raises(ValueError, match="invalid bases: N"):
        normalize_rna("AUGN")


def test_current_dna_normalizers_have_distinct_ambiguity_contracts() -> None:
    assert ncbi_blast.normalize_dna("AUGC") == "ATGC"
    with pytest.raises(ValueError, match="invalid bases: N"):
        ncbi_blast.normalize_dna("ACGN")

    assert normalize_panel_dna("a c g u n") == "ACGTN"


def test_current_fasta_adapters_preserve_record_boundaries() -> None:
    queries = ncbi_blast.parse_fasta_records(
        ">AS_one description\nAUGC\n>AS_two\nCCGA\n"
    )

    assert [query.name for query in queries] == ["AS_one description", "AS_two"]
    assert [query.sequence_5to3 for query in queries] == ["AUGC", "CCGA"]

    with pytest.raises(ValueError, match="exactly one transcript record"):
        ncbi_blast.validate_single_transcript_record(
            ">first\nAUGC\n>second\nCCGA\n",
            "contract test",
        )


def test_query_fasta_adapter_preserves_names_and_blank_header_defaults() -> None:
    queries = ncbi_blast.parse_fasta_records(
        ">AS_one   description\nAUGC\n>\nCCGA\n"
    )

    assert [query.name for query in queries] == ["AS_one   description", "AS_2"]


def test_equal_length_mismatch_positions_are_one_based() -> None:
    assert ncbi_blast.mismatch_positions("AUGC", "AUAC") == (3,)

    with pytest.raises(ValueError, match="same length; received 4 and 3"):
        ncbi_blast.mismatch_positions("AUGC", "AUG")
