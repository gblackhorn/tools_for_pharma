from __future__ import annotations

import ast
from dataclasses import fields
from pathlib import Path

from tools_for_pharma.oligo import ncbi_blast
from tools_for_pharma.oligo.transcript_scan import models
from tools_for_pharma.oligo.transcript_scan import queries
from tools_for_pharma.oligo.transcript_scan import scanner


def test_legacy_facade_reexports_domain_models_and_operations() -> None:
    reexports = {
        "AntisenseQuery": models.AntisenseQuery,
        "AntisenseRegion": models.AntisenseRegion,
        "TranscriptMatch": models.TranscriptMatch,
        "TranscriptTargetResult": models.TranscriptTargetResult,
        "QueryTargetSummary": models.QueryTargetSummary,
        "ComparisonResult": models.ComparisonResult,
        "PrivatePanelScanResult": models.PrivatePanelScanResult,
        "clean_text_for_id": queries.clean_text_for_id,
        "sanitize_fasta_name": queries.sanitize_fasta_name,
        "assign_unique_blast_query_ids": queries.assign_unique_blast_query_ids,
        "normalize_sequence_type": queries.normalize_sequence_type,
        "default_query_name": queries.default_query_name,
        "batch_antisense_queries": queries.batch_antisense_queries,
        "parse_fasta_records": queries.parse_fasta_records,
        "parse_plain_antisense_lines": queries.parse_plain_antisense_lines,
        "read_antisense_file": queries.read_antisense_file,
        "duplicate_sequence_groups": queries.duplicate_sequence_groups,
        "parse_scan_region": queries.parse_scan_region,
        "parse_scan_regions": queries.parse_scan_regions,
        "antisense_region_sequence": queries.antisense_region_sequence,
        "mismatch_positions": scanner.mismatch_positions,
        "scan_antisense_against_transcript": (
            scanner.scan_antisense_against_transcript
        ),
        "scan_sense_against_transcript": scanner.scan_sense_against_transcript,
        "closest_transcript_matches": scanner.closest_transcript_matches,
        "comparison_result_for_region": scanner.comparison_result_for_region,
    }

    assert {
        name: getattr(ncbi_blast, name) is implementation
        for name, implementation in reexports.items()
    } == {name: True for name in reexports}


def test_query_and_match_dataclass_field_order_remains_stable() -> None:
    assert [field.name for field in fields(models.AntisenseQuery)] == [
        "name",
        "sequence_5to3",
        "target_accession",
        "target_gene",
        "species",
        "notes",
        "sequence_type",
        "blast_query_id",
        "source_fields",
    ]
    assert [field.name for field in fields(models.TranscriptMatch)] == [
        "transcript_name",
        "antisense_name",
        "scan_region",
        "as_region_start",
        "as_region_end",
        "antisense_5to3",
        "antisense_region_5to3",
        "target_5to3",
        "transcript_start",
        "transcript_end",
        "mismatches",
        "transcript_window_5to3",
        "transcript_match_as_5to3",
        "mismatch_positions_1based",
        "as_mismatch_positions_1based",
        "sequence_type",
    ]


def test_direct_domain_scanner_preserves_orientation_and_coordinates() -> None:
    as_matches = scanner.scan_antisense_against_transcript(
        "AUGC",
        "AAAGCAUCCC",
        transcript_name="NM_test.1",
        antisense_name="AS_test",
        max_mismatches=0,
    )
    ss_matches = scanner.scan_sense_against_transcript(
        "AUGC",
        "CCCAUGCAAA",
        transcript_name="NM_test.1",
        sense_name="SS_test",
        max_mismatches=0,
    )

    assert [
        (
            match.sequence_type,
            match.transcript_start,
            match.transcript_end,
            match.transcript_window_5to3,
            match.transcript_match_as_5to3,
        )
        for match in [*as_matches, *ss_matches]
    ] == [
        ("AS", 4, 7, "GCAU", "AUGC"),
        ("SS", 4, 7, "AUGC", "AUGC"),
    ]


def test_direct_best_match_reports_full_query_mismatch_positions() -> None:
    query = models.AntisenseQuery("AS_seed", "A" * 21)
    region = queries.parse_scan_region("seed:2-8")
    transcript = "UUUUUUG"
    all_matches = scanner.scan_antisense_against_transcript(
        query.sequence_5to3,
        transcript,
        scan_region=region,
        max_mismatches=None,
    )

    result = scanner.comparison_result_for_region(
        input_order=1,
        query=query,
        target_accession="NM_test.1",
        scan_region=region,
        all_matches=all_matches,
    )

    assert result.region_start == 2
    assert result.region_end == 8
    assert result.best_mismatches == 1
    assert result.mismatch_positions_in_query_1based == (2,)
    assert result.differences == "2:A>C"


def test_domain_modules_do_not_import_interface_or_transport_dependencies() -> None:
    forbidden_roots = {
        "argparse",
        "openpyxl",
        "pandas",
        "tkinter",
        "tools_for_pharma.oligo.ncbi_transport",
        "tools_for_pharma.shared.excel_utils",
    }
    module_paths = [
        Path(models.__file__),
        Path(queries.__file__),
        Path(scanner.__file__),
    ]

    imported_names: set[str] = set()
    for module_path in module_paths:
        tree = ast.parse(module_path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported_names.update(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imported_names.add(node.module)

    assert imported_names.isdisjoint(forbidden_roots)
