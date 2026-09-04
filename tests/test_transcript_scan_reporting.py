"""Contracts for transcript-scan text and workbook reporting."""

from __future__ import annotations

import ast
from pathlib import Path
from types import SimpleNamespace

import pandas as pd

from tools_for_pharma.oligo import ncbi_blast
from tools_for_pharma.oligo.transcript_scan import reporting
from tools_for_pharma.oligo.transcript_scan.models import (
    AntisenseQuery,
    AntisenseRegion,
    ComparisonResult,
    PrivatePanelScanResult,
    QueryTargetSummary,
    TranscriptMatch,
    TranscriptTargetResult,
)


REPORTING_EXPORTS = (
    "comparison_result_rows",
    "default_gui_result_workbook",
    "default_private_panel_workbook",
    "default_result_workbook",
    "format_closest_transcript_matches_for_terminal",
    "format_single_sequence_scan_result",
    "format_transcript_matches_for_terminal",
    "input_query_rows",
    "metadata_rows",
    "transcript_match_rows",
    "transcript_matches_to_csv",
    "write_excel_workbook",
    "write_result_workbook",
)


def _match() -> TranscriptMatch:
    return TranscriptMatch(
        transcript_name="NM_000001.1 target",
        antisense_name="AS_demo",
        scan_region="full",
        as_region_start=1,
        as_region_end=21,
        antisense_5to3="A" * 21,
        antisense_region_5to3="A" * 21,
        target_5to3="U" * 21,
        transcript_start=10,
        transcript_end=30,
        mismatches=2,
        transcript_window_5to3="A" + "U" * 19 + "A",
        transcript_match_as_5to3="U" + "A" * 19 + "U",
        mismatch_positions_1based=(1, 21),
        as_mismatch_positions_1based=(1, 21),
    )


def test_facade_reexports_reporting_functions() -> None:
    assert {
        name: getattr(ncbi_blast, name) is getattr(reporting, name)
        for name in REPORTING_EXPORTS
    } == {name: True for name in REPORTING_EXPORTS}


def test_generated_workbook_widths_keep_headers_and_short_values_visible(
    tmp_path: Path,
) -> None:
    from openpyxl import load_workbook

    workbook_path = tmp_path / "readable.xlsx"
    reporting.write_excel_workbook(
        workbook_path,
        {
            "comparison_results": [
                {
                    "query_name": "AS_demo",
                    "mismatch_positions_in_query_1based": "1;21",
                    "differences": "1:A>U; 21:A>U",
                }
            ]
        },
    )

    worksheet = load_workbook(workbook_path)["comparison_results"]

    assert worksheet.column_dimensions["A"].width >= len("query_name") + 2
    assert worksheet.column_dimensions["B"].width >= len(
        "mismatch_positions_in_query_1based"
    ) + 2
    assert worksheet.column_dimensions["C"].width >= len("1:A>U; 21:A>U") + 2


def test_input_query_rows_preserve_noncanonical_source_fields() -> None:
    query = AntisenseQuery(
        name="AS_demo",
        sequence_5to3="AUGC",
        target_accession="NM_000001.1",
        source_fields={"Pos20": 20, "antisense_name": "do not overwrite"},
    )

    row = reporting.input_query_rows([query])[0]

    assert row["Pos20"] == 20
    assert row["antisense_name"] == "AS_demo"
    assert row["target_accession"] == "NM_000001.1"


def test_reporting_rows_keep_one_based_mismatch_coordinates() -> None:
    match_row = reporting.transcript_match_rows([_match()])[0]
    comparison = ComparisonResult(
        input_order=1,
        query_name="AS_demo",
        target_accession="NM_000001.1",
        scan_region="full",
        region_start=1,
        region_end=21,
        result="mismatch",
        sites_within_threshold=1,
        best_mismatches=2,
        mismatch_positions_in_query_1based=(1, 21),
        best_transcript_start=10,
        best_transcript_end=30,
        query_region_5to3="A" * 21,
        best_match_in_query_orientation_5to3="U" + "A" * 19 + "U",
        differences="1:A>U; 21:A>U",
    )
    comparison_row = reporting.comparison_result_rows([comparison])[0]

    assert match_row["mismatch_positions_1based"] == "1;21"
    assert match_row["as_mismatch_positions_1based"] == "1;21"
    assert comparison_row["region_start"] == 1
    assert comparison_row["region_end"] == 21
    assert comparison_row["mismatch_positions_in_query_1based"] == "1;21"


def test_comparison_workbook_schema_and_sheet_order_are_stable(tmp_path) -> None:
    query = AntisenseQuery(
        name="AS_demo",
        sequence_5to3="A" * 21,
        target_accession="NM_000001.1",
        source_fields={"Pos20": 20},
    )
    comparison = ComparisonResult(
        input_order=1,
        query_name=query.name,
        target_accession=query.target_accession,
        scan_region="full",
        region_start=1,
        region_end=21,
        result="mismatch",
        sites_within_threshold=1,
        best_mismatches=2,
        mismatch_positions_in_query_1based=(1, 21),
        best_transcript_start=10,
        best_transcript_end=30,
        query_region_5to3=query.sequence_5to3,
        best_match_in_query_orientation_5to3="U" + "A" * 19 + "U",
        differences="1:A>U; 21:A>U",
    )
    target = TranscriptTargetResult(
        requested_accession=query.target_accession,
        retrieved_accession=query.target_accession,
        status="ready",
    )
    args = ncbi_blast.build_parser().parse_args(
        ["--as-sequence", query.sequence_5to3, "--target-sequence", "U" * 21]
    )
    workbook_path = tmp_path / "comparison.xlsx"

    reporting.write_result_workbook(
        workbook_path,
        args,
        [query],
        [AntisenseRegion("full")],
        [_match()],
        [],
        "start",
        "end",
        include_blast_sheets=False,
        comparison_results=[comparison],
        transcript_targets=[target],
    )

    assert pd.ExcelFile(workbook_path).sheet_names == [
        "input_queries",
        "comparison_results",
        "local_transcript_scan",
        "transcript_targets",
        "run_metadata",
    ]
    input_table = pd.read_excel(workbook_path, sheet_name="input_queries")
    comparison_table = pd.read_excel(
        workbook_path,
        sheet_name="comparison_results",
    )
    detail_table = pd.read_excel(
        workbook_path,
        sheet_name="local_transcript_scan",
    )
    assert input_table.loc[0, "Pos20"] == 20
    assert comparison_table.loc[0, "region_start"] == 1
    assert comparison_table.loc[0, "region_end"] == 21
    assert comparison_table.loc[0, "mismatch_positions_in_query_1based"] == "1;21"
    assert "transcript_match_as_5to3" not in detail_table.columns
    assert list(detail_table.columns) == list(
        reporting.transcript_match_rows(
            [_match()],
            include_as_oriented_match=False,
        )[0]
    )


def test_single_sequence_text_output_contract() -> None:
    query = AntisenseQuery("AS_demo", "AUGC")
    region = AntisenseRegion("full")
    match = TranscriptMatch(
        transcript_name="manual transcript",
        antisense_name=query.name,
        scan_region=region.name,
        as_region_start=1,
        as_region_end=4,
        antisense_5to3=query.sequence_5to3,
        antisense_region_5to3=query.sequence_5to3,
        target_5to3="GCAU",
        transcript_start=4,
        transcript_end=7,
        mismatches=0,
        transcript_window_5to3="GCAU",
        transcript_match_as_5to3="AUGC",
        mismatch_positions_1based=(),
        as_mismatch_positions_1based=(),
    )
    target = TranscriptTargetResult(
        requested_accession="manual transcript",
        transcript_name="manual transcript",
        sequence_length_nt=10,
        cache_status="pasted sequence",
        status="ready",
    )
    summary = QueryTargetSummary(
        query_name=query.name,
        sequence_type="AS",
        requested_accession="manual transcript",
        retrieved_accession="",
        target_status="ready",
        scan_status="matched",
        scan_regions="full",
        match_count=1,
        exact_match_count=1,
        best_mismatches=0,
    )
    result = PrivatePanelScanResult(
        targets=(target,),
        matches=(match,),
        summaries=(summary,),
        closest_matches=(match,),
    )

    text = reporting.format_single_sequence_scan_result(
        SimpleNamespace(max_mismatches=0),
        [query],
        [region],
        result,
    )

    assert text.startswith("LOCAL SINGLE-SEQUENCE TRANSCRIPT SCAN\n")
    assert "Transcript source: pasted sequence\n" in text
    assert "NCBI transcript retrieval: Not used\n" in text
    assert "Guide sequence sent to NCBI: No\n" in text
    assert "Matches within threshold: 1\n" in text
    assert "Region: full\n" in text
    assert text.rstrip().splitlines()[-1].split() == [
        "AS_demo",
        "full",
        "4",
        "7",
        "0",
        "GCAU",
        "AUGC",
        "-",
        "-",
    ]


def test_reporting_module_does_not_import_gui_or_workflow_facade() -> None:
    source_path = Path(reporting.__file__)
    tree = ast.parse(source_path.read_text(encoding="utf-8"))
    imports = {
        alias.name
        for node in ast.walk(tree)
        if isinstance(node, ast.Import)
        for alias in node.names
    }
    imports.update(
        node.module or ""
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom)
    )

    assert "tkinter" not in imports
    assert "tools_for_pharma.oligo.ncbi_blast" not in imports
