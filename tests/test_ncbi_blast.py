from __future__ import annotations

import pandas as pd
import pytest
from urllib.error import URLError

import tools_for_pharma.oligo.ncbi_blast as ncbi_blast
from tools_for_pharma.oligo.transcript_scan import app_services

from tools_for_pharma.oligo.ncbi_blast import (
    AntisenseRegion,
    AntisenseQuery,
    BlastSubmission,
    NcbiBlastClient,
    NcbiHttpClient,
    TranscriptTargetResult,
    assign_unique_blast_query_ids,
    batch_antisense_queries,
    build_parser,
    closest_transcript_matches,
    comparison_result_for_region,
    comparison_result_rows,
    default_gui_result_workbook,
    default_result_workbook,
    format_closest_transcript_matches_for_terminal,
    format_single_sequence_scan_result,
    format_transcript_matches_for_terminal,
    gui_args,
    input_query_rows,
    multi_fasta,
    parse_blast_field,
    panel_accessions_from_args,
    parse_plain_antisense_lines,
    parse_scan_region,
    read_antisense_queries,
    read_target_accession_table,
    read_transcript_input,
    run_blast_batches,
    run_local_scan,
    run_local_scan_with_comparison,
    run_private_panel_scan,
    run_private_panel_workflow,
    run_single_sequence_scan,
    retrieve_transcript_targets,
    scan_antisense_against_transcript,
    scan_sense_against_transcript,
    shared_gui_transcript_cache_dir,
    single_sequence_gui_draft,
    single_sequence_gui_args,
    transcript_match_rows,
    transcript_matches_to_csv,
    validate_runtime_args,
    write_result_workbook,
)


def test_parse_blast_field_reads_rid_and_rtoe() -> None:
    text = "    RID = ABC123\n    RTOE = 42\n"

    assert parse_blast_field(text, "RID") == "ABC123"
    assert parse_blast_field(text, "RTOE") == "42"


def test_ncbi_contact_email_has_no_default_and_requires_valid_input() -> None:
    with pytest.raises(ValueError, match="valid contact email"):
        ncbi_blast.require_email("")
    with pytest.raises(ValueError, match="valid contact email"):
        ncbi_blast.require_email("not-an-email")

    assert ncbi_blast.require_email(" user@example.com ") == "user@example.com"


def test_nc_genomic_accession_error_explains_local_target_alternatives() -> None:
    with pytest.raises(ValueError) as error_info:
        ncbi_blast.normalize_versioned_refseq_accession("NC_000005.10")

    message = str(error_info.value)
    assert "genomic RefSeq accession" in message
    assert "not a transcript accession" in message
    assert "paste one transcript sequence" in message
    assert "one-record FASTA/text file" in message


def test_portable_gui_settings_and_cache_use_local_data_subfolder(
    tmp_path,
    monkeypatch,
) -> None:
    data_dir = tmp_path / "TranscriptScanData"
    monkeypatch.setattr(app_services, "application_data_dir", lambda: data_dir)

    saved_path = app_services.save_gui_settings({"ncbi_email": "user@example.com"})

    assert saved_path == data_dir / "settings.json"
    assert app_services.load_gui_settings() == {"ncbi_email": "user@example.com"}
    assert app_services.shared_gui_transcript_cache_dir() == data_dir / "transcript_cache"
    assert app_services.gui_log_path() == data_dir / "logs" / "transcript_scan.log"


def test_scan_antisense_against_transcript_finds_reverse_complement() -> None:
    matches = scan_antisense_against_transcript(
        antisense_5to3="AUGC",
        transcript_sequence="GGGCAUTTT",
        transcript_name="test_transcript",
        max_mismatches=0,
    )

    assert len(matches) == 1
    assert matches[0].target_5to3 == "GCAU"
    assert matches[0].transcript_start == 3
    assert matches[0].transcript_end == 6
    assert matches[0].mismatches == 0
    assert matches[0].transcript_window_5to3 == "GCAU"
    assert matches[0].transcript_match_as_5to3 == "AUGC"


def test_scan_sense_against_transcript_finds_direct_sequence() -> None:
    matches = scan_sense_against_transcript(
        sense_5to3="AUGC",
        transcript_sequence="GGGAUGCUUU",
        transcript_name="test_transcript",
        max_mismatches=0,
    )

    assert len(matches) == 1
    assert matches[0].sequence_type == "SS"
    assert matches[0].target_5to3 == "AUGC"
    assert matches[0].transcript_start == 4
    assert matches[0].transcript_end == 7
    assert matches[0].mismatches == 0
    assert matches[0].transcript_window_5to3 == "AUGC"
    assert matches[0].transcript_match_as_5to3 == "AUGC"


def test_scan_antisense_against_transcript_reports_mismatches() -> None:
    matches = scan_antisense_against_transcript(
        antisense_5to3="AUGC",
        transcript_sequence="GGGCAAUUU",
        transcript_name="test_transcript",
        max_mismatches=1,
    )

    assert len(matches) == 1
    assert matches[0].transcript_window_5to3 == "GCAA"
    assert matches[0].transcript_match_as_5to3 == "UUGC"
    assert matches[0].mismatches == 1
    assert matches[0].mismatch_positions_1based == (4,)
    assert matches[0].as_mismatch_positions_1based == (1,)


def test_transcript_matches_to_csv_includes_expected_columns() -> None:
    matches = scan_antisense_against_transcript(
        antisense_5to3="AUGC",
        transcript_sequence="GGGCAUTTT",
        transcript_name="test_transcript",
        max_mismatches=0,
    )

    text = transcript_matches_to_csv(matches)

    assert "transcript_name,antisense_name,scan_region,as_region_start" in text
    assert "transcript_window_5to3,transcript_match_as_5to3" in text
    assert "test_transcript,antisense_query,full,1,4,AUGC,AUGC,GCAU,3,6,0,GCAU,AUGC,," in text


def test_workbook_rows_can_omit_as_oriented_match_without_changing_default_schema() -> None:
    matches = scan_antisense_against_transcript(
        antisense_5to3="AUGC",
        transcript_sequence="GGGCAUTTT",
        max_mismatches=0,
    )

    default_rows = transcript_match_rows(matches)
    workbook_rows = transcript_match_rows(matches, include_as_oriented_match=False)

    assert default_rows[0]["transcript_match_as_5to3"] == "AUGC"
    assert "transcript_match_as_5to3" not in workbook_rows[0]


def test_gui_workbook_omits_blast_sheets_but_default_workbook_keeps_them(tmp_path) -> None:
    args = build_parser().parse_args(
        ["--as-sequence", "AUGC", "--target-sequence", "GCAU"]
    )
    queries = [AntisenseQuery("AS_demo", "AUGC")]
    regions = [AntisenseRegion("full")]
    matches = scan_antisense_against_transcript(
        antisense_5to3="AUGC",
        transcript_sequence="GCAU",
        max_mismatches=0,
    )
    gui_path = tmp_path / "gui.xlsx"
    full_path = tmp_path / "full.xlsx"

    write_result_workbook(
        gui_path,
        args,
        queries,
        regions,
        matches,
        [],
        "start",
        "end",
        include_blast_sheets=False,
    )
    write_result_workbook(
        full_path,
        args,
        queries,
        regions,
        matches,
        [],
        "start",
        "end",
    )

    assert pd.ExcelFile(gui_path).sheet_names == [
        "input_queries",
        "local_transcript_scan",
        "run_metadata",
    ]
    assert pd.ExcelFile(full_path).sheet_names == [
        "input_queries",
        "local_transcript_scan",
        "blast_hits_raw",
        "blast_hits_filtered",
        "blast_batches",
        "run_metadata",
    ]


def test_gui_comparison_results_are_second_and_detail_sheet_is_unchanged(tmp_path) -> None:
    args = build_parser().parse_args(
        ["--as-sequence", "AUGC", "--target-sequence", "GCAU"]
    )
    query = AntisenseQuery("AS_demo", "AUGC")
    region = AntisenseRegion("full")
    matches = scan_antisense_against_transcript(
        query.sequence_5to3,
        "GCAU",
        transcript_name="NM_000001.1 target",
        antisense_name=query.name,
        scan_region=region,
        max_mismatches=None,
    )
    comparison = comparison_result_for_region(
        input_order=1,
        query=query,
        target_accession="NM_000001.1",
        scan_region=region,
        passing_matches=matches,
        all_matches=matches,
    )
    target = TranscriptTargetResult(
        requested_accession="NM_000001.1",
        retrieved_accession="NM_000001.1",
        status="ready",
    )
    workbook = tmp_path / "gui_comparison.xlsx"

    write_result_workbook(
        workbook,
        args,
        [query],
        [region],
        matches,
        [],
        "start",
        "end",
        include_blast_sheets=False,
        comparison_results=[comparison],
        transcript_targets=[target],
    )

    assert pd.ExcelFile(workbook).sheet_names == [
        "input_queries",
        "comparison_results",
        "local_transcript_scan",
        "transcript_targets",
        "run_metadata",
    ]
    expected_detail_columns = list(
        transcript_match_rows(matches, include_as_oriented_match=False)[0]
    )
    assert list(pd.read_excel(workbook, sheet_name="local_transcript_scan").columns) == (
        expected_detail_columns
    )
    assert list(pd.read_excel(workbook, sheet_name="comparison_results").columns) == [
        "input_order",
        "query_name",
        "target_accession",
        "scan_region",
        "region_start",
        "region_end",
        "result",
        "sites_within_threshold",
        "best_mismatches",
        "mismatch_positions_in_query_1based",
        "best_transcript_start",
        "best_transcript_end",
        "query_region_5to3",
        "best_match_in_query_orientation_5to3",
        "differences",
    ]


def test_format_transcript_matches_for_terminal_shows_quick_summary() -> None:
    query = AntisenseQuery("AS_demo", "AUGC")
    region = AntisenseRegion("full")
    matches = scan_antisense_against_transcript(
        antisense_5to3=query.sequence_5to3,
        transcript_sequence="GGGCAUTTT",
        transcript_name="test_transcript",
        antisense_name=query.name,
        scan_region=region,
        max_mismatches=0,
    )

    text = format_transcript_matches_for_terminal(matches, [query], [region], max_mismatches=0)

    assert "Local transcript scan" in text
    assert "AS name: AS_demo" in text
    assert "Transcript: test_transcript" in text
    assert "Matches: 1" in text
    assert "start" in text
    assert "matched_as_5to3" in text
    assert "3" in text
    assert "AUGC" in text


def test_closest_transcript_matches_can_ignore_mismatch_cutoff() -> None:
    matches = scan_antisense_against_transcript(
        antisense_5to3="AUGC",
        transcript_sequence="GGGAAAUGGCAUTTT",
        transcript_name="test_transcript",
        max_mismatches=None,
    )

    closest = closest_transcript_matches(matches, 2)
    text = format_closest_transcript_matches_for_terminal(
        closest,
        closest=2,
        max_mismatches=0,
    )

    assert len(closest) == 2
    assert closest[0].mismatches == 0
    assert closest[1].mismatches > 0
    assert "Closest transcript windows" in text
    assert "not filtered by --max-mismatches 0" in text


def test_scan_antisense_against_transcript_supports_subregions() -> None:
    matches = scan_antisense_against_transcript(
        antisense_5to3="AUGC",
        transcript_sequence="GGCAUU",
        transcript_name="test_transcript",
        scan_region=AntisenseRegion("seed", 2, 3),
        max_mismatches=0,
    )

    assert len(matches) == 1
    assert matches[0].scan_region == "seed"
    assert matches[0].as_region_start == 2
    assert matches[0].as_region_end == 3
    assert matches[0].antisense_region_5to3 == "UG"
    assert matches[0].target_5to3 == "CA"
    assert matches[0].transcript_start == 3
    assert matches[0].transcript_match_as_5to3 == "UG"


def test_comparison_result_reports_full_query_positions_in_query_orientation() -> None:
    query = AntisenseQuery("AS_21", "A" * 21)
    matched_in_query_orientation = "C" + "A" * 19 + "G"
    transcript = ncbi_blast.get_complementary_sequence(
        matched_in_query_orientation,
        reverse=True,
    )
    all_matches = scan_antisense_against_transcript(
        query.sequence_5to3,
        transcript,
        transcript_name="NM_000001.1 target",
        scan_region=AntisenseRegion("full"),
        max_mismatches=None,
    )

    result = comparison_result_for_region(
        input_order=1,
        query=query,
        target_accession="NM_000001.1",
        scan_region=AntisenseRegion("full"),
        passing_matches=all_matches,
        all_matches=all_matches,
    )
    row = comparison_result_rows([result])[0]

    assert row["region_start"] == 1
    assert row["region_end"] == 21
    assert row["mismatch_positions_in_query_1based"] == "1;21"
    assert row["query_region_5to3"] == "A" * 21
    assert row["best_match_in_query_orientation_5to3"] == matched_in_query_orientation
    assert row["differences"] == "1:A>C; 21:A>G"


def test_comparison_result_offsets_seed_mismatches_to_full_query_coordinates() -> None:
    query = AntisenseQuery("AS_seed", "A" * 21)
    region = AntisenseRegion("seed", 2, 8)
    matched_seed = "C" + "A" * 6
    transcript = ncbi_blast.get_complementary_sequence(matched_seed, reverse=True)
    all_matches = scan_antisense_against_transcript(
        query.sequence_5to3,
        transcript,
        transcript_name="NM_000001.1 target",
        scan_region=region,
        max_mismatches=None,
    )

    result = comparison_result_for_region(
        input_order=1,
        query=query,
        target_accession="NM_000001.1",
        scan_region=region,
        passing_matches=all_matches,
        all_matches=all_matches,
    )
    row = comparison_result_rows([result])[0]

    assert row["region_start"] == 2
    assert row["region_end"] == 8
    assert row["mismatch_positions_in_query_1based"] == "2"
    assert row["differences"] == "2:A>C"


def test_parse_plain_antisense_lines_accepts_named_and_unnamed_sequences() -> None:
    records = parse_plain_antisense_lines("AS_001,AUGC\nCCGA\nAS_003\tUUAA\n")

    assert records == [
        AntisenseQuery("AS_001", "AUGC"),
        AntisenseQuery("AS_2", "CCGA"),
        AntisenseQuery("AS_003", "UUAA"),
    ]


def test_read_antisense_queries_from_table(tmp_path) -> None:
    table_path = tmp_path / "as_list.csv"
    table_path.write_text(
        "id,antisense,target_accession,species,notes\n"
        "AS_A,AUGC,NM_001,human,lead\n"
        "AS_B,CCGA,NM_002,mouse,backup\n",
        encoding="utf-8",
    )

    records = read_antisense_queries(as_table=table_path, as_name_column="id")

    assert records == [
        AntisenseQuery("AS_A", "AUGC", target_accession="NM_001", species="human", notes="lead"),
        AntisenseQuery("AS_B", "CCGA", target_accession="NM_002", species="mouse", notes="backup"),
    ]


def test_input_query_rows_include_additional_columns_from_source_excel(tmp_path) -> None:
    table_path = tmp_path / "as_list.xlsx"
    pd.DataFrame(
        {
            "id": ["AS_A", "AS_B"],
            "antisense": ["AUGC", "CCGA"],
            "Pos20": ["A", "G"],
            "project_code": [101, 102],
        }
    ).to_excel(table_path, index=False)

    records = read_antisense_queries(as_table=table_path, as_name_column="id")
    rows = input_query_rows(records)

    assert rows[0]["Pos20"] == "A"
    assert rows[1]["Pos20"] == "G"
    assert rows[0]["project_code"] == 101
    assert rows[0]["id"] == "AS_A"
    assert rows[0]["antisense"] == "AUGC"


def test_read_antisense_queries_accepts_ss_sequence() -> None:
    records = read_antisense_queries(ss_sequence="AUGC", ss_name="SS_A")

    assert records == [AntisenseQuery("SS_A", "AUGC", sequence_type="SS")]


def test_read_antisense_queries_from_ss_table(tmp_path) -> None:
    table_path = tmp_path / "ss_list.csv"
    table_path.write_text(
        "id,sense,target_accession\n"
        "SS_A,AUGC,NM_001\n",
        encoding="utf-8",
    )

    records = read_antisense_queries(
        ss_table=table_path,
        ss_name_column="id",
    )

    assert records == [
        AntisenseQuery("SS_A", "AUGC", target_accession="NM_001", sequence_type="SS"),
    ]


def test_batch_antisense_queries_groups_by_total_bases() -> None:
    records = [
        AntisenseQuery("AS_1", "AAAA"),
        AntisenseQuery("AS_2", "CCCC"),
        AntisenseQuery("AS_3", "GGGG"),
    ]

    batches = batch_antisense_queries(records, max_batch_bases=8)

    assert batches == [records[:2], records[2:]]


def test_parse_scan_region_accepts_named_range() -> None:
    assert parse_scan_region("seed:2-8") == AntisenseRegion("seed", 2, 8)


def test_input_query_rows_flags_duplicate_sequences() -> None:
    rows = input_query_rows(
        [
            AntisenseQuery("AS_A", "AUGC"),
            AntisenseQuery("AS_B", "CCGA"),
            AntisenseQuery("AS_C", "ATGC"),
        ]
    )

    assert rows[0]["is_duplicate_sequence"] is True
    assert rows[2]["is_duplicate_sequence"] is True
    assert rows[0]["duplicate_group_names"] == "AS_A;AS_C"
    assert rows[1]["is_duplicate_sequence"] is False


def test_blast_query_ids_remain_unique_after_sanitizing_collisions() -> None:
    records = [
        AntisenseQuery("A B", "AUGC"),
        AntisenseQuery("A_B", "CCGA"),
        AntisenseQuery("A/B", "UUAA"),
    ]

    prepared = assign_unique_blast_query_ids(records)

    assert [record.blast_query_id for record in prepared] == ["A_B", "A_B_2", "A_B_3"]
    assert [row["blast_query_id"] for row in input_query_rows(records)] == [
        "A_B",
        "A_B_2",
        "A_B_3",
    ]
    fasta_headers = [
        line for line in multi_fasta(records).splitlines() if line.startswith(">")
    ]
    assert fasta_headers == [
        ">A_B",
        ">A_B_2",
        ">A_B_3",
    ]


def test_local_scan_rejects_multi_record_target_fasta(tmp_path) -> None:
    target_path = tmp_path / "panel.fasta"
    target_path.write_text(">tx1\nAAAA\n>tx2\nCCCC\n", encoding="utf-8")

    with pytest.raises(ValueError, match="contains 2 FASTA records"):
        read_transcript_input(transcript_file=target_path)


def test_local_scan_reports_all_queries_missing_target_accessions() -> None:
    args = build_parser().parse_args(
        ["--as-sequence", "AUGC", "--target-accession-column", "target_accession"]
    )
    queries = [AntisenseQuery("AS_A", "AUGC"), AntisenseQuery("AS_B", "CCGA")]

    with pytest.raises(ValueError, match="AS_A, AS_B"):
        run_local_scan(args, queries, [AntisenseRegion("full")])


@pytest.mark.parametrize(
    ("option", "value", "message"),
    [
        ("--max-mismatches", "-1", "--max-mismatches"),
        ("--filter-max-mismatches", "-1", "--filter-max-mismatches"),
        ("--filter-max-gap-opens", "-1", "--filter-max-gap-opens"),
        ("--filter-min-alignment-fraction", "1.1", "--filter-min-alignment-fraction"),
        ("--hitlist-size", "0", "--hitlist-size"),
        ("--max-batch-bases", "0", "--max-batch-bases"),
        ("--timeout-seconds", "0", "--timeout-seconds"),
        ("--request-seconds", "-1", "--request-seconds"),
        ("--poll-seconds", "-1", "--poll-seconds"),
        ("--expect", "0", "--expect"),
    ],
)
def test_runtime_validation_rejects_invalid_numeric_options(
    option: str,
    value: str,
    message: str,
) -> None:
    args = build_parser().parse_args(
        ["--as-sequence", "AUGC", "--target-sequence", "GCAU", option, value]
    )

    with pytest.raises(ValueError, match=message):
        validate_runtime_args(args)


def test_megablast_uses_compatible_default_word_size() -> None:
    args = build_parser().parse_args(["--as-sequence", "AUGC", "--blast-only", "--megablast"])

    validate_runtime_args(args)

    assert args.word_size == 28


def test_megablast_rejects_other_incompatible_word_sizes() -> None:
    args = build_parser().parse_args(
        [
            "--as-sequence",
            "AUGC",
            "--blast-only",
            "--megablast",
            "--word-size",
            "11",
        ]
    )

    with pytest.raises(ValueError, match="Invalid --word-size 11 for megablast"):
        validate_runtime_args(args)


def test_remote_blast_warns_before_mocked_submission(capsys) -> None:
    class FakeBlastClient:
        def __init__(self, **_kwargs) -> None:
            pass

        def submit_blastn(self, **_kwargs) -> BlastSubmission:
            return BlastSubmission("RID_TEST", None)

        def wait_for_result(self, *_args, **_kwargs) -> None:
            pass

        def fetch_csv(self, *_args, **_kwargs) -> str:
            return ""

    args = build_parser().parse_args(
        [
            "--as-sequence",
            "AUGC",
            "--blast-only",
            "--email",
            "test@example.com",
        ]
    )
    validate_runtime_args(args)

    results = run_blast_batches(
        args,
        [AntisenseQuery("A B", "AUGC"), AntisenseQuery("A_B", "CCGA")],
        client_factory=FakeBlastClient,
    )

    captured = capsys.readouterr()
    assert "Privacy warning" in captured.err
    assert "transmit 2 oligo sequence(s)" in captured.err
    assert [query.blast_query_id for query in results[0].queries] == ["A_B", "A_B_2"]


def test_mocked_blast_client_submit_poll_fetch_lifecycle() -> None:
    class StubBlastClient(NcbiBlastClient):
        def __init__(self) -> None:
            super().__init__(email="test@example.com", request_seconds=0)
            self.requests: list[dict[str, object]] = []

        def get_text(self, _url: str, params: dict[str, object]) -> str:
            self.requests.append(params)
            if params.get("CMD") == "Put":
                return "RID = RID_TEST\n"
            if params.get("FORMAT_OBJECT") == "SearchInfo":
                return "Status = READY\n"
            return "query,subject,100,4,0,0,1,4,10,13,0.1,8\n"

    client = StubBlastClient()

    submission, csv_text = client.run_blastn(query_sequence="AUGC", timeout_seconds=1)

    assert submission.rid == "RID_TEST"
    assert csv_text.startswith("query,subject")
    assert [request["CMD"] for request in client.requests] == ["Put", "Get", "Get"]
    assert client.requests[0]["QUERY"] == ">oligo_query\nATGC"


def test_ncbi_http_failure_is_reported_without_real_network() -> None:
    def fail_urlopen(*_args, **_kwargs):
        raise URLError("offline")

    client = NcbiHttpClient(
        email="test@example.com",
        request_seconds=0,
        opener=fail_urlopen,
    )

    with pytest.raises(ValueError, match="NCBI request failed: offline"):
        client.get_text("https://example.invalid", {"id": "NM_TEST"})


def test_wait_for_blast_result_times_out_without_real_network() -> None:
    monotonic_values = iter([0.0, 2.0])
    client = NcbiBlastClient(
        email="test@example.com",
        request_seconds=0,
        monotonic=lambda: next(monotonic_values),
    )

    with pytest.raises(TimeoutError, match="Timed out waiting"):
        client.wait_for_result("RID_TEST", timeout_seconds=1)


def test_workbook_default_names_distinguish_local_and_remote_workflows(tmp_path) -> None:
    source = tmp_path / "queries.csv"
    local_args = build_parser().parse_args(
        ["--as-table", str(source), "--target-sequence", "GCAU"]
    )
    remote_args = build_parser().parse_args(
        ["--as-table", str(source), "--blast-only"]
    )

    assert default_result_workbook(local_args).name == "queries_ncbi_transcript_scan_results.xlsx"
    assert default_result_workbook(remote_args).name == "queries_ncbi_blast_results.xlsx"
    assert default_gui_result_workbook(source).name == "queries_ncbi_transcript_scan_results.xlsx"


def test_local_cli_default_workbook_omits_remote_blast_sheets(
    tmp_path,
    monkeypatch,
) -> None:
    source = tmp_path / "queries.csv"
    source.write_text("id,antisense\nAS_A,AUGC\n", encoding="utf-8")
    monkeypatch.setattr(
        ncbi_blast.sys,
        "argv",
        [
            "ncbi_blast",
            "--as-table",
            str(source),
            "--as-name-column",
            "id",
            "--target-sequence",
            "GCAU",
            "--max-mismatches",
            "0",
        ],
    )

    exit_code = ncbi_blast.main()

    workbook = tmp_path / "queries_ncbi_transcript_scan_results.xlsx"
    assert exit_code == 0
    assert pd.ExcelFile(workbook).sheet_names == [
        "input_queries",
        "local_transcript_scan",
        "run_metadata",
    ]


def test_panel_accessions_combine_repeated_flags_and_target_table(tmp_path) -> None:
    target_table = tmp_path / "targets.csv"
    target_table.write_text(
        "accession\nNM_000002.2\nNM_000001.1\n",
        encoding="utf-8",
    )
    args = build_parser().parse_args(
        [
            "--private-panel",
            "--target-accession",
            "NM_000001.1",
            "--target-table",
            str(target_table),
        ]
    )

    validate_runtime_args(args)

    assert panel_accessions_from_args(args) == ["NM_000001.1", "NM_000002.2"]
    assert read_target_accession_table(target_table) == [
        "NM_000002.2",
        "NM_000001.1",
    ]


def test_plain_target_list_can_include_a_header(tmp_path) -> None:
    target_list = tmp_path / "targets.txt"
    target_list.write_text(
        "accession\nNM_000001.1\nNM_000002.2\n",
        encoding="utf-8",
    )

    assert read_target_accession_table(target_list) == [
        "NM_000001.1",
        "NM_000002.2",
    ]


def test_missing_target_table_has_clear_validation_error(tmp_path) -> None:
    with pytest.raises(ValueError, match="does not exist"):
        read_target_accession_table(tmp_path / "missing.csv")


def test_private_panel_requires_exact_versioned_refseq_accessions() -> None:
    args = build_parser().parse_args(
        ["--private-panel", "--target-accession", "NM_000001"]
    )
    validate_runtime_args(args)

    with pytest.raises(ValueError, match="must include an exact RefSeq transcript version"):
        panel_accessions_from_args(args)


def test_single_target_accession_remains_backward_compatible_with_repeatable_parser(tmp_path) -> None:
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    (cache_dir / "NM_000001.1.fasta").write_text(
        ">NM_000001.1 cached transcript\nGCAU\n",
        encoding="utf-8",
    )
    args = build_parser().parse_args(
        [
            "--as-sequence",
            "AUGC",
            "--target-accession",
            "NM_000001.1",
            "--cache-dir",
            str(cache_dir),
            "--max-mismatches",
            "0",
        ]
    )
    validate_runtime_args(args)

    matches = run_local_scan(
        args,
        [AntisenseQuery("AS_A", "AUGC")],
        [AntisenseRegion("full")],
    )

    assert len(matches) == 1
    assert matches[0].transcript_name.startswith("NM_000001.1")


def test_gui_local_comparison_keeps_no_match_query_and_shows_closest_window() -> None:
    args = build_parser().parse_args(
        [
            "--as-sequence",
            "AAAA",
            "--target-sequence",
            "AAAA",
            "--max-mismatches",
            "0",
        ]
    )
    query = AntisenseQuery("AS_no_match", "AAAA")

    matches, comparisons = run_local_scan_with_comparison(
        args,
        [query],
        [AntisenseRegion("full")],
    )

    assert matches == []
    assert len(comparisons) == 1
    assert comparisons[0].result == "no_match"
    assert comparisons[0].sites_within_threshold == 0
    assert comparisons[0].best_mismatches == 4
    assert comparisons[0].best_transcript_start == 1
    assert comparisons[0].mismatch_positions_in_query_1based == (1, 2, 3, 4)


def test_private_panel_cannot_submit_guides_to_remote_blast() -> None:
    args = build_parser().parse_args(
        [
            "--as-sequence",
            "AUGC",
            "--private-panel",
            "--target-accession",
            "NM_000001.1",
            "--blast",
        ]
    )

    with pytest.raises(ValueError, match="cannot be combined with --blast"):
        validate_runtime_args(args)


def test_panel_retrieval_sends_accessions_only_and_writes_separate_cache_files(tmp_path) -> None:
    class RecordingClient:
        def __init__(self) -> None:
            self.requests: list[dict[str, object]] = []

        def get_text(self, _url: str, params: dict[str, object]) -> str:
            self.requests.append(params)
            accession = str(params["id"])
            sequence = "GCAU" if accession == "NM_000001.1" else "AUGC"
            return f">{accession} public transcript\n{sequence}\n"

    cache_dir = tmp_path / "cache"
    client = RecordingClient()

    targets = retrieve_transcript_targets(
        ["NM_000001.1", "NM_000002.2"],
        email="test@example.com",
        cache_dir=cache_dir,
        request_seconds=0,
        client=client,
    )

    assert [target.status for target in targets] == ["ready", "ready"]
    assert [target.cache_status for target in targets] == ["downloaded", "downloaded"]
    assert (cache_dir / "NM_000001.1.fasta").exists()
    assert (cache_dir / "NM_000002.2.fasta").exists()
    assert [request["id"] for request in client.requests] == [
        "NM_000001.1",
        "NM_000002.2",
    ]
    assert all("QUERY" not in request for request in client.requests)
    assert all("AUGC" not in request.values() for request in client.requests)


def test_panel_retrieval_reuses_valid_cache_without_network(tmp_path) -> None:
    class UnexpectedClient:
        def get_text(self, *_args, **_kwargs) -> str:
            raise AssertionError("A valid cached transcript must not be downloaded again.")

    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    (cache_dir / "NM_000001.1.fasta").write_text(
        ">NM_000001.1 cached transcript\nGCAU\n",
        encoding="utf-8",
    )

    targets = retrieve_transcript_targets(
        ["NM_000001.1"],
        email="test@example.com",
        cache_dir=cache_dir,
        request_seconds=0,
        client=UnexpectedClient(),
    )

    assert targets[0].status == "ready"
    assert targets[0].cache_status == "cache"


def test_panel_retrieval_refreshes_only_when_requested(tmp_path) -> None:
    class RefreshedClient:
        def get_text(self, *_args, **_kwargs) -> str:
            return ">NM_000001.1 refreshed transcript\nAAAA\n"

    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    cache_path = cache_dir / "NM_000001.1.fasta"
    cache_path.write_text(
        ">NM_000001.1 cached transcript\nGCAU\n",
        encoding="utf-8",
    )

    targets = retrieve_transcript_targets(
        ["NM_000001.1"],
        email="test@example.com",
        cache_dir=cache_dir,
        refresh=True,
        request_seconds=0,
        client=RefreshedClient(),
    )

    assert targets[0].status == "ready"
    assert targets[0].cache_status == "refreshed"
    assert "AAAA" in cache_path.read_text(encoding="utf-8")


def test_panel_retrieval_rejects_returned_version_mismatch_without_caching(tmp_path) -> None:
    class WrongVersionClient:
        def get_text(self, _url: str, _params: dict[str, object]) -> str:
            return ">NM_000001.2 wrong version\nGCAU\n"

    cache_dir = tmp_path / "cache"

    targets = retrieve_transcript_targets(
        ["NM_000001.1"],
        email="test@example.com",
        cache_dir=cache_dir,
        request_seconds=0,
        client=WrongVersionClient(),
    )

    assert targets[0].status == "error"
    assert "retrieved NM_000001.2" in targets[0].error
    assert not (cache_dir / "NM_000001.1.fasta").exists()


def test_panel_offline_mode_uses_cache_and_reports_missing_targets(tmp_path) -> None:
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    (cache_dir / "NM_000001.1.fasta").write_text(
        ">NM_000001.1 cached transcript\nGCAU\n",
        encoding="utf-8",
    )

    targets = retrieve_transcript_targets(
        ["NM_000001.1", "NM_000002.2"],
        email="test@example.com",
        cache_dir=cache_dir,
        offline=True,
    )

    assert targets[0].status == "ready"
    assert targets[0].cache_status == "cache"
    assert targets[1].status == "error"
    assert "Offline mode requires cached transcript" in targets[1].error


def test_panel_retrieval_cancellation_records_remaining_target_statuses(tmp_path) -> None:
    progress_events = []

    targets = retrieve_transcript_targets(
        ["NM_000001.1", "NM_000002.2"],
        email="test@example.com",
        cache_dir=tmp_path / "cache",
        offline=True,
        progress_callback=lambda *event: progress_events.append(event),
        cancel_check=lambda: True,
    )

    assert len(targets) == 2
    assert {target.error for target in targets} == {
        "Transcript retrieval cancelled by user."
    }
    assert progress_events[-1] == (2, 2, "NM_000002.2", "cancelled")


def test_private_panel_scan_builds_complete_as_ss_cross_product() -> None:
    queries = [
        AntisenseQuery("AS_A", "AUGC", sequence_type="AS"),
        AntisenseQuery("SS_A", "AUGC", sequence_type="SS"),
    ]
    targets = [
        TranscriptTargetResult(
            requested_accession="NM_000001.1",
            retrieved_accession="NM_000001.1",
            transcript_name="NM_000001.1 AS target",
            sequence_5to3="GCAU",
            sequence_length_nt=4,
            status="ready",
            exact_version_match=True,
        ),
        TranscriptTargetResult(
            requested_accession="NM_000002.2",
            retrieved_accession="NM_000002.2",
            transcript_name="NM_000002.2 SS target",
            sequence_5to3="AUGC",
            sequence_length_nt=4,
            status="ready",
            exact_version_match=True,
        ),
    ]

    result = run_private_panel_scan(
        queries,
        targets,
        [AntisenseRegion("full")],
        max_mismatches=0,
    )

    assert len(result.summaries) == 4
    summary_by_pair = {
        (summary.query_name, summary.requested_accession): summary
        for summary in result.summaries
    }
    assert summary_by_pair[("AS_A", "NM_000001.1")].scan_status == "matched"
    assert summary_by_pair[("AS_A", "NM_000002.2")].scan_status == "no_match"
    assert summary_by_pair[("SS_A", "NM_000001.1")].scan_status == "no_match"
    assert summary_by_pair[("SS_A", "NM_000002.2")].scan_status == "matched"
    assert len(result.matches) == 2


def test_private_panel_scan_records_target_errors_for_every_query() -> None:
    queries = [AntisenseQuery("AS_A", "AUGC"), AntisenseQuery("AS_B", "CCGA")]
    target = TranscriptTargetResult(
        requested_accession="NM_000001.1",
        status="error",
        error="cache missing",
    )

    result = run_private_panel_scan(
        queries,
        [target],
        [AntisenseRegion("full")],
        max_mismatches=0,
    )

    assert len(result.summaries) == 2
    assert {summary.scan_status for summary in result.summaries} == {"target_error"}
    assert {summary.error for summary in result.summaries} == {"cache missing"}


def test_private_panel_closest_windows_ignore_match_cutoff() -> None:
    query = AntisenseQuery("AS_A", "AUGC")
    target = TranscriptTargetResult(
        requested_accession="NM_000001.1",
        retrieved_accession="NM_000001.1",
        transcript_name="NM_000001.1 near target",
        sequence_5to3="GCAAUU",
        sequence_length_nt=6,
        status="ready",
        exact_version_match=True,
    )

    result = run_private_panel_scan(
        [query],
        [target],
        [AntisenseRegion("full")],
        max_mismatches=0,
        closest=2,
    )

    assert result.summaries[0].scan_status == "no_match"
    assert result.summaries[0].best_mismatches == 1
    assert len(result.closest_matches) == 2
    assert result.closest_matches[0].mismatches == 1


def test_private_panel_closest_windows_are_ranked_per_scan_region() -> None:
    query = AntisenseQuery("AS_A", "A" * 18)
    target = TranscriptTargetResult(
        requested_accession="NM_000001.1",
        retrieved_accession="NM_000001.1",
        transcript_name="NM_000001.1 target",
        sequence_5to3="A" * 25,
        sequence_length_nt=25,
        status="ready",
        exact_version_match=True,
    )

    result = run_private_panel_scan(
        [query],
        [target],
        [
            AntisenseRegion("full"),
            AntisenseRegion("seed", 2, 8),
            AntisenseRegion("core", 2, 18),
        ],
        max_mismatches=0,
        closest=2,
    )

    assert len(result.closest_matches) == 6
    assert {
        region: sum(match.scan_region == region for match in result.closest_matches)
        for region in ("full", "seed", "core")
    } == {"full": 2, "seed": 2, "core": 2}


def test_private_panel_cli_offline_workbook_contains_complete_status(tmp_path, monkeypatch) -> None:
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    (cache_dir / "NM_000001.1.fasta").write_text(
        ">NM_000001.1 AS target\nGCAU\n",
        encoding="utf-8",
    )
    (cache_dir / "NM_000002.2.fasta").write_text(
        ">NM_000002.2 no match\nAAAA\n",
        encoding="utf-8",
    )
    workbook = tmp_path / "panel.xlsx"
    monkeypatch.setattr(
        ncbi_blast.sys,
        "argv",
        [
            "ncbi_blast",
            "--as-sequence",
            "AUGC",
            "--private-panel",
            "--target-accession",
            "NM_000001.1",
            "--target-accession",
            "NM_000002.2",
            "--cache-dir",
            str(cache_dir),
            "--offline",
            "--max-mismatches",
            "0",
            "--result-workbook",
            str(workbook),
        ],
    )

    exit_code = ncbi_blast.main()

    assert exit_code == 0
    assert pd.ExcelFile(workbook).sheet_names == [
        "input_queries",
        "transcript_targets",
        "local_transcript_scan",
        "query_target_summary",
        "run_metadata",
    ]
    summaries = pd.read_excel(workbook, sheet_name="query_target_summary")
    assert summaries["scan_status"].tolist() == ["matched", "no_match"]
    metadata = pd.read_excel(workbook, sheet_name="run_metadata").set_index("key")["value"]
    assert metadata["privacy_mode"] == "local_guide_scan"
    assert metadata["guide_sequence_transmitted_to_ncbi"] in {False, "False", 0}


def test_gui_private_panel_puts_compact_results_second_and_omits_redundant_summary(
    tmp_path,
) -> None:
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    (cache_dir / "NM_000001.1.fasta").write_text(
        ">NM_000001.1 target\nGCAU\n",
        encoding="utf-8",
    )
    workbook = tmp_path / "gui_panel.xlsx"
    args = build_parser().parse_args(
        [
            "--as-sequence",
            "AUGC",
            "--private-panel",
            "--target-accession",
            "NM_000001.1",
            "--cache-dir",
            str(cache_dir),
            "--offline",
            "--max-mismatches",
            "0",
            "--result-workbook",
            str(workbook),
        ]
    )
    validate_runtime_args(args)

    exit_code = run_private_panel_workflow(
        args,
        "2026-08-12T00:00:00+00:00",
        include_comparison_results=True,
    )

    assert exit_code == 0
    assert pd.ExcelFile(workbook).sheet_names == [
        "input_queries",
        "comparison_results",
        "local_transcript_scan",
        "transcript_targets",
        "run_metadata",
    ]
    result = pd.read_excel(
        workbook,
        sheet_name="comparison_results",
        keep_default_na=False,
    ).iloc[0]
    assert result["result"] == "exact_match"
    assert result["region_start"] == 1
    assert result["region_end"] == 4
    assert result["differences"] == "None"


def test_private_panel_cli_returns_partial_status_when_one_target_fails(
    tmp_path,
    monkeypatch,
) -> None:
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    (cache_dir / "NM_000001.1.fasta").write_text(
        ">NM_000001.1 ready target\nGCAU\n",
        encoding="utf-8",
    )
    workbook = tmp_path / "partial.xlsx"
    monkeypatch.setattr(
        ncbi_blast.sys,
        "argv",
        [
            "ncbi_blast",
            "--as-sequence",
            "AUGC",
            "--private-panel",
            "--target-accession",
            "NM_000001.1",
            "--target-accession",
            "NM_000002.2",
            "--cache-dir",
            str(cache_dir),
            "--offline",
            "--result-workbook",
            str(workbook),
        ],
    )

    exit_code = ncbi_blast.main()

    assert exit_code == 2
    summaries = pd.read_excel(workbook, sheet_name="query_target_summary")
    assert summaries["scan_status"].tolist() == ["matched", "target_error"]


def test_download_targets_only_does_not_require_or_read_a_guide(tmp_path, monkeypatch) -> None:
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    (cache_dir / "NM_000001.1.fasta").write_text(
        ">NM_000001.1 cached target\nGCAU\n",
        encoding="utf-8",
    )
    workbook = tmp_path / "targets_only.xlsx"
    monkeypatch.setattr(
        ncbi_blast.sys,
        "argv",
        [
            "ncbi_blast",
            "--download-targets-only",
            "--target-accession",
            "NM_000001.1",
            "--cache-dir",
            str(cache_dir),
            "--offline",
            "--result-workbook",
            str(workbook),
        ],
    )

    exit_code = ncbi_blast.main()

    assert exit_code == 0
    assert pd.ExcelFile(workbook).sheet_names == [
        "input_queries",
        "transcript_targets",
        "local_transcript_scan",
        "run_metadata",
    ]
    targets = pd.read_excel(workbook, sheet_name="transcript_targets")
    assert targets.loc[0, "status"] == "ready"


def test_gui_args_support_private_ss_panel_without_remote_blast(tmp_path) -> None:
    input_file = tmp_path / "guides.xlsx"
    settings = {
        "sequence_type": "SS",
        "as_column": "sense",
        "as_name_column": "id",
        "target_accession": ["NM_000001.1", "NM_000002.2"],
        "target_accession_column": None,
        "target_file": None,
        "target_table": None,
        "private_panel": True,
        "offline": True,
        "email": "colleague@example.com",
        "scan_regions": ["full"],
        "max_mismatches": 0,
    }

    args = gui_args(input_file, "Sheet1", settings)
    validate_runtime_args(args)

    assert args.as_table is None
    assert args.ss_table == input_file
    assert args.ss_column == "sense"
    assert args.private_panel is True
    assert args.offline is True
    assert args.email == "colleague@example.com"
    assert args.blast is False
    assert args.result_workbook.name == "guides_private_transcript_panel_results.xlsx"
    assert args.cache_dir == shared_gui_transcript_cache_dir()


def test_single_sequence_gui_args_use_presets_shared_cache_and_local_workflow(tmp_path) -> None:
    settings = {
        "sequence_type": "AS",
        "sequence_name": "AS_1",
        "sequence": "A" * 18,
        "target_accession": "NM_000001.1",
        "scan_regions": ["full", "seed:2-8", "core:2-18"],
        "max_mismatches": 2,
        "closest": 5,
        "email": "colleague@example.com",
        "refresh_targets": False,
    }

    args = single_sequence_gui_args(settings)
    validate_runtime_args(args)

    assert args.as_sequence == "A" * 18
    assert args.as_name == "AS_1"
    assert args.ss_sequence is None
    assert args.target_accession == ["NM_000001.1"]
    assert args.scan_region == ["full", "seed:2-8", "core:2-18"]
    assert args.closest == 5
    assert args.email == "colleague@example.com"
    assert args.cache_dir == shared_gui_transcript_cache_dir()
    assert args.private_panel is True
    assert args.offline is False
    assert args.blast is False
    assert args.result_workbook is None


def test_single_sequence_gui_draft_retains_session_values_without_mutating_source(
    tmp_path,
) -> None:
    previous = {
        "sequence_type": "SS",
        "sequence_name": "SS_lead",
        "sequence": "AUGC",
        "target_mode": "paste",
        "target_accession": "NM_000001.1",
        "target_name": "manual transcript",
        "target_sequence": ">manual\nAUGC",
        "target_file": str(tmp_path / "target.fasta"),
        "scan_regions": ["full", "seed:2-8"],
        "max_mismatches": 1,
        "closest": 3,
        "refresh_targets": True,
    }

    draft = single_sequence_gui_draft(previous)
    draft["scan_regions"].append("core:2-18")

    assert draft["sequence_type"] == "SS"
    assert draft["sequence"] == "AUGC"
    assert draft["target_mode"] == "paste"
    assert draft["target_accession"] == "NM_000001.1"
    assert draft["target_name"] == "manual transcript"
    assert draft["target_sequence"] == ">manual\nAUGC"
    assert draft["target_file"] == str(tmp_path / "target.fasta")
    assert previous["scan_regions"] == ["full", "seed:2-8"]


def test_single_sequence_pasted_target_runs_locally_without_email_or_cache(tmp_path) -> None:
    args = single_sequence_gui_args(
        {
            "sequence_type": "SS",
            "sequence_name": "SS_1",
            "sequence": "AUGC",
            "target_mode": "paste",
            "target_name": "manual transcript",
            "target_sequence": "CCCAUGCUUU",
            "scan_regions": ["full"],
            "max_mismatches": 0,
            "closest": 2,
            "refresh_targets": True,
            "cache_dir": tmp_path / "cache",
            "email": "",
        }
    )
    validate_runtime_args(args)

    queries, scan_regions, result = run_single_sequence_scan(args)
    text = format_single_sequence_scan_result(args, queries, scan_regions, result)

    assert args.private_panel is False
    assert args.refresh_targets is False
    assert args.target_accession is None
    assert args.target_file is None
    assert result.targets[0].transcript_name == "manual transcript"
    assert result.targets[0].cache_status == "pasted sequence"
    assert result.targets[0].cache_path == ""
    assert result.summaries[0].exact_match_count == 1
    assert "Transcript source: pasted sequence" in text
    assert "NCBI transcript retrieval: Not used" in text
    assert "Transcript accession:" not in text
    assert not (tmp_path / "cache").exists()


def test_single_sequence_file_target_runs_locally_and_reports_file(tmp_path) -> None:
    target_file = tmp_path / "manual_target.fasta"
    target_file.write_text(">local target\nGCAU\n", encoding="utf-8")
    args = single_sequence_gui_args(
        {
            "sequence_type": "AS",
            "sequence_name": "AS_1",
            "sequence": "AUGC",
            "target_mode": "file",
            "target_file": target_file,
            "scan_regions": ["full"],
            "max_mismatches": 0,
            "closest": 1,
            "cache_dir": tmp_path / "cache",
            "email": "",
        }
    )
    validate_runtime_args(args)

    queries, scan_regions, result = run_single_sequence_scan(args)
    text = format_single_sequence_scan_result(args, queries, scan_regions, result)

    assert args.private_panel is False
    assert result.targets[0].transcript_name == "local target"
    assert result.targets[0].cache_status == "local file"
    assert result.targets[0].cache_path == str(target_file)
    assert result.summaries[0].exact_match_count == 1
    assert f"Target file: {target_file}" in text
    assert "NCBI transcript retrieval: Not used" in text


def test_single_sequence_gui_scan_reuses_cache_and_formats_direct_text_output(
    tmp_path,
) -> None:
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    (cache_dir / "NM_000001.1.fasta").write_text(
        ">NM_000001.1 cached transcript\n" + "A" * 25 + "\n",
        encoding="utf-8",
    )
    args = single_sequence_gui_args(
        {
            "sequence_type": "AS",
            "sequence_name": "AS_1",
            "sequence": "A" * 18,
            "target_accession": "NM_000001.1",
            "scan_regions": ["full", "seed:2-8", "core:2-18"],
            "max_mismatches": 0,
            "closest": 2,
            "refresh_targets": False,
            "cache_dir": cache_dir,
            "email": "test@example.com",
        }
    )
    validate_runtime_args(args)

    queries, scan_regions, result = run_single_sequence_scan(args)
    text = format_single_sequence_scan_result(args, queries, scan_regions, result)

    assert result.targets[0].cache_status == "cache"
    assert {
        region: sum(match.scan_region == region for match in result.closest_matches)
        for region in ("full", "seed", "core")
    } == {"full": 2, "seed": 2, "core": 2}
    assert "LOCAL SINGLE-SEQUENCE TRANSCRIPT SCAN" in text
    assert "Transcript source: cache" in text
    assert "Guide sequence sent to NCBI: No" in text
    assert "Region: full" in text
    assert "Region: seed" in text
    assert "Region: core" in text


def test_single_sequence_gui_scan_downloads_when_cache_is_missing(tmp_path) -> None:
    class RecordingClient:
        def __init__(self) -> None:
            self.requests: list[dict[str, object]] = []

        def get_text(self, _url: str, params: dict[str, object]) -> str:
            self.requests.append(params)
            return ">NM_000001.1 downloaded transcript\nGCAU\n"

    cache_dir = tmp_path / "cache"
    args = single_sequence_gui_args(
        {
            "sequence_type": "AS",
            "sequence_name": "AS_1",
            "sequence": "AUGC",
            "target_accession": "NM_000001.1",
            "scan_regions": ["full"],
            "max_mismatches": 0,
            "closest": 2,
            "refresh_targets": False,
            "cache_dir": cache_dir,
            "email": "test@example.com",
        }
    )
    client = RecordingClient()

    _queries, _scan_regions, result = run_single_sequence_scan(args, client=client)

    assert result.targets[0].cache_status == "downloaded"
    assert (cache_dir / "NM_000001.1.fasta").exists()
    assert [request["id"] for request in client.requests] == ["NM_000001.1"]
    assert all("QUERY" not in request for request in client.requests)


def test_refresh_targets_cannot_be_combined_with_offline_mode() -> None:
    args = build_parser().parse_args(
        [
            "--as-sequence",
            "AUGC",
            "--private-panel",
            "--target-accession",
            "NM_000001.1",
            "--offline",
            "--refresh-targets",
        ]
    )

    with pytest.raises(ValueError, match="cannot be combined"):
        validate_runtime_args(args)
