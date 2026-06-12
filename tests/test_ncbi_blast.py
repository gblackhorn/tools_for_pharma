from __future__ import annotations

from tools_for_pharma.oligo.ncbi_blast import (
    AntisenseRegion,
    AntisenseQuery,
    batch_antisense_queries,
    input_query_rows,
    parse_blast_field,
    parse_plain_antisense_lines,
    parse_scan_region,
    read_antisense_queries,
    scan_antisense_against_transcript,
    transcript_matches_to_csv,
)


def test_parse_blast_field_reads_rid_and_rtoe() -> None:
    text = "    RID = ABC123\n    RTOE = 42\n"

    assert parse_blast_field(text, "RID") == "ABC123"
    assert parse_blast_field(text, "RTOE") == "42"


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
