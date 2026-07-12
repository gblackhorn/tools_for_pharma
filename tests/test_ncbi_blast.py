from __future__ import annotations

from tools_for_pharma.oligo.ncbi_blast import (
    AntisenseRegion,
    AntisenseQuery,
    batch_antisense_queries,
    closest_transcript_matches,
    format_closest_transcript_matches_for_terminal,
    format_transcript_matches_for_terminal,
    input_query_rows,
    parse_blast_field,
    parse_plain_antisense_lines,
    parse_scan_region,
    read_antisense_queries,
    scan_antisense_against_transcript,
    scan_sense_against_transcript,
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
