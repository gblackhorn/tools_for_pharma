from pathlib import Path

import pandas as pd

from tools_for_pharma.qpcr.common import (
    ANIMAL_ID_COLUMN,
    INDIVIDUAL_RQ_COLUMN,
    MEAN_RQ_COLUMN,
    REFERENCE_SOURCE_COLUMN,
    SAMPLE_ID_COLUMN,
    SAMPLE_SIZE_COLUMN,
    SEM_COLUMN,
)
from tools_for_pharma.qpcr.extract import (
    COLUMN_ROLE_ANIMAL_ID,
    COLUMN_ROLE_COMPOUND,
    COLUMN_ROLE_GROUP,
    COLUMN_ROLE_INDIVIDUAL,
    COLUMN_ROLE_MEAN,
    COLUMN_ROLE_SAMPLE_ID,
    COLUMN_ROLE_SEM,
    column_signature,
    detected_mapping_columns,
    extraction_preview,
    mapping_table_from_cells,
    summarize_aggregate_table,
    summarize_full_table,
    summarize_qpcr_region,
)
from tools_for_pharma.qpcr.plot import (
    AXIS_SPINE_OFFSET_POINTS,
    PLOT_MODE_SPLIT,
    bar_summary,
    create_plots,
    finish_plot,
    get_pyplot,
    group_label,
    identifier_label,
    jitter_positions,
)


def make_full_qpcr_table(include_individual_values: bool = True) -> pd.DataFrame:
    values = [
        [
            "Sample ID",
            "Group",
            "Animal ID",
            "Compound ID",
            "HTT1a",
            "GAPDH",
            "Relative to control group",
            "MEAN RQ",
            "Sample size (n)",
            "SEM",
            "MEAN RQ ± SEM",
        ],
        [None, None, None, None, "Mean CT", "Mean CT", None, None, None, None, None],
    ]
    individual_values = [0.9, 1.0, 1.1, 1.0] if include_individual_values else [None] * 4
    for index, individual_value in enumerate(individual_values, start=1):
        values.append(
            [
                index,
                "G1" if index == 1 else None,
                f"A{index}",
                "saline" if index == 1 else None,
                26.0,
                23.0,
                individual_value,
                1.0 if index == 1 else None,
                4 if index == 1 else None,
                0.04 if index == 1 else None,
                "1.00 ± 0.04" if index == 1 else None,
            ]
        )
    return pd.DataFrame(
        values,
        index=range(1, len(values) + 1),
        columns=range(1, len(values[0]) + 1),
    )


def make_alternate_qpcr_table() -> pd.DataFrame:
    top_headers = [
        "Sample ID",
        "Group",
        "Animal ID",
        "Compound ID",
        "APOE",
        "APOE",
        "APOE",
        "GAPDH",
        "GAPDH",
        "GAPDH",
        "ΔCT",
        "Relative RQ",
        "Mean control RQ",
        "Normalized RQ",
        "Mean",
        "SEM",
    ]
    bottom_headers = [
        "Sample ID",
        "Group",
        "Animal ID",
        "Compound ID",
        "CT (Duplicate well 1)",
        "CT (Duplicate well 2)",
        "Mean CT",
        "CT (Duplicate well 1)",
        "CT (Duplicate well 2)",
        "Mean CT",
        "ΔCT",
        "Relative RQ",
        "Mean control RQ",
        "Normalized RQ",
        "Mean",
        "SEM",
    ]
    values = [top_headers, bottom_headers]
    normalized_values = [0.92, 1.03, 1.08, 0.97]
    for index, normalized_value in enumerate(normalized_values, start=1):
        values.append(
            [
                index,
                "G1" if index == 1 else None,
                f"A{index}",
                "saline" if index == 1 else None,
                26.0,
                26.1,
                26.05,
                23.0,
                23.1,
                23.05,
                3.0,
                normalized_value * 1.1,
                1.1,
                normalized_value,
                1.0 if index == 1 else None,
                0.035 if index == 1 else None,
            ]
        )
    return pd.DataFrame(
        values,
        index=range(1, len(values) + 1),
        columns=range(1, len(top_headers) + 1),
    )


def test_full_table_extraction_keeps_one_row_per_animal() -> None:
    summary = summarize_full_table(make_full_qpcr_table(), table_number=1)

    assert len(summary) == 4
    assert summary[REFERENCE_SOURCE_COLUMN].tolist() == ["GAPDH"] * 4
    assert summary[SAMPLE_ID_COLUMN].tolist() == [1, 2, 3, 4]
    assert summary[ANIMAL_ID_COLUMN].tolist() == ["A1", "A2", "A3", "A4"]
    assert summary[INDIVIDUAL_RQ_COLUMN].tolist() == ["0.9", "1", "1.1", "1"]
    assert summary[MEAN_RQ_COLUMN].tolist() == ["1"] * 4
    assert summary[SEM_COLUMN].tolist() == ["0.04"] * 4


def test_full_table_without_individual_values_stays_one_summary_row() -> None:
    summary = summarize_full_table(
        make_full_qpcr_table(include_individual_values=False),
        table_number=1,
    )

    assert len(summary) == 1
    assert summary.iloc[0][MEAN_RQ_COLUMN] == "1"
    assert summary.iloc[0][INDIVIDUAL_RQ_COLUMN] == ""


def test_alternate_headers_auto_detect_normalized_rq_and_final_mean() -> None:
    table = make_alternate_qpcr_table()

    summary = summarize_full_table(table, table_number=1)
    detected = detected_mapping_columns(table)

    assert len(summary) == 4
    assert summary[REFERENCE_SOURCE_COLUMN].tolist() == ["GAPDH"] * 4
    assert summary[INDIVIDUAL_RQ_COLUMN].tolist() == ["0.92", "1.03", "1.08", "0.97"]
    assert summary[MEAN_RQ_COLUMN].tolist() == ["1"] * 4
    assert summary[SEM_COLUMN].tolist() == ["0.035"] * 4
    assert summary[SAMPLE_SIZE_COLUMN].tolist() == [4] * 4
    assert detected[COLUMN_ROLE_INDIVIDUAL] == 14
    assert detected[COLUMN_ROLE_MEAN] == 15
    assert detected[COLUMN_ROLE_SEM] == 16


def test_manual_mapping_supports_nonstandard_column_names() -> None:
    table = make_alternate_qpcr_table().copy()
    rename = {
        1: "Specimen",
        2: "Treatment cohort",
        3: "Subject",
        4: "Molecule",
        14: "Dot value",
        15: "Average value",
        16: "Std error",
    }
    for column, header in rename.items():
        table.loc[1, column] = header
        table.loc[2, column] = header

    mapping = {
        COLUMN_ROLE_SAMPLE_ID: column_signature(table, 1),
        COLUMN_ROLE_GROUP: column_signature(table, 2),
        COLUMN_ROLE_ANIMAL_ID: column_signature(table, 3),
        COLUMN_ROLE_COMPOUND: column_signature(table, 4),
        COLUMN_ROLE_INDIVIDUAL: column_signature(table, 14),
        COLUMN_ROLE_MEAN: column_signature(table, 15),
        COLUMN_ROLE_SEM: column_signature(table, 16),
    }

    summary = summarize_qpcr_region(table, mapping)

    assert len(summary) == 4
    assert summary["Group"].tolist() == ["G1"] * 4
    assert summary[INDIVIDUAL_RQ_COLUMN].tolist() == ["0.92", "1.03", "1.08", "0.97"]


def test_mapping_table_finds_embedded_alternate_header() -> None:
    table = make_alternate_qpcr_table()
    blank_rows = pd.DataFrame(
        [[None] * len(table.columns)] * 4,
        index=range(1, 5),
        columns=table.columns,
    )
    embedded = pd.concat([blank_rows, table.set_axis(range(5, 11), axis=0)])

    detected_table, header_row, table_count = mapping_table_from_cells(embedded)

    assert header_row == 5
    assert table_count == 1
    assert detected_mapping_columns(detected_table)[COLUMN_ROLE_MEAN] == 15


def test_extraction_preview_counts_values_and_checks_group_mean() -> None:
    summary = summarize_full_table(make_alternate_qpcr_table(), table_number=1)

    bar_count, individual_count, issues = extraction_preview(summary)

    assert bar_count == 1
    assert individual_count == 4
    assert issues == []

    summary.loc[0, INDIVIDUAL_RQ_COLUMN] = "2"
    _bar_count, _individual_count, issues = extraction_preview(summary)
    assert len(issues) == 1


def test_aggregate_extraction_keeps_geomean_per_animal() -> None:
    metadata = make_full_qpcr_table()
    aggregate = pd.DataFrame(
        [
            ["Geomean", "MEAN RQ", "Sample size (n)", "SEM", "MEAN RQ ± SEM"],
            [None, None, None, None, None],
            [0.95, 1.0, 4, 0.03, "1.00 ± 0.03"],
            [1.02, None, None, None, None],
            [1.08, None, None, None, None],
            [0.95, None, None, None, None],
        ],
        index=range(1, 7),
        columns=range(1, 6),
    )

    summary = summarize_aggregate_table(
        aggregate,
        metadata,
        "Geomean (GAPDH)",
        table_number=2,
    )

    assert len(summary) == 4
    assert summary[SAMPLE_ID_COLUMN].tolist() == [1, 2, 3, 4]
    assert summary[ANIMAL_ID_COLUMN].tolist() == ["A1", "A2", "A3", "A4"]
    assert summary[INDIVIDUAL_RQ_COLUMN].tolist() == ["0.95", "1.02", "1.08", "0.95"]


def test_plotter_collapses_animal_rows_and_creates_scatter_overlay(tmp_path: Path) -> None:
    rows = []
    for group, compound, mean, values in [
        ("G1", "saline", 1.0, [0.9, 1.1]),
        ("G2", "AD001", 0.7, [0.65, 0.75]),
    ]:
        for sample_index, value in enumerate(values, start=1):
            rows.append(
                {
                    "Group": group,
                    "Compound ID": compound,
                    REFERENCE_SOURCE_COLUMN: "GAPDH",
                    SAMPLE_SIZE_COLUMN: 2,
                    MEAN_RQ_COLUMN: mean,
                    SEM_COLUMN: 0.05,
                    SAMPLE_ID_COLUMN: sample_index,
                    ANIMAL_ID_COLUMN: f"{group}-A{sample_index}",
                    INDIVIDUAL_RQ_COLUMN: value,
                }
            )
    extracted = pd.DataFrame(rows)
    extracted["Plot label"] = extracted.apply(
        lambda row: f"{row['Group']} | {row['Compound ID']}",
        axis=1,
    )
    assert len(bar_summary(extracted)) == 2
    assert jitter_positions(2.0, 4, 0.12) == [1.88, 1.96, 2.04, 2.12]

    input_file = tmp_path / "qpcr_points.xlsx"
    extracted.drop(columns="Plot label").to_excel(
        input_file,
        sheet_name="plotdata-qPCR",
        index=False,
    )
    paths = create_plots(
        input_file,
        output_dir=tmp_path / "plots",
        plot_mode=PLOT_MODE_SPLIT,
        sheet_name="plotdata-qPCR",
    )

    assert len(paths) == 2
    assert all(path.exists() and path.stat().st_size > 0 for path in paths)


def test_identifier_labels_drop_only_numeric_trailing_zero() -> None:
    assert identifier_label(68.0) == "68"
    assert identifier_label(68.5) == "68.5"
    assert identifier_label("68.0") == "68.0"
    assert identifier_label("068") == "068"
    assert identifier_label("G68") == "G68"
    assert group_label(pd.Series({"Group": 68.0, "Compound ID": "AD001"})) == "68 | AD001"


def test_finished_plot_has_detached_spines_and_vertical_labels(tmp_path: Path) -> None:
    plt = get_pyplot()
    _figure, axis = plt.subplots()
    axis.bar([0, 1], [1.0, 0.8])

    finish_plot(axis, "Test", ["68 | saline", "G2 | AD001"], tmp_path / "spines.png")

    assert axis.spines["left"].get_position() == (
        "outward",
        AXIS_SPINE_OFFSET_POINTS,
    )
    assert axis.spines["bottom"].get_position() == (
        "outward",
        AXIS_SPINE_OFFSET_POINTS,
    )
    assert [label.get_rotation() for label in axis.get_xticklabels()] == [90.0, 90.0]
