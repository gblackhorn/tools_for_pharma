from __future__ import annotations

import pandas as pd

from tools_for_pharma.plotting.curve import (
    DEFAULT_IC_MARKER_TEXT,
    DEFAULT_IC_MARKERS,
    ERROR_COLUMN,
    FIT_METHOD_4PL,
    PLOT_MODE_BOTH,
    PLOT_MODE_COMBINED,
    calculate_curve_descriptor,
    calculate_fit_ic_markers,
    calculate_ic_markers,
    create_curve_plot,
    fit_curve,
    four_parameter_logistic,
    normalize_raw_table,
    parse_ic_markers,
    parse_response_value,
    prepare_curve_data,
)
from tools_for_pharma.plotting.ic50_summary import descriptor_lines


def test_default_ic_markers_are_drug_development_focused() -> None:
    assert DEFAULT_IC_MARKERS == [50.0, 75.0, 90.0]
    assert DEFAULT_IC_MARKER_TEXT == "IC50,IC75,IC90"
    assert parse_ic_markers("IC50,75,90") == [50.0, 75.0, 90.0]


def test_parse_response_value_accepts_mean_error_text() -> None:
    assert parse_response_value("94.50\u00b10.81") == (94.5, 0.81)
    assert parse_response_value("94.50 +/- 0.81") == (94.5, 0.81)
    assert parse_response_value("94.50 +- 0.81") == (94.5, 0.81)
    assert parse_response_value("94.50") == (94.5, None)


def test_prepare_curve_data_detects_columns_and_series() -> None:
    table = pd.DataFrame(
        {
            "Compound": ["A", "A", "B", "B"],
            "Concentration (nM)": [1, 10, 1, 10],
            "Inhibition rate (%)": [20, 80, 10, 55],
        }
    )

    summary, series_column, x_column, y_column = prepare_curve_data(table)

    assert series_column == "Compound"
    assert x_column == "Concentration (nM)"
    assert y_column == "Inhibition rate (%)"
    assert summary["Compound"].tolist() == ["A", "A", "B", "B"]


def test_prepare_curve_data_uses_wide_table_by_default() -> None:
    table = pd.DataFrame(
        {
            "Concentration (nM)": [0.1, 1, 10, 100],
            "AD-001": [5, 28, 61, 92],
            "AD-002": [3, 18, 49, 88],
        }
    )

    summary, series_column, x_column, y_column = prepare_curve_data(table)

    assert series_column == "Compound"
    assert x_column == "Concentration (nM)"
    assert y_column == "Inhibition rate (%)"
    assert summary["Compound"].tolist() == [
        "AD-001",
        "AD-001",
        "AD-001",
        "AD-001",
        "AD-002",
        "AD-002",
        "AD-002",
        "AD-002",
    ]
    assert summary["Inhibition rate (%)"].tolist() == [5, 28, 61, 92, 3, 18, 49, 88]


def test_prepare_curve_data_uses_error_from_wide_mean_error_cells() -> None:
    table = pd.DataFrame(
        {
            "Concentration (nM)": [0.1, 1, 10],
            "AD-001": ["5.00\u00b10.50", "28.0\u00b11.2", "61.0\u00b13.4"],
            "AD-002": ["3 +/- 0.2", "18 +/- 0.9", "49 +/- 2.1"],
        }
    )

    summary, series_column, x_column, y_column = prepare_curve_data(table)

    assert series_column == "Compound"
    assert x_column == "Concentration (nM)"
    assert y_column == "Inhibition rate (%)"
    assert summary["Inhibition rate (%)"].tolist() == [5.0, 28.0, 61.0, 3.0, 18.0, 49.0]
    assert summary[ERROR_COLUMN].tolist() == [0.5, 1.2, 3.4, 0.2, 0.9, 2.1]


def test_prepare_curve_data_uses_error_from_long_mean_error_cells() -> None:
    table = pd.DataFrame(
        {
            "Compound": ["A", "A", "B", "B"],
            "Concentration (nM)": [1, 10, 1, 10],
            "Inhibition rate (%)": ["20\u00b12", "80\u00b15", "10\u00b11", "55\u00b14"],
        }
    )

    summary, series_column, _x_column, y_column = prepare_curve_data(
        table,
        y_column="Inhibition rate (%)",
        series_column="Compound",
    )

    assert series_column == "Compound"
    assert y_column == "Inhibition rate (%)"
    assert summary["Inhibition rate (%)"].tolist() == [20.0, 80.0, 10.0, 55.0]
    assert summary[ERROR_COLUMN].tolist() == [2.0, 5.0, 1.0, 4.0]


def test_normalize_raw_table_treats_numeric_first_row_as_data() -> None:
    raw = pd.DataFrame(
        [
            [0.1, 5, 3],
            [1, 28, 18],
            [10, 61, 49],
        ]
    )

    table = normalize_raw_table(raw)
    summary, series_column, x_column, y_column = prepare_curve_data(table)

    assert table.columns.tolist() == ["Concentration", "Compound 1", "Compound 2"]
    assert series_column == "Compound"
    assert x_column == "Concentration"
    assert y_column == "Inhibition rate (%)"
    assert summary["Compound"].tolist() == [
        "Compound 1",
        "Compound 1",
        "Compound 1",
        "Compound 2",
        "Compound 2",
        "Compound 2",
    ]


def test_calculate_ic_markers_interpolates_on_log_concentration() -> None:
    summary = pd.DataFrame(
        {
            "Concentration (nM)": [1, 10, 100],
            "Inhibition rate (%)": [0, 50, 100],
        }
    )

    markers = calculate_ic_markers(
        summary,
        "Concentration (nM)",
        "Inhibition rate (%)",
        [25, 50, 75],
    )

    assert [marker.label for marker in markers] == ["IC25", "IC50", "IC75"]
    assert [round(marker.x, 3) if marker.x else None for marker in markers] == [
        3.162,
        10.0,
        31.623,
    ]
    assert [round(marker.slope, 1) if marker.slope else None for marker in markers] == [
        50.0,
        50.0,
        50.0,
    ]


def test_curve_descriptor_reports_range_auc_and_marker_slope() -> None:
    summary = pd.DataFrame(
        {
            "Concentration (nM)": [1, 10, 100],
            "Inhibition rate (%)": [0, 50, 100],
        }
    )

    descriptor = calculate_curve_descriptor(
        summary,
        "A",
        "Concentration (nM)",
        "Inhibition rate (%)",
        [50],
    )
    lines = descriptor_lines([descriptor], log_x=True)

    assert descriptor.series_label == "A"
    assert descriptor.min_response == 0
    assert descriptor.max_response == 100
    assert round(descriptor.auc, 1) == 100.0
    assert any("IC50: 10, slope 50 per log10 dose" in line for line in lines)


def test_four_pl_fit_recovers_smooth_dose_response_parameters() -> None:
    concentrations = pd.Series([0.1, 0.3, 1, 3, 10, 30, 100])
    responses = four_parameter_logistic(
        concentrations.to_numpy(dtype=float),
        bottom=2,
        top=98,
        log_ic50=1,
        hill_slope=1.2,
    )
    summary = pd.DataFrame(
        {
            "Concentration (nM)": concentrations,
            "Inhibition rate (%)": responses,
        }
    )

    fitted_curve = fit_curve(
        summary,
        "Concentration (nM)",
        "Inhibition rate (%)",
        FIT_METHOD_4PL,
    )
    markers = calculate_fit_ic_markers(
        summary,
        "Concentration (nM)",
        "Inhibition rate (%)",
        fitted_curve,
        [50],
    )
    descriptor = calculate_curve_descriptor(
        summary,
        "A",
        "Concentration (nM)",
        "Inhibition rate (%)",
        [50],
        fit_method=FIT_METHOD_4PL,
        curve_fit=fitted_curve,
    )

    assert fitted_curve.method == FIT_METHOD_4PL
    assert round(fitted_curve.hill_slope, 1) == 1.2
    assert round(markers[0].x, 1) == 10.0
    assert descriptor.r_squared is not None
    assert descriptor.r_squared > 0.999


def test_create_curve_plot_can_write_combined_and_single_outputs(tmp_path) -> None:
    input_file = tmp_path / "inhibition.csv"
    input_file.write_text(
        "\n".join(
            [
                "Concentration (nM),A,B",
                "0.1,5,3",
                "1,28,18",
                "10,61,49",
                "100,92,88",
            ]
        )
    )
    output_dir = tmp_path / "plots"

    paths = create_curve_plot(
        input_file,
        output_dir=output_dir,
        plot_mode=PLOT_MODE_BOTH,
    )

    assert len(paths) == 6
    assert (output_dir / "inhibition_dose_response_curve.png") in paths
    assert (output_dir / "inhibition_A_dose_response_curve.png") in paths
    assert (output_dir / "inhibition_B_dose_response_curve.png") in paths


def test_combined_summary_renderer_is_skipped_when_too_many_curves(tmp_path) -> None:
    input_file = tmp_path / "many.csv"
    input_file.write_text(
        "\n".join(
            [
                "Concentration (nM),A,B,C,D",
                "0.1,5,4,3,2",
                "1,28,25,22,19",
                "10,61,58,55,52",
                "100,92,90,88,86",
            ]
        )
    )
    calls = []

    def renderer(_axis, descriptors, _log_x, _colors) -> None:
        calls.append(len(descriptors))

    create_curve_plot(
        input_file,
        output_dir=tmp_path / "plots",
        plot_mode=PLOT_MODE_COMBINED,
        summary_renderer=renderer,
        summary_max_curves=3,
    )

    assert calls == []
