from __future__ import annotations

import pandas as pd

from tools_for_pharma.qpcr.simple_group_plot import (
    parse_mean_sem,
    prepare_plot_data,
    prepare_wide_time_data,
    read_group_table,
)


def test_parse_mean_sem_accepts_plus_minus_sign() -> None:
    assert parse_mean_sem("0.72 \u00b1 0.13") == (0.72, 0.13)


def test_prepare_plot_data_groups_by_prefix() -> None:
    table = pd.DataFrame(
        {
            "Group": [
                "G1-baseline",
                "G1-2mpk D33",
                "G1-5mpk D33",
                "G2-baseline",
                "G2-2mpk D33",
                "G2-5mpk D33",
            ],
            "(Mean SEM)": [
                "1 +/- 0",
                "0.72 +/- 0.13",
                "0.68 +/- 0.04",
                "1 +/- 0",
                "0.74 +/- 0.06",
                "0.58 +/- 0.1",
            ],
        }
    )

    summary = prepare_plot_data(table)

    assert summary["Outer group"].tolist() == ["G1", "G1", "G1", "G2", "G2", "G2"]
    assert summary["Condition"].tolist() == [
        "baseline",
        "2mpk D33",
        "5mpk D33",
        "baseline",
        "2mpk D33",
        "5mpk D33",
    ]
    assert summary["Mean"].tolist() == [1.0, 0.72, 0.68, 1.0, 0.74, 0.58]
    assert summary["SEM"].tolist() == [0.0, 0.13, 0.04, 0.0, 0.06, 0.1]


def test_read_group_table_detects_title_and_header_rows(tmp_path) -> None:
    input_file = tmp_path / "group_plot.xlsx"
    raw = pd.DataFrame(
        [
            ["MSH3 remaining on D33 relative to baseline in Liver", None],
            ["Group", "(Mean SEM)"],
            ["G1-baseline", "1 +/- 0"],
            ["G1-2mpk D33", "0.72 +/- 0.13"],
        ]
    )
    raw.to_excel(input_file, index=False, header=False)

    table = read_group_table(input_file)

    assert table.attrs["title"] == "MSH3 remaining on D33 relative to baseline in Liver"
    assert table.columns.tolist() == ["Group", "(Mean SEM)"]
    assert table["Group"].tolist() == ["G1-baseline", "G1-2mpk D33"]


def test_prepare_wide_time_data_unpivots_dose_group_and_timepoints() -> None:
    table = pd.DataFrame(
        {
            "Dose (mpk)": [2, 5],
            "Group": ["G1 AD12111", "G1 AD12111"],
            "Time-baseline": ["1 +/- 0", "1 +/- 0"],
            "Time-D8": ["0.4 +/- 0.06", "0.44 +/- 0.02"],
            "Time-D29": ["0.42 +/- 0.03", "0.34 +/- 0.02"],
        }
    )

    summary = prepare_wide_time_data(table)

    assert summary["Compound"].tolist() == [
        "G1 AD12111",
        "G1 AD12111",
        "G1 AD12111",
        "G1 AD12111",
        "G1 AD12111",
        "G1 AD12111",
    ]
    assert summary["Dose"].tolist() == [
        "2 mpk",
        "2 mpk",
        "2 mpk",
        "5 mpk",
        "5 mpk",
        "5 mpk",
    ]
    assert summary["Timepoint"].tolist() == [
        "baseline",
        "D8",
        "D29",
        "baseline",
        "D8",
        "D29",
    ]
    assert summary["Compound dose"].tolist()[0] == "G1 AD12111 | 2 mpk"
    assert summary["Mean"].tolist() == [1.0, 0.4, 0.42, 1.0, 0.44, 0.34]
    assert summary["SEM"].tolist() == [0.0, 0.06, 0.03, 0.0, 0.02, 0.02]
