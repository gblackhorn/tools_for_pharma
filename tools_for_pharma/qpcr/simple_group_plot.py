"""Create grouped bar plots from simple Group and mean/SEM tables."""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import sys

import pandas as pd

from tools_for_pharma.qpcr.common import clean_text
from tools_for_pharma.shared.excel_utils import list_excel_sheets


INPUT_FILE = Path("group_plot.xlsx")
PLOT_MODE_ALL = "all"
PLOT_MODE_ALL_VARIABLES = "all-variables"
PLOT_MODE_BY_COMPOUND = "by-compound"
PLOT_MODE_BY_DOSE = "by-dose"
PLOT_MODE_LEGACY = "legacy"
PLOT_MODE_CHOICES = [
    PLOT_MODE_ALL,
    PLOT_MODE_ALL_VARIABLES,
    PLOT_MODE_BY_COMPOUND,
    PLOT_MODE_BY_DOSE,
    PLOT_MODE_LEGACY,
]

FONT_FAMILY = "Arial"
FONT_SIZE = 11
TITLE_FONT_SIZE = 14
AXIS_LABEL_FONT_SIZE = 13
TICK_LABEL_FONT_SIZE = 11
LEGEND_FONT_SIZE = 11
FIGURE_HEIGHT = 6.5
LAYOUT_TOP = 0.86
LEGEND_TOP = 0.93
MAX_LEGEND_COLUMNS = 5

TEXT_COLOR = "#222222"
AXIS_COLOR = "#444444"
GRID_COLOR = "#D9D9D9"
BAR_EDGE_COLOR = "#303030"
ERROR_BAR_COLOR = "#2A2A2A"
BAR_PALETTE = [
    "#4E79A7",
    "#59A14F",
    "#F28E2B",
    "#B07AA1",
    "#E15759",
    "#76B7B2",
]

NUMBER_PATTERN = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
MEAN_SEM_PATTERN = re.compile(
    rf"^\s*(?P<mean>{NUMBER_PATTERN})\s*(?:\u00b1|\+/-|\+-)\s*"
    rf"(?P<sem>{NUMBER_PATTERN})\s*$"
)


def get_pyplot():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": [FONT_FAMILY, "DejaVu Sans", "Calibri", "Arial"],
            "font.size": FONT_SIZE,
            "axes.titlesize": TITLE_FONT_SIZE,
            "axes.labelsize": AXIS_LABEL_FONT_SIZE,
            "axes.titleweight": "bold",
            "axes.edgecolor": AXIS_COLOR,
            "axes.labelcolor": TEXT_COLOR,
            "xtick.color": TEXT_COLOR,
            "ytick.color": TEXT_COLOR,
            "legend.fontsize": LEGEND_FONT_SIZE,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )
    return plt


def default_plot_dir(input_file: Path) -> Path:
    return input_file.with_name(input_file.stem)


def read_raw_table(input_file: Path, sheet_name: str | int | None) -> pd.DataFrame:
    suffix = input_file.suffix.lower()
    if suffix in {".xlsx", ".xlsm", ".xls"}:
        return pd.read_excel(input_file, sheet_name=sheet_name or 0, header=None)
    if suffix in {".csv", ".txt"}:
        return pd.read_csv(input_file, header=None)
    raise ValueError("Input file must be an Excel workbook or CSV/text table.")


def find_header_row(raw: pd.DataFrame, header_row: int | None = None) -> int:
    if header_row is not None:
        return header_row

    max_rows = min(len(raw), 20)
    for row_index in range(max_rows):
        values = [clean_text(value).lower() for value in raw.iloc[row_index].tolist()]
        has_group = any(value == "group" for value in values)
        has_dose = any(value.startswith("dose") for value in values)
        has_time = any(value.startswith("time") for value in values)
        has_mean_sem = any("mean" in value and "sem" in value for value in values)
        if has_group and (has_mean_sem or has_time or has_dose):
            return row_index
    return 0


def detect_title(raw: pd.DataFrame, header_index: int) -> str | None:
    for row_index in range(header_index - 1, -1, -1):
        values = [
            clean_text(value)
            for value in raw.iloc[row_index].tolist()
            if clean_text(value)
        ]
        if values:
            return " ".join(dict.fromkeys(values))
    return None


def read_group_table(
    input_file: Path,
    sheet_name: str | int | None = None,
    header_row: int | None = None,
) -> pd.DataFrame:
    raw = read_raw_table(input_file, sheet_name)
    header_index = find_header_row(raw, header_row)
    if header_index >= len(raw):
        raise ValueError(f"Header row {header_index} is outside the input table.")

    columns = [clean_text(value) or f"Column {index + 1}" for index, value in enumerate(raw.iloc[header_index])]
    table = raw.iloc[header_index + 1 :].copy()
    table.columns = columns
    table = table.dropna(how="all")
    if len(table.columns) < 2:
        raise ValueError("Input table must contain at least two columns.")
    table.attrs["title"] = detect_title(raw, header_index)
    return table


def parse_mean_sem(value: object) -> tuple[float, float]:
    text = clean_text(value)
    match = MEAN_SEM_PATTERN.match(text)
    if not match:
        raise ValueError(
            f"Could not parse mean/SEM value '{text}'. Expected formats like "
            "'0.72 +/- 0.13', '0.72 +- 0.13', or the plus-minus symbol."
        )
    return float(match.group("mean")), float(match.group("sem"))


def split_group_label(label: object, delimiter: str = "-") -> tuple[str, str]:
    text = clean_text(label)
    if delimiter and delimiter in text:
        prefix, rest = text.split(delimiter, 1)
        return clean_text(prefix), clean_text(rest) or text

    parts = text.split(maxsplit=1)
    if len(parts) == 2:
        return parts[0], parts[1]
    return text, text


def prepare_plot_data(
    table: pd.DataFrame,
    group_column: str | None = None,
    mean_sem_column: str | None = None,
    delimiter: str = "-",
) -> pd.DataFrame:
    group_column = group_column or str(table.columns[0])
    mean_sem_column = mean_sem_column or str(table.columns[1])
    missing = [
        column
        for column in [group_column, mean_sem_column]
        if column not in table.columns
    ]
    if missing:
        raise ValueError(f"Input table is missing required columns: {', '.join(missing)}.")

    records = []
    for _, row in table.iterrows():
        raw_group = clean_text(row[group_column])
        raw_value = clean_text(row[mean_sem_column])
        if not raw_group or not raw_value:
            continue
        mean, sem = parse_mean_sem(raw_value)
        outer_group, condition = split_group_label(raw_group, delimiter)
        records.append(
            {
                "Group": raw_group,
                "Outer group": outer_group,
                "Condition": condition,
                "Mean": mean,
                "SEM": sem,
            }
        )

    summary = pd.DataFrame.from_records(records)
    if summary.empty:
        raise ValueError("No plottable rows were found in the input table.")
    return summary


def is_wide_time_table(table: pd.DataFrame) -> bool:
    columns = [clean_text(column).lower() for column in table.columns]
    return (
        any(column.startswith("dose") for column in columns)
        and any(column == "group" or "compound" in column for column in columns)
        and any(column.startswith("time") for column in columns)
    )


def detect_wide_column(
    table: pd.DataFrame,
    requested_column: str | None,
    label: str,
    predicate,
) -> str:
    if requested_column:
        if requested_column not in table.columns:
            raise ValueError(f"Input table is missing {label} column: {requested_column}.")
        return requested_column

    for column in table.columns:
        if predicate(clean_text(column).lower()):
            return str(column)
    raise ValueError(f"Could not detect the {label} column.")


def dose_unit_from_column(column: str) -> str:
    match = re.search(r"\(([^)]+)\)", column)
    return clean_text(match.group(1)) if match else ""


def format_dose(value: object, dose_column: str) -> str:
    text = clean_text(value)
    unit = dose_unit_from_column(dose_column)
    if unit and unit.lower() not in text.lower():
        return f"{text} {unit}"
    return text


def timepoint_from_column(column: str) -> str:
    text = clean_text(column)
    if "-" in text:
        return clean_text(text.split("-", 1)[1]) or text
    if ":" in text:
        return clean_text(text.split(":", 1)[1]) or text
    parts = text.split(maxsplit=1)
    return parts[1] if len(parts) == 2 and parts[0].lower() == "time" else text


def prepare_wide_time_data(
    table: pd.DataFrame,
    dose_column: str | None = None,
    compound_column: str | None = None,
) -> pd.DataFrame:
    dose_column = detect_wide_column(
        table,
        dose_column,
        "dose",
        lambda value: value.startswith("dose"),
    )
    compound_column = detect_wide_column(
        table,
        compound_column,
        "compound/group",
        lambda value: value == "group" or "compound" in value,
    )
    time_columns = [
        str(column)
        for column in table.columns
        if clean_text(column).lower().startswith("time")
    ]
    if not time_columns:
        raise ValueError("Could not detect any timepoint columns such as 'Time-D8'.")

    records = []
    for _, row in table.iterrows():
        compound = clean_text(row[compound_column])
        dose = format_dose(row[dose_column], dose_column)
        if not compound or not dose:
            continue

        for time_column in time_columns:
            raw_value = clean_text(row[time_column])
            if not raw_value:
                continue
            mean, sem = parse_mean_sem(raw_value)
            records.append(
                {
                    "Compound": compound,
                    "Dose": dose,
                    "Timepoint": timepoint_from_column(time_column),
                    "Compound dose": f"{compound} | {dose}",
                    "Mean": mean,
                    "SEM": sem,
                }
            )

    summary = pd.DataFrame.from_records(records)
    if summary.empty:
        raise ValueError("No plottable timepoint rows were found in the input table.")
    return summary


def finish_plot(axis, title: str, y_label: str, output_path: Path) -> list[Path]:
    plt = get_pyplot()
    png_path = output_path.with_suffix(".png")
    svg_path = output_path.with_suffix(".svg")
    axis.figure.suptitle(
        title,
        y=0.99,
        color=TEXT_COLOR,
        fontsize=TITLE_FONT_SIZE,
        fontweight="bold",
    )
    axis.set_ylabel(y_label)
    axis.yaxis.grid(True, color=GRID_COLOR, linewidth=0.7)
    axis.set_axisbelow(True)
    axis.tick_params(axis="both", length=3, width=0.8)
    axis.margins(x=0.08)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_linewidth(0.8)
    axis.spines["bottom"].set_linewidth(0.8)
    axis.set_ylim(bottom=0)
    axis.figure.tight_layout(rect=(0, 0, 1, LAYOUT_TOP))
    png_path.parent.mkdir(parents=True, exist_ok=True)
    axis.figure.savefig(png_path, dpi=300, bbox_inches="tight")
    axis.figure.savefig(svg_path, bbox_inches="tight")
    plt.close(axis.figure)
    return [png_path, svg_path]


def plot_grouped_summary(
    summary: pd.DataFrame,
    x_column: str,
    series_column: str,
    output_path: Path,
    title: str,
    y_label: str,
    x_label: str | None = None,
) -> list[Path]:
    plt = get_pyplot()
    x_labels = list(dict.fromkeys(summary[x_column].tolist()))
    series_labels = list(dict.fromkeys(summary[series_column].tolist()))
    if not x_labels or not series_labels:
        return []

    figure_width = max(7.5, len(x_labels) * max(len(series_labels), 1) * 0.42)
    figure, axis = plt.subplots(figsize=(figure_width, FIGURE_HEIGHT))
    bar_width = min(0.2, 0.82 / max(len(series_labels), 1))
    first_offset = -bar_width * (len(series_labels) - 1) / 2

    for series_index, series_label in enumerate(series_labels):
        subset = summary[summary[series_column] == series_label]
        by_x_label = {
            row[x_column]: row
            for _, row in subset.iterrows()
        }
        x_values = [
            x_index + first_offset + series_index * bar_width
            for x_index, x_label in enumerate(x_labels)
            if x_label in by_x_label
        ]
        y_values = [
            by_x_label[x_label]["Mean"]
            for x_label in x_labels
            if x_label in by_x_label
        ]
        y_errors = [
            by_x_label[x_label]["SEM"]
            for x_label in x_labels
            if x_label in by_x_label
        ]
        axis.bar(
            x_values,
            y_values,
            width=bar_width,
            yerr=y_errors,
            capsize=3,
            label=series_label,
            color=BAR_PALETTE[series_index % len(BAR_PALETTE)],
            edgecolor=BAR_EDGE_COLOR,
            error_kw={
                "elinewidth": 0.9,
                "ecolor": ERROR_BAR_COLOR,
                "capthick": 0.9,
            },
            linewidth=0.35,
        )

    axis.set_xlabel(x_label or x_column)
    axis.set_xticks(range(len(x_labels)))
    axis.set_xticklabels(x_labels, fontsize=TICK_LABEL_FONT_SIZE)
    handles, legend_labels = axis.get_legend_handles_labels()
    axis.figure.legend(
        handles,
        legend_labels,
        frameon=False,
        ncols=min(len(series_labels), MAX_LEGEND_COLUMNS),
        loc="upper center",
        bbox_to_anchor=(0.5, LEGEND_TOP),
        handlelength=1.8,
        columnspacing=1.2,
    )
    return finish_plot(axis, title, y_label, output_path)


def plot_grouped_bars(
    summary: pd.DataFrame,
    output_path: Path,
    title: str,
    y_label: str,
) -> list[Path]:
    plt = get_pyplot()
    outer_groups = list(dict.fromkeys(summary["Outer group"].tolist()))
    conditions = list(dict.fromkeys(summary["Condition"].tolist()))
    if not outer_groups or not conditions:
        return []

    figure_width = max(7.5, len(outer_groups) * max(len(conditions), 1) * 0.48)
    figure, axis = plt.subplots(figsize=(figure_width, FIGURE_HEIGHT))
    bar_width = min(0.22, 0.8 / max(len(conditions), 1))
    first_offset = -bar_width * (len(conditions) - 1) / 2

    for condition_index, condition in enumerate(conditions):
        subset = summary[summary["Condition"] == condition]
        by_group = {
            row["Outer group"]: row
            for _, row in subset.iterrows()
        }
        x_values = [
            group_index + first_offset + condition_index * bar_width
            for group_index, group in enumerate(outer_groups)
            if group in by_group
        ]
        y_values = [
            by_group[group]["Mean"]
            for group in outer_groups
            if group in by_group
        ]
        y_errors = [
            by_group[group]["SEM"]
            for group in outer_groups
            if group in by_group
        ]
        axis.bar(
            x_values,
            y_values,
            width=bar_width,
            yerr=y_errors,
            capsize=3,
            label=condition,
            color=BAR_PALETTE[condition_index % len(BAR_PALETTE)],
            edgecolor=BAR_EDGE_COLOR,
            error_kw={
                "elinewidth": 0.9,
                "ecolor": ERROR_BAR_COLOR,
                "capthick": 0.9,
            },
            linewidth=0.35,
        )

    axis.set_xticks(range(len(outer_groups)))
    axis.set_xticklabels(outer_groups, fontsize=TICK_LABEL_FONT_SIZE)
    handles, legend_labels = axis.get_legend_handles_labels()
    figure.legend(
        handles,
        legend_labels,
        frameon=False,
        ncols=min(len(conditions), MAX_LEGEND_COLUMNS),
        loc="upper center",
        bbox_to_anchor=(0.5, LEGEND_TOP),
        handlelength=1.8,
        columnspacing=1.2,
    )
    return finish_plot(axis, title, y_label, output_path)


def create_wide_time_plots(
    summary: pd.DataFrame,
    plot_dir: Path,
    title: str,
    y_label: str,
    plot_mode: str,
) -> list[Path]:
    output_paths = []

    if plot_mode in {PLOT_MODE_ALL, PLOT_MODE_ALL_VARIABLES}:
        output_paths.extend(
            plot_grouped_summary(
                summary,
                "Timepoint",
                "Compound dose",
                plot_dir / "time_compound_dose_grouped_bar_plot.png",
                f"{title} - all compounds and doses",
                y_label,
                "Timepoint",
            )
        )

    if plot_mode in {PLOT_MODE_ALL, PLOT_MODE_BY_COMPOUND}:
        for compound, subset in summary.groupby("Compound", sort=False):
            output_paths.extend(
                plot_grouped_summary(
                    subset,
                    "Timepoint",
                    "Dose",
                    plot_dir / f"{sanitize_filename(compound)}_dose_by_time.png",
                    f"{title} - {compound}: dose by time",
                    y_label,
                    "Timepoint",
                )
            )

    if plot_mode in {PLOT_MODE_ALL, PLOT_MODE_BY_DOSE}:
        for dose, subset in summary.groupby("Dose", sort=False):
            output_paths.extend(
                plot_grouped_summary(
                    subset,
                    "Timepoint",
                    "Compound",
                    plot_dir / f"{sanitize_filename(dose)}_compound_by_time.png",
                    f"{title} - {dose}: compound by time",
                    y_label,
                    "Timepoint",
                )
            )

    return output_paths


def sanitize_filename(value: object) -> str:
    name = re.sub(r"[^A-Za-z0-9._-]+", "_", clean_text(value)).strip("_")
    return name or "plot"


def create_group_plot(
    input_file: Path,
    output_dir: Path | None = None,
    sheet_name: str | int | None = None,
    header_row: int | None = None,
    group_column: str | None = None,
    mean_sem_column: str | None = None,
    delimiter: str = "-",
    title: str | None = None,
    y_label: str = "Mean +/- SEM",
    plot_mode: str = PLOT_MODE_ALL,
    dose_column: str | None = None,
) -> list[Path]:
    table = read_group_table(input_file, sheet_name, header_row)
    plot_dir = output_dir if output_dir else default_plot_dir(input_file)
    plot_title = title or table.attrs.get("title") or input_file.stem

    if plot_mode != PLOT_MODE_LEGACY and is_wide_time_table(table):
        summary = prepare_wide_time_data(table, dose_column, group_column)
        return create_wide_time_plots(summary, plot_dir, plot_title, y_label, plot_mode)

    summary = prepare_plot_data(table, group_column, mean_sem_column, delimiter)
    return plot_grouped_bars(
        summary,
        plot_dir / "grouped_bar_plot.png",
        plot_title,
        y_label,
    )


def choose_sheet_gui(root, input_file: Path) -> str | None:
    import tkinter as tk
    from tkinter import ttk

    if input_file.suffix.lower() not in {".xlsx", ".xlsm", ".xls"}:
        return None

    sheets = list_excel_sheets(input_file)
    if len(sheets) <= 1:
        return None

    selected = {"value": sheets[0]}
    window = tk.Toplevel(root)
    window.title("Select data sheet")
    window.resizable(False, False)
    window.columnconfigure(1, weight=1)

    ttk.Label(window, text="Worksheet").grid(
        row=0, column=0, padx=16, pady=(16, 8), sticky="w"
    )
    sheet_var = tk.StringVar(value=sheets[0])
    sheet_box = ttk.Combobox(
        window,
        textvariable=sheet_var,
        values=sheets,
        state="readonly",
        width=max(30, min(60, max(len(sheet) for sheet in sheets) + 2)),
    )
    sheet_box.grid(row=0, column=1, padx=16, pady=(16, 8), sticky="ew")

    buttons = ttk.Frame(window)
    buttons.grid(row=1, column=0, columnspan=2, padx=16, pady=(8, 16), sticky="e")

    def use_sheet() -> None:
        selected["value"] = sheet_var.get()
        window.destroy()

    def cancel() -> None:
        selected["value"] = None
        window.destroy()

    ttk.Button(buttons, text="Cancel", command=cancel).grid(row=0, column=0, padx=(0, 8))
    ttk.Button(buttons, text="Continue", command=use_sheet).grid(row=0, column=1)
    window.protocol("WM_DELETE_WINDOW", cancel)
    window.bind("<Return>", lambda _event: use_sheet())
    window.bind("<Escape>", lambda _event: cancel())
    window.grab_set()
    sheet_box.focus_set()
    window.wait_window()
    return selected["value"]


def choose_plot_mode_gui(root) -> str | None:
    import tkinter as tk
    from tkinter import ttk

    labels = {
        PLOT_MODE_ALL: "All comparison plots",
        PLOT_MODE_ALL_VARIABLES: "All variables together",
        PLOT_MODE_BY_COMPOUND: "Same compound: doses by time",
        PLOT_MODE_BY_DOSE: "Same dose: compounds by time",
        PLOT_MODE_LEGACY: "Old two-column group plot",
    }
    selected = {"value": PLOT_MODE_ALL}
    window = tk.Toplevel(root)
    window.title("Select plot style")
    window.resizable(False, False)
    window.columnconfigure(1, weight=1)

    ttk.Label(window, text="Plot style").grid(
        row=0, column=0, padx=16, pady=(16, 8), sticky="w"
    )
    mode_var = tk.StringVar(value=labels[PLOT_MODE_ALL])
    mode_box = ttk.Combobox(
        window,
        textvariable=mode_var,
        values=[labels[mode] for mode in PLOT_MODE_CHOICES],
        state="readonly",
        width=36,
    )
    mode_box.grid(row=0, column=1, padx=16, pady=(16, 8), sticky="ew")

    buttons = ttk.Frame(window)
    buttons.grid(row=1, column=0, columnspan=2, padx=16, pady=(8, 16), sticky="e")

    def use_mode() -> None:
        reverse_labels = {label: mode for mode, label in labels.items()}
        selected["value"] = reverse_labels[mode_var.get()]
        window.destroy()

    def cancel() -> None:
        selected["value"] = None
        window.destroy()

    ttk.Button(buttons, text="Cancel", command=cancel).grid(row=0, column=0, padx=(0, 8))
    ttk.Button(buttons, text="Continue", command=use_mode).grid(row=0, column=1)
    window.protocol("WM_DELETE_WINDOW", cancel)
    window.bind("<Return>", lambda _event: use_mode())
    window.bind("<Escape>", lambda _event: cancel())
    window.grab_set()
    mode_box.focus_set()
    window.wait_window()
    return selected["value"]


def run_gui() -> int:
    import tkinter as tk
    from tkinter import filedialog, messagebox

    root = tk.Tk()
    root.withdraw()

    try:
        input_path = filedialog.askopenfilename(
            title="Select simple grouped-bar data file",
            filetypes=[
                ("Excel and CSV files", "*.xlsx *.xlsm *.xls *.csv *.txt"),
                ("Excel files", "*.xlsx *.xlsm *.xls"),
                ("CSV files", "*.csv"),
                ("All files", "*.*"),
            ],
        )
        if not input_path:
            return 0

        input_file = Path(input_path)
        sheet_name = choose_sheet_gui(root, input_file)
        is_excel = input_file.suffix.lower() in {".xlsx", ".xlsm", ".xls"}
        if sheet_name is None and is_excel:
            sheets = list_excel_sheets(input_file)
            if len(sheets) > 1:
                return 0

        plot_mode = choose_plot_mode_gui(root)
        if plot_mode is None:
            return 0

        plot_paths = create_group_plot(
            input_file,
            sheet_name=sheet_name,
            plot_mode=plot_mode,
        )
        output_dir = default_plot_dir(input_file)
        message = f"Created {len(plot_paths)} grouped bar plot files."
        if plot_paths:
            message += f"\n\nSaved to:\n{output_dir}"
        messagebox.showinfo("Done", message)
        return 0
    except Exception as error:
        messagebox.showerror("Grouped bar plot failed", str(error))
        return 1
    finally:
        root.destroy()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Create a grouped bar plot from Group and mean/SEM columns."
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        default=INPUT_FILE,
        help=f"Input .xlsx/.csv table. Defaults to {INPUT_FILE}.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        help="Folder for plot PNG and SVG files. Defaults to a folder beside the input.",
    )
    parser.add_argument("--sheet", help="Excel worksheet name. Defaults to the first sheet.")
    parser.add_argument(
        "--header-row",
        type=int,
        help="Zero-based header row index. Defaults to auto-detecting Group and SEM headers.",
    )
    parser.add_argument(
        "--group-column",
        help=(
            "Group/compound column name. Defaults to the first detected table "
            "column for old two-column tables or the 'Group' column for wide tables."
        ),
    )
    parser.add_argument(
        "--dose-column",
        help="Dose column name for wide timepoint tables. Defaults to the 'Dose...' column.",
    )
    parser.add_argument(
        "--mean-sem-column",
        help="Mean/SEM column name. Defaults to the second detected table column.",
    )
    parser.add_argument(
        "--delimiter",
        default="-",
        help="Separator between outer group and condition. Defaults to '-'.",
    )
    parser.add_argument(
        "--title",
        help="Plot title. Defaults to the input filename.",
    )
    parser.add_argument(
        "--y-label",
        default="Mean +/- SEM",
        help="Y-axis label. Defaults to 'Mean +/- SEM'.",
    )
    parser.add_argument(
        "--plot-mode",
        choices=PLOT_MODE_CHOICES,
        default=PLOT_MODE_ALL,
        help=(
            "For wide tables: all = create every comparison; all-variables = "
            "timepoints with compound+dose bars; by-compound = one plot per "
            "compound with dose bars; by-dose = one plot per dose with compound "
            "bars; legacy = force old two-column plotting."
        ),
    )
    parser.add_argument(
        "--gui",
        action="store_true",
        help="Choose the data Excel/CSV file with dialogs.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.gui:
        return run_gui()

    try:
        plot_paths = create_group_plot(
            args.input,
            args.output_dir,
            args.sheet,
            args.header_row,
            args.group_column,
            args.mean_sem_column,
            args.delimiter,
            args.title,
            args.y_label,
            args.plot_mode,
            args.dose_column,
        )
    except Exception as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1

    if not plot_paths:
        print("No grouped bar plots were created.")
        return 0

    print("Grouped bar plot:")
    for plot_path in plot_paths:
        print(f"  {plot_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
