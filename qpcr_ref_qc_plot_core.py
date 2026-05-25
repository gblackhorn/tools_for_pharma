"""Plot reference-gene Mean CT QC summaries from extracted qPCR sheets."""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import sys

import pandas as pd

from excel_utils import list_excel_sheets
from qpcr_common import (
    MEAN_CT_COLUMN,
    REFERENCE_GENE_COLUMN,
    REFQC_SHEET_PREFIX,
    SAMPLE_SIZE_COLUMN,
    SEM_COLUMN,
    clean_text,
    is_refqc_sheet_name,
)


INPUT_FILE = Path("qpcr_result.xlsx")

REQUIRED_COLUMNS = [
    "Group",
    "Compound ID",
    REFERENCE_GENE_COLUMN,
    SAMPLE_SIZE_COLUMN,
    MEAN_CT_COLUMN,
    SEM_COLUMN,
]

FONT_FAMILY = "Arial"
FONT_SIZE = 11
TITLE_FONT_SIZE = 14
AXIS_LABEL_FONT_SIZE = 13
TICK_LABEL_FONT_SIZE = 11
LEGEND_FONT_SIZE = 11
FIGURE_HEIGHT = 7.2
LAYOUT_TOP = 0.84
LEGEND_TOP = 0.91
MAX_LEGEND_COLUMNS = 4

TEXT_COLOR = "#222222"
AXIS_COLOR = "#444444"
GRID_COLOR = "#D9D9D9"
ERROR_BAR_COLOR = "#2A2A2A"
REFERENCE_PALETTE = [
    "#4E79A7",
    "#59A14F",
    "#F28E2B",
    "#B07AA1",
    "#E15759",
    "#76B7B2",
]


def group_label(row: pd.Series) -> str:
    return f"{clean_text(row['Group'])} | {clean_text(row['Compound ID'])}"


def default_plot_dir(input_file: Path) -> Path:
    return input_file.with_name(input_file.stem)


def refqc_sheets(input_file: Path) -> list[str]:
    return [
        sheet
        for sheet in list_excel_sheets(input_file)
        if is_refqc_sheet_name(sheet)
    ]


def resolve_refqc_sheet(input_file: Path, sheet_name: str | None) -> str:
    if sheet_name is not None:
        if not is_refqc_sheet_name(sheet_name):
            raise ValueError(
                f"Sheet '{sheet_name}' is not a qPCR reference-QC sheet. "
                f"Choose a sheet whose name starts with '{REFQC_SHEET_PREFIX}'."
            )
        return sheet_name

    candidates = refqc_sheets(input_file)
    if not candidates:
        raise ValueError(
            f"No qPCR reference-QC sheet found in {input_file}. "
            f"Run qpcr_ref_qc_extract.py first, then choose a sheet whose name "
            f"starts with '{REFQC_SHEET_PREFIX}'."
        )
    return candidates[-1]


def validate_refqc_columns(summary: pd.DataFrame, sheet_name: str) -> None:
    missing_columns = [
        column for column in REQUIRED_COLUMNS
        if column not in summary.columns
    ]
    if missing_columns:
        missing = ", ".join(missing_columns)
        raise ValueError(
            f"Reference-QC sheet '{sheet_name}' is missing required columns: "
            f"{missing}."
        )


def read_refqc_summary(input_file: Path, sheet_name: str | None = None) -> pd.DataFrame:
    sheet = resolve_refqc_sheet(input_file, sheet_name)
    summary = pd.read_excel(input_file, sheet_name=sheet)
    validate_refqc_columns(summary, sheet)
    summary[MEAN_CT_COLUMN] = pd.to_numeric(summary[MEAN_CT_COLUMN], errors="coerce")
    summary[SEM_COLUMN] = pd.to_numeric(summary[SEM_COLUMN], errors="coerce")
    summary = summary.dropna(subset=[MEAN_CT_COLUMN, SEM_COLUMN])
    if summary.empty:
        raise ValueError("The selected reference-QC sheet has no plottable Mean CT rows.")
    summary["Plot label"] = summary.apply(group_label, axis=1)
    return summary


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


def finish_plot(axis, title: str, labels: list[str], output_path: Path) -> list[Path]:
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
    axis.set_ylabel(MEAN_CT_COLUMN)
    axis.yaxis.grid(True, color=GRID_COLOR, linewidth=0.7)
    axis.set_axisbelow(True)
    axis.set_xticks(range(len(labels)))
    axis.set_xticklabels(
        labels,
        rotation=90,
        ha="center",
        fontsize=TICK_LABEL_FONT_SIZE,
    )
    axis.tick_params(axis="both", length=3, width=0.8)
    axis.margins(x=0.02)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_linewidth(0.8)
    axis.spines["bottom"].set_linewidth(0.8)
    axis.figure.tight_layout(rect=(0, 0, 1, LAYOUT_TOP))
    png_path.parent.mkdir(parents=True, exist_ok=True)
    axis.figure.savefig(png_path, dpi=300, bbox_inches="tight")
    axis.figure.savefig(svg_path, bbox_inches="tight")
    plt.close(axis.figure)
    return [png_path, svg_path]


def plot_reference_mean_ct(summary: pd.DataFrame, plot_dir: Path) -> list[Path]:
    plt = get_pyplot()
    labels = list(dict.fromkeys(summary["Plot label"].tolist()))
    reference_genes = list(dict.fromkeys(summary[REFERENCE_GENE_COLUMN].tolist()))
    if not labels or not reference_genes:
        return []

    figure_width = max(12, len(labels) * 0.56)
    figure, axis = plt.subplots(figsize=(figure_width, FIGURE_HEIGHT))
    offset_step = min(0.18, 0.8 / max(len(reference_genes), 1))
    first_offset = -offset_step * (len(reference_genes) - 1) / 2

    plotted_y_values = []
    plotted_y_errors = []
    for gene_index, reference_gene in enumerate(reference_genes):
        subset = summary[summary[REFERENCE_GENE_COLUMN] == reference_gene]
        by_label = {
            label: row
            for label, row in subset.set_index("Plot label").iterrows()
        }
        x_values = [
            index + first_offset + gene_index * offset_step
            for index, label in enumerate(labels)
            if label in by_label
        ]
        y_values = [
            by_label[label][MEAN_CT_COLUMN]
            for label in labels
            if label in by_label
        ]
        y_errors = [
            by_label[label][SEM_COLUMN]
            for label in labels
            if label in by_label
        ]
        plotted_y_values.extend(y_values)
        plotted_y_errors.extend(y_errors)
        axis.errorbar(
            x_values,
            y_values,
            yerr=y_errors,
            fmt="o",
            markersize=5.5,
            capsize=3,
            elinewidth=0.9,
            capthick=0.9,
            linestyle="none",
            label=reference_gene,
            color=REFERENCE_PALETTE[gene_index % len(REFERENCE_PALETTE)],
            ecolor=ERROR_BAR_COLOR,
        )

    if plotted_y_values:
        lows = [
            value - error
            for value, error in zip(plotted_y_values, plotted_y_errors)
        ]
        highs = [
            value + error
            for value, error in zip(plotted_y_values, plotted_y_errors)
        ]
        axis.set_ylim(min(lows) - 0.5, max(highs) + 0.5)

    handles, legend_labels = axis.get_legend_handles_labels()
    figure.legend(
        handles,
        legend_labels,
        frameon=False,
        ncols=min(len(reference_genes), MAX_LEGEND_COLUMNS),
        loc="upper center",
        bbox_to_anchor=(0.5, LEGEND_TOP),
        handlelength=1.8,
        columnspacing=1.2,
    )

    return finish_plot(
        axis,
        "Reference gene Mean CT by group",
        labels,
        plot_dir / "reference_gene_mean_ct_qc.png",
    )


def create_refqc_plots(
    input_file: Path,
    output_dir: Path | None = None,
    sheet_name: str | None = None,
) -> list[Path]:
    summary = read_refqc_summary(input_file, sheet_name)
    plot_dir = output_dir if output_dir else default_plot_dir(input_file)
    return plot_reference_mean_ct(summary, plot_dir)


def choose_sheet_gui(root, input_file: Path) -> str | None:
    import tkinter as tk
    from tkinter import ttk

    sheets = list_excel_sheets(input_file)
    if len(sheets) <= 1:
        return None

    refqc = [sheet for sheet in sheets if is_refqc_sheet_name(sheet)]
    default_sheet = refqc[-1] if refqc else sheets[-1]
    selected = {"value": default_sheet}
    window = tk.Toplevel(root)
    window.title("Select reference QC sheet")
    window.resizable(False, False)
    window.columnconfigure(1, weight=1)

    ttk.Label(window, text="Worksheet").grid(
        row=0, column=0, padx=16, pady=(16, 8), sticky="w"
    )
    sheet_var = tk.StringVar(value=default_sheet)
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


def run_gui() -> int:
    import tkinter as tk
    from tkinter import filedialog, messagebox

    root = tk.Tk()
    root.withdraw()

    try:
        input_path = filedialog.askopenfilename(
            title="Select qPCR reference-QC Excel file",
            filetypes=[
                ("Excel files", "*.xlsx *.xlsm"),
                ("All files", "*.*"),
            ],
        )
        if not input_path:
            return 0

        input_file = Path(input_path)
        sheet_name = choose_sheet_gui(root, input_file)
        if sheet_name is None and len(list_excel_sheets(input_file)) > 1:
            return 0

        plot_paths = create_refqc_plots(input_file, None, sheet_name)
        message = f"Created {len(plot_paths)} reference-QC plot files."
        if plot_paths:
            message += f"\n\nSaved to:\n{default_plot_dir(input_file)}"
        messagebox.showinfo("Done", message)
        return 0
    except Exception as error:
        messagebox.showerror("qPCR reference-QC plotting failed", str(error))
        return 1
    finally:
        root.destroy()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Plot qPCR reference-gene Mean CT QC summaries."
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        default=INPUT_FILE,
        help=f"Input .xlsx file with a refqc sheet. Defaults to {INPUT_FILE}.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        help=(
            "Folder for plot PNG and SVG files. Defaults to a subfolder beside "
            "the input workbook with the same name as the workbook stem."
        ),
    )
    parser.add_argument(
        "--sheet",
        help=(
            "Worksheet containing extracted reference QC data. Defaults to the last "
            f"sheet whose name starts with '{REFQC_SHEET_PREFIX}'."
        ),
    )
    parser.add_argument(
        "--gui",
        action="store_true",
        help="Choose the reference-QC file and sheet with dialogs.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.gui:
        return run_gui()

    try:
        plot_paths = create_refqc_plots(args.input, args.output_dir, args.sheet)
    except Exception as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1
    if not plot_paths:
        print("No reference-QC plots were created.")
        return 0

    print("Reference-QC plots:")
    for plot_path in plot_paths:
        print(f"  {plot_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
