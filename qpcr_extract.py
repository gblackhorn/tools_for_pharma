"""qPCR Excel table extraction logic."""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import sys
from typing import Optional

import pandas as pd

from excel_utils import (
    append_or_replace_sheet,
    list_excel_sheets,
    read_excel_cells_with_merged_values,
    sanitize_sheet_name,
    write_table,
)
from qpcr_common import (
    CONTROL_COMPOUNDS,
    MEAN_RQ_COLUMN,
    PLOTDATA_SHEET_PREFIX,
    REFERENCE_SOURCE_COLUMN,
    SAMPLE_SIZE_COLUMN,
    SEM_COLUMN,
    clean_text,
    is_blank,
    result_columns,
)


INPUT_FILE = Path("qpcr_result.xlsx")
OUTPUT_FILE = Path("qpcr_summary.xlsx")

SHEET_NAME = None
HEADER_ROWS = [1, 2]
TOP_GENE_ROW = 1
DATA_START_ROW = 3


def parse_number(value: object) -> Optional[float]:
    if is_blank(value):
        return None
    try:
        return float(clean_text(value).replace(",", ""))
    except ValueError:
        return None


def format_number(value: object) -> str:
    number = parse_number(value)
    if number is None:
        return clean_text(value)
    return f"{number:g}"


def split_mean_rq_sem(value: object) -> tuple[Optional[str], Optional[str]]:
    numbers = re.findall(
        r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?",
        clean_text(value),
    )
    if len(numbers) < 2:
        return None, None
    return format_number(numbers[0]), format_number(numbers[1])


def header_text(data: pd.DataFrame, column: int) -> str:
    return " ".join(
        clean_text(data.loc[row, column])
        for row in HEADER_ROWS
        if row in data.index
    ).lower()


def find_column(
    data: pd.DataFrame,
    include_any: Optional[list[str]] = None,
    include_all: Optional[list[str]] = None,
    exclude_any: Optional[list[str]] = None,
) -> Optional[int]:
    include_any = [item.lower() for item in include_any or []]
    include_all = [item.lower() for item in include_all or []]
    exclude_any = [item.lower() for item in exclude_any or []]

    for column in data.columns:
        text = header_text(data, column)
        if include_all and not all(item in text for item in include_all):
            continue
        if include_any and not any(item in text for item in include_any):
            continue
        if exclude_any and any(item in text for item in exclude_any):
            continue
        return column
    return None


def text_tokens(text: str) -> set[str]:
    return set(re.findall(r"[a-z0-9]+", text.lower()))


def find_sample_size_column(data: pd.DataFrame) -> Optional[int]:
    for column in data.columns:
        text = header_text(data, column)
        tokens = text_tokens(text)
        if "sample size" in text or "\u7ec4\u5185\u6837\u672c\u91cf" in text or "n" in tokens:
            if "sample id" not in text:
                return column
    return None


def find_geomean_column(data: pd.DataFrame) -> Optional[int]:
    return find_column(data, include_any=["geomean"])


def require_column(name: str, column: Optional[int]) -> int:
    if column is None:
        raise ValueError(f"Could not find required column: {name}")
    return column


def is_empty_column(data: pd.DataFrame, column: int) -> bool:
    return all(is_blank(value) for value in data[column])


def split_side_by_side_tables(data: pd.DataFrame) -> list[pd.DataFrame]:
    """Split one worksheet into table blocks separated by fully empty columns."""
    tables = []
    current_columns = []
    for column in data.columns:
        if is_empty_column(data, column):
            if current_columns:
                tables.append(data.loc[:, current_columns])
                current_columns = []
            continue
        current_columns.append(column)
    if current_columns:
        tables.append(data.loc[:, current_columns])
    return tables


def row_text(data: pd.DataFrame, row: int) -> str:
    return " ".join(clean_text(data.loc[row, column]) for column in data.columns).lower()


def row_has_exact_cell(data: pd.DataFrame, row: int, text: str) -> bool:
    text = text.lower()
    return any(clean_text(data.loc[row, column]).lower() == text for column in data.columns)


def find_embedded_qpcr_header_row(data: pd.DataFrame) -> Optional[int]:
    """Find a qPCR table header row embedded inside a larger worksheet."""
    for row in data.index:
        text = row_text(data, row)
        if not row_has_exact_cell(data, row, "group"):
            continue
        if "compound id" not in text and "compound" not in text:
            continue
        if "mean rq" not in text or "sem" not in text:
            continue
        return row
    return None


def normalize_embedded_qpcr_region(data: pd.DataFrame) -> pd.DataFrame:
    header_row = find_embedded_qpcr_header_row(data)
    if header_row is None:
        return pd.DataFrame()

    region = data.loc[header_row:].copy()
    region.index = range(1, len(region) + 1)
    return region


def looks_like_gene_name(value: object) -> bool:
    text = clean_text(value)
    return bool(re.fullmatch(r"[A-Za-z][A-Za-z0-9_.-]{0,19}", text))


def detect_gene_names(data: pd.DataFrame) -> list[str]:
    genes = []
    seen = set()
    for column in data.columns:
        gene = clean_text(data.loc[TOP_GENE_ROW, column])
        if not looks_like_gene_name(gene):
            continue
        nearby_header = " ".join(
            clean_text(data.loc[row, column])
            for row in HEADER_ROWS
            if row != TOP_GENE_ROW
        ).lower()
        if "ct" not in nearby_header:
            continue
        if gene not in seen:
            genes.append(gene)
            seen.add(gene)
    return genes


def control_or_blank_row(
    row: pd.Series,
    group_col: Optional[int],
    compound_col: Optional[int],
) -> bool:
    group = row[group_col] if group_col is not None else ""
    compound = row[compound_col] if compound_col is not None else ""
    return (
        is_blank(group)
        or is_blank(compound)
        or clean_text(compound).upper() in CONTROL_COMPOUNDS
    )


def build_mean_rq(
    row: pd.Series,
    mean_rq_sem_col: Optional[int],
    mean_rq_col: Optional[int],
) -> Optional[str]:
    if mean_rq_col is not None and not is_blank(row[mean_rq_col]):
        return format_number(row[mean_rq_col])
    if mean_rq_sem_col is None or is_blank(row[mean_rq_sem_col]):
        return None
    mean_rq, _sem = split_mean_rq_sem(row[mean_rq_sem_col])
    return mean_rq


def build_sem(
    row: pd.Series,
    mean_rq_sem_col: Optional[int],
    sem_col: Optional[int],
) -> Optional[str]:
    if sem_col is not None and not is_blank(row[sem_col]):
        return format_number(row[sem_col])
    if mean_rq_sem_col is None or is_blank(row[mean_rq_sem_col]):
        return None
    _mean_rq, sem = split_mean_rq_sem(row[mean_rq_sem_col])
    return sem


def empty_summary() -> pd.DataFrame:
    return pd.DataFrame(columns=result_columns())


def summarize_full_table(data: pd.DataFrame, table_number: int) -> pd.DataFrame:
    genes = detect_gene_names(data)
    reference_source = "; ".join(genes[1:]) if len(genes) >= 2 else ""
    group_col = require_column(
        f"Group in table {table_number}",
        find_column(data, include_any=["group"]),
    )
    compound_col = find_column(data, include_any=["compound id", "compound"])
    sample_size_col = find_sample_size_column(data)
    mean_rq_sem_col = find_column(data, include_all=["mean", "rq", "sem"])
    mean_rq_col = find_column(
        data,
        include_all=["mean", "rq"],
        exclude_any=["sem", "control"],
    )
    sem_col = find_column(data, include_any=["sem"], exclude_any=["mean rq"])

    if mean_rq_sem_col is None and (mean_rq_col is None or sem_col is None):
        raise ValueError(
            "Could not find MEAN RQ +/- SEM, or separate MEAN RQ and SEM "
            f"columns in table {table_number}."
        )

    rows = data.loc[DATA_START_ROW:].copy()
    for column in [
        group_col,
        compound_col,
        sample_size_col,
        mean_rq_sem_col,
        mean_rq_col,
        sem_col,
    ]:
        if column is not None:
            rows[column] = rows[column].ffill()

    rows = rows[
        ~rows.apply(
            control_or_blank_row,
            axis=1,
            args=(group_col, compound_col),
        )
    ]

    summary = pd.DataFrame(index=rows.index)
    summary["Group"] = rows[group_col]
    summary["Compound ID"] = rows[compound_col] if compound_col is not None else ""
    summary[REFERENCE_SOURCE_COLUMN] = reference_source
    summary[SAMPLE_SIZE_COLUMN] = (
        rows[sample_size_col] if sample_size_col is not None else ""
    )
    summary[MEAN_RQ_COLUMN] = rows.apply(
        build_mean_rq,
        axis=1,
        args=(mean_rq_sem_col, mean_rq_col),
    )
    summary[SEM_COLUMN] = rows.apply(
        build_sem,
        axis=1,
        args=(mean_rq_sem_col, sem_col),
    )

    summary = summary[
        summary["Group"].map(lambda value: not is_blank(value))
        & summary[MEAN_RQ_COLUMN].map(lambda value: not is_blank(value))
        & summary[SEM_COLUMN].map(lambda value: not is_blank(value))
    ]
    return summary.drop_duplicates(subset=result_columns()).reset_index(drop=True)


def summarize_aggregate_table(
    data: pd.DataFrame,
    metadata_table: pd.DataFrame,
    reference_source: str,
    table_number: int,
) -> pd.DataFrame:
    group_col = find_column(metadata_table, include_any=["group"])
    compound_col = find_column(metadata_table, include_any=["compound id", "compound"])
    sample_size_col = find_sample_size_column(data)
    mean_rq_sem_col = find_column(data, include_all=["mean", "rq", "sem"])
    mean_rq_col = find_column(
        data,
        include_all=["mean", "rq"],
        exclude_any=["sem", "control"],
    )
    sem_col = find_column(data, include_any=["sem"], exclude_any=["mean rq"])

    if group_col is None or compound_col is None:
        raise ValueError(
            "Could not use aggregate table because the metadata table is missing "
            "Group or Compound ID."
        )
    if mean_rq_sem_col is None and (mean_rq_col is None or sem_col is None):
        raise ValueError(
            "Could not find MEAN RQ +/- SEM, or separate MEAN RQ and SEM "
            f"columns in aggregate table {table_number}."
        )

    rows = data.loc[DATA_START_ROW:].copy()
    metadata_rows = metadata_table.loc[DATA_START_ROW:].copy()
    summary = pd.DataFrame(index=rows.index)
    summary["Group"] = metadata_rows[group_col]
    summary["Compound ID"] = metadata_rows[compound_col]
    summary[REFERENCE_SOURCE_COLUMN] = reference_source
    summary[SAMPLE_SIZE_COLUMN] = (
        rows[sample_size_col] if sample_size_col is not None else ""
    )
    summary[MEAN_RQ_COLUMN] = rows.apply(
        build_mean_rq,
        axis=1,
        args=(mean_rq_sem_col, mean_rq_col),
    )
    summary[SEM_COLUMN] = rows.apply(
        build_sem,
        axis=1,
        args=(mean_rq_sem_col, sem_col),
    )

    summary = summary[
        summary["Group"].map(lambda value: not is_blank(value))
        & summary["Compound ID"].map(lambda value: not is_blank(value))
        & summary["Compound ID"].map(lambda value: clean_text(value).upper() not in CONTROL_COMPOUNDS)
        & summary[MEAN_RQ_COLUMN].map(lambda value: not is_blank(value))
        & summary[SEM_COLUMN].map(lambda value: not is_blank(value))
    ]
    return summary.drop_duplicates(subset=result_columns()).reset_index(drop=True)


def summarize_qpcr_region(data: pd.DataFrame) -> pd.DataFrame:
    tables = split_side_by_side_tables(data)
    if not tables:
        return empty_summary()

    full_tables = [
        (table_number, table)
        for table_number, table in enumerate(tables, start=1)
        if find_column(table, include_any=["group"]) is not None
    ]
    aggregate_tables = [
        (table_number, table)
        for table_number, table in enumerate(tables, start=1)
        if find_column(table, include_any=["group"]) is None
        and find_geomean_column(table) is not None
    ]

    summaries = [
        summarize_full_table(table, table_number)
        for table_number, table in full_tables
    ]
    if aggregate_tables and full_tables:
        reference_genes = []
        for _table_number, table in full_tables:
            reference_genes.extend(detect_gene_names(table)[1:])
        aggregate_source = f"Geomean ({'; '.join(dict.fromkeys(reference_genes))})"
        metadata_table = full_tables[0][1]
        summaries.extend(
            summarize_aggregate_table(
                table,
                metadata_table,
                aggregate_source,
                table_number,
            )
            for table_number, table in aggregate_tables
        )

    return pd.concat(summaries, ignore_index=True) if summaries else empty_summary()


def build_summary(
    input_file: Path,
    sheet_name: Optional[str] = SHEET_NAME,
) -> pd.DataFrame:
    data = read_excel_cells_with_merged_values(input_file, sheet_name)
    summary = summarize_qpcr_region(data)
    if not summary.empty:
        return summary

    embedded_region = normalize_embedded_qpcr_region(data)
    if embedded_region.empty:
        return empty_summary()
    return summarize_qpcr_region(embedded_region)


def default_output_path(input_file: Path) -> Path:
    return input_file.with_name(f"{input_file.stem}_qpcr_summary.xlsx")


def resolved_sheet_name(input_file: Path, sheet_name: Optional[str]) -> str:
    sheets = list_excel_sheets(input_file)
    if not sheets:
        raise ValueError(f"Workbook has no sheets: {input_file}")
    return sheet_name if sheet_name is not None else sheets[0]


def default_plotdata_sheet_name(input_file: Path, source_sheet_name: Optional[str]) -> str:
    return sanitize_sheet_name(
        f"{PLOTDATA_SHEET_PREFIX}{resolved_sheet_name(input_file, source_sheet_name)}"
    )


def default_summary_sheet_name(input_file: Path, source_sheet_name: Optional[str]) -> str:
    return default_plotdata_sheet_name(input_file, source_sheet_name)


def save_summary(
    summary: pd.DataFrame,
    input_file: Path,
    output_file: Optional[Path] = None,
    append_sheet: Optional[str] = None,
    output_sheet: Optional[str] = None,
) -> Path:
    if append_sheet:
        append_or_replace_sheet(input_file, summary, append_sheet)
        if output_file:
            write_table(summary, str(output_file), sheet_name=append_sheet)
            return output_file
        return input_file

    destination = output_file if output_file else default_output_path(input_file)
    return write_table(summary, str(destination), sheet_name=output_sheet)


def extract_summary(
    input_file: Path = INPUT_FILE,
    output_file: Optional[Path] = None,
    sheet_name: Optional[str] = SHEET_NAME,
    append_sheet: Optional[str] = "",
) -> pd.DataFrame:
    summary = build_summary(input_file, sheet_name)
    output_sheet = default_plotdata_sheet_name(input_file, sheet_name)
    if append_sheet == "":
        append_sheet = output_sheet
    save_summary(summary, input_file, output_file, append_sheet, output_sheet)
    return summary


def choose_sheet_gui(root, input_file: Path) -> Optional[str]:
    import tkinter as tk
    from tkinter import ttk

    sheets = list_excel_sheets(input_file)
    if len(sheets) <= 1:
        return None

    selected = {"value": sheets[0]}
    window = tk.Toplevel(root)
    window.title("Select sheet")
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


def run_gui() -> int:
    import tkinter as tk
    from tkinter import filedialog, messagebox, simpledialog

    root = tk.Tk()
    root.withdraw()

    try:
        input_path = filedialog.askopenfilename(
            title="Select qPCR Excel file",
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

        default_sheet_name = default_plotdata_sheet_name(input_file, sheet_name)
        append_sheet = simpledialog.askstring(
            "Plot data sheet name",
            "Sheet name:",
            initialvalue=default_sheet_name,
            parent=root,
        )
        if not append_sheet:
            return 0

        summary = extract_summary(input_file, None, sheet_name, append_sheet)
        messagebox.showinfo(
            "Done",
            f"Extracted {len(summary)} rows.\n\nSaved to:\n{input_file}\nSheet: {append_sheet}",
        )
        return 0
    except Exception as error:
        messagebox.showerror("qPCR extraction failed", str(error))
        return 1
    finally:
        root.destroy()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Extract qPCR plot summary values from an Excel table."
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        default=INPUT_FILE,
        help=f"Input .xlsx file. Defaults to {INPUT_FILE}.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help=(
            "Optional separate output .xlsx file. By default, extraction is saved "
            "as a plotdata sheet appended to the input workbook."
        ),
    )
    parser.add_argument(
        "--sheet",
        default=SHEET_NAME,
        help="Excel sheet name. Defaults to the first sheet.",
    )
    parser.add_argument(
        "--append-sheet",
        nargs="?",
        const="",
        default="",
        help=(
            "Append/replace the extracted summary as a sheet in the input workbook. "
            "Defaults to 'plotdata-[source sheet]' when no sheet name is provided."
        ),
    )
    parser.add_argument(
        "--separate-output",
        action="store_true",
        help="Save only to a separate output workbook instead of appending to input.",
    )
    parser.add_argument(
        "--gui",
        action="store_true",
        help="Choose the input and output Excel files with file dialogs.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.gui:
        return run_gui()

    append_sheet = None if args.separate_output else args.append_sheet
    output = args.output if args.output or args.separate_output else None
    try:
        summary = extract_summary(args.input, output, args.sheet, append_sheet)
    except Exception as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1
    if append_sheet is not None and not output:
        appended_sheet = append_sheet if append_sheet else default_plotdata_sheet_name(args.input, args.sheet)
        destination = f"{args.input} [{sanitize_sheet_name(appended_sheet)}]"
    else:
        destination = output if output else default_output_path(args.input)
    print(summary)
    print(f"Saved to: {destination}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
