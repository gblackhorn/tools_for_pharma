"""qPCR Excel table extraction logic."""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import sys
from typing import Optional

import pandas as pd

from tools_for_pharma.qpcr.common import (
    ANIMAL_ID_COLUMN,
    CONTROL_COMPOUNDS,
    INDIVIDUAL_RQ_COLUMN,
    MEAN_RQ_COLUMN,
    PLOTDATA_SHEET_PREFIX,
    REFERENCE_SOURCE_COLUMN,
    SAMPLE_ID_COLUMN,
    SAMPLE_SIZE_COLUMN,
    SEM_COLUMN,
    clean_text,
    is_blank,
    result_columns,
    summary_columns,
)
from tools_for_pharma.shared.excel_utils import (
    append_or_replace_sheet,
    list_excel_sheets,
    read_excel_cells_with_merged_values,
    sanitize_sheet_name,
    write_table,
)


INPUT_FILE = Path("qpcr_result.xlsx")
OUTPUT_FILE = Path("qpcr_summary.xlsx")

SHEET_NAME = None
HEADER_ROWS = [1, 2]
TOP_GENE_ROW = 1
DATA_START_ROW = 3

COLUMN_ROLE_GROUP = "group"
COLUMN_ROLE_COMPOUND = "compound"
COLUMN_ROLE_SAMPLE_ID = "sample_id"
COLUMN_ROLE_ANIMAL_ID = "animal_id"
COLUMN_ROLE_SAMPLE_SIZE = "sample_size"
COLUMN_ROLE_MEAN = "mean"
COLUMN_ROLE_SEM = "sem"
COLUMN_ROLE_INDIVIDUAL = "individual"

ColumnSignature = tuple[str, str]
ColumnMapping = dict[str, Optional[ColumnSignature]]


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


def column_signature(data: pd.DataFrame, column: int) -> ColumnSignature:
    values = [
        clean_text(data.loc[row, column]) if row in data.index else ""
        for row in HEADER_ROWS
    ]
    return values[0], values[1]


def normalized_column_signature(signature: ColumnSignature) -> ColumnSignature:
    return tuple(clean_text(value).lower() for value in signature)  # type: ignore[return-value]


def find_exact_header_column(
    data: pd.DataFrame,
    names: list[str],
) -> Optional[int]:
    expected = {clean_text(name).lower() for name in names}
    for column in data.columns:
        if any(
            clean_text(data.loc[row, column]).lower() in expected
            for row in HEADER_ROWS
            if row in data.index
        ):
            return column
    return None


def find_summary_mean_column(data: pd.DataFrame) -> Optional[int]:
    mean_rq_column = find_column(
        data,
        include_all=["mean", "rq"],
        exclude_any=["sem", "control", "ct"],
    )
    if mean_rq_column is not None:
        return mean_rq_column
    return find_exact_header_column(data, ["mean"])


def find_column_by_signature(
    data: pd.DataFrame,
    signature: ColumnSignature,
) -> Optional[int]:
    expected = normalized_column_signature(signature)
    exact_matches = [
        column
        for column in data.columns
        if normalized_column_signature(column_signature(data, column)) == expected
    ]
    if exact_matches:
        return exact_matches[0]

    meaningful = [value for value in expected if value]
    if not meaningful:
        return None
    fallback_matches = [
        column
        for column in data.columns
        if meaningful[-1]
        in normalized_column_signature(column_signature(data, column))
    ]
    return fallback_matches[0] if len(fallback_matches) == 1 else None


def resolve_mapped_column(
    data: pd.DataFrame,
    mapping: Optional[ColumnMapping],
    role: str,
    detected_column: Optional[int],
) -> Optional[int]:
    if mapping is None or role not in mapping:
        return detected_column
    signature = mapping[role]
    if signature is None:
        return None
    return find_column_by_signature(data, signature)


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


def find_individual_rq_column(data: pd.DataFrame) -> Optional[int]:
    relative_to_control = find_column(
        data,
        include_all=["relative", "control", "group"],
    )
    if relative_to_control is not None:
        return relative_to_control
    normalized_rq = find_exact_header_column(data, ["normalized rq", "normalised rq"])
    if normalized_rq is not None:
        return normalized_rq
    return find_exact_header_column(data, ["remaining", "individual rq"])


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
        has_summary_mean = "mean rq" in text or row_has_exact_cell(data, row, "mean")
        if not has_summary_mean or not row_has_exact_cell(data, row, "sem"):
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


def finalize_summary(summary: pd.DataFrame) -> pd.DataFrame:
    """Keep one row per animal, or one summary row when raw values are unavailable."""
    if summary.empty:
        return empty_summary()

    groups = []
    for _key, subset in summary.groupby(
        summary_columns(),
        sort=False,
        dropna=False,
    ):
        with_individual_values = subset[
            subset[INDIVIDUAL_RQ_COLUMN].map(lambda value: not is_blank(value))
        ]
        selected_rows = (
            with_individual_values
            if not with_individual_values.empty
            else subset.iloc[:1]
        )
        if (
            not with_individual_values.empty
            and selected_rows[SAMPLE_SIZE_COLUMN].map(is_blank).all()
        ):
            selected_rows = selected_rows.copy()
            selected_rows[SAMPLE_SIZE_COLUMN] = len(with_individual_values)
        groups.append(selected_rows)

    return (
        pd.concat(groups, ignore_index=True)[result_columns()]
        .drop_duplicates(subset=result_columns())
        .reset_index(drop=True)
    )


def summarize_full_table(
    data: pd.DataFrame,
    table_number: int,
    column_mapping: Optional[ColumnMapping] = None,
) -> pd.DataFrame:
    genes = detect_gene_names(data)
    reference_source = "; ".join(genes[1:]) if len(genes) >= 2 else ""
    group_col = require_column(
        f"Group in table {table_number}",
        resolve_mapped_column(
            data,
            column_mapping,
            COLUMN_ROLE_GROUP,
            find_exact_header_column(data, ["group"]),
        ),
    )
    compound_col = resolve_mapped_column(
        data,
        column_mapping,
        COLUMN_ROLE_COMPOUND,
        find_column(data, include_any=["compound id", "compound"]),
    )
    sample_id_col = resolve_mapped_column(
        data,
        column_mapping,
        COLUMN_ROLE_SAMPLE_ID,
        find_column(data, include_any=["sample id"]),
    )
    animal_id_col = resolve_mapped_column(
        data,
        column_mapping,
        COLUMN_ROLE_ANIMAL_ID,
        find_column(data, include_any=["animal id"]),
    )
    individual_rq_col = resolve_mapped_column(
        data,
        column_mapping,
        COLUMN_ROLE_INDIVIDUAL,
        find_individual_rq_column(data),
    )
    sample_size_col = resolve_mapped_column(
        data,
        column_mapping,
        COLUMN_ROLE_SAMPLE_SIZE,
        find_sample_size_column(data),
    )
    mean_rq_col = resolve_mapped_column(
        data,
        column_mapping,
        COLUMN_ROLE_MEAN,
        find_summary_mean_column(data),
    )
    sem_col = resolve_mapped_column(
        data,
        column_mapping,
        COLUMN_ROLE_SEM,
        find_exact_header_column(data, ["sem"]),
    )
    mean_rq_sem_col = find_column(data, include_all=["mean", "rq", "sem"])

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
    summary[SAMPLE_ID_COLUMN] = (
        rows[sample_id_col] if sample_id_col is not None else ""
    )
    summary[ANIMAL_ID_COLUMN] = (
        rows[animal_id_col] if animal_id_col is not None else ""
    )
    summary[INDIVIDUAL_RQ_COLUMN] = (
        rows[individual_rq_col].map(format_number)
        if individual_rq_col is not None
        else ""
    )

    summary = summary[
        summary["Group"].map(lambda value: not is_blank(value))
        & summary[MEAN_RQ_COLUMN].map(lambda value: not is_blank(value))
        & summary[SEM_COLUMN].map(lambda value: not is_blank(value))
    ]
    return finalize_summary(summary)


def summarize_aggregate_table(
    data: pd.DataFrame,
    metadata_table: pd.DataFrame,
    reference_source: str,
    table_number: int,
) -> pd.DataFrame:
    group_col = find_column(metadata_table, include_any=["group"])
    compound_col = find_column(metadata_table, include_any=["compound id", "compound"])
    sample_id_col = find_column(metadata_table, include_any=["sample id"])
    animal_id_col = find_column(metadata_table, include_any=["animal id"])
    individual_rq_col = find_geomean_column(data)
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
    for column in [sample_size_col, mean_rq_sem_col, mean_rq_col, sem_col]:
        if column is not None:
            rows[column] = rows[column].ffill()
    for column in [group_col, compound_col]:
        metadata_rows[column] = metadata_rows[column].ffill()
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
    summary[SAMPLE_ID_COLUMN] = (
        metadata_rows[sample_id_col] if sample_id_col is not None else ""
    )
    summary[ANIMAL_ID_COLUMN] = (
        metadata_rows[animal_id_col] if animal_id_col is not None else ""
    )
    summary[INDIVIDUAL_RQ_COLUMN] = (
        rows[individual_rq_col].map(format_number)
        if individual_rq_col is not None
        else ""
    )

    summary = summary[
        summary["Group"].map(lambda value: not is_blank(value))
        & summary["Compound ID"].map(lambda value: not is_blank(value))
        & summary["Compound ID"].map(lambda value: clean_text(value).upper() not in CONTROL_COMPOUNDS)
        & summary[MEAN_RQ_COLUMN].map(lambda value: not is_blank(value))
        & summary[SEM_COLUMN].map(lambda value: not is_blank(value))
    ]
    return finalize_summary(summary)


def summarize_qpcr_region(
    data: pd.DataFrame,
    column_mapping: Optional[ColumnMapping] = None,
) -> pd.DataFrame:
    tables = split_side_by_side_tables(data)
    if not tables:
        return empty_summary()

    full_tables = [
        (table_number, table)
        for table_number, table in enumerate(tables, start=1)
        if resolve_mapped_column(
            table,
            column_mapping,
            COLUMN_ROLE_GROUP,
            find_exact_header_column(table, ["group"]),
        )
        is not None
    ]
    aggregate_tables = [
        (table_number, table)
        for table_number, table in enumerate(tables, start=1)
        if resolve_mapped_column(
            table,
            column_mapping,
            COLUMN_ROLE_GROUP,
            find_exact_header_column(table, ["group"]),
        )
        is None
        and find_geomean_column(table) is not None
    ]

    summaries = [
        summarize_full_table(table, table_number, column_mapping)
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
    column_mapping: Optional[ColumnMapping] = None,
) -> pd.DataFrame:
    data = read_excel_cells_with_merged_values(input_file, sheet_name)
    summary = summarize_qpcr_region(data, column_mapping)
    if not summary.empty:
        return summary

    embedded_region = normalize_embedded_qpcr_region(data)
    if embedded_region.empty:
        return empty_summary()
    return summarize_qpcr_region(embedded_region, column_mapping)


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
    column_mapping: Optional[ColumnMapping] = None,
) -> pd.DataFrame:
    summary = build_summary(input_file, sheet_name, column_mapping)
    output_sheet = default_plotdata_sheet_name(input_file, sheet_name)
    if append_sheet == "":
        append_sheet = output_sheet
    save_summary(summary, input_file, output_file, append_sheet, output_sheet)
    return summary


def mapping_table_from_cells(
    data: pd.DataFrame,
) -> tuple[pd.DataFrame, int, int]:
    """Return a representative table, original header row, and table-block count."""
    direct_tables = [
        table
        for table in split_side_by_side_tables(data)
        if find_exact_header_column(table, ["group"]) is not None
    ]
    if direct_tables:
        return direct_tables[0], 1, len(direct_tables)

    header_row = find_embedded_qpcr_header_row(data)
    if header_row is None:
        raise ValueError(
            "Could not locate a qPCR header containing Group, Compound ID, Mean, and SEM."
        )
    region = normalize_embedded_qpcr_region(data)
    embedded_tables = [
        table
        for table in split_side_by_side_tables(region)
        if find_exact_header_column(table, ["group"]) is not None
    ]
    if not embedded_tables:
        raise ValueError("Could not identify a table block for column mapping.")
    return embedded_tables[0], int(header_row), len(embedded_tables)


def detected_mapping_columns(table: pd.DataFrame) -> dict[str, Optional[int]]:
    return {
        COLUMN_ROLE_GROUP: find_exact_header_column(table, ["group"]),
        COLUMN_ROLE_COMPOUND: find_column(
            table,
            include_any=["compound id", "compound"],
        ),
        COLUMN_ROLE_SAMPLE_ID: find_column(table, include_any=["sample id"]),
        COLUMN_ROLE_ANIMAL_ID: find_column(table, include_any=["animal id"]),
        COLUMN_ROLE_INDIVIDUAL: find_individual_rq_column(table),
        COLUMN_ROLE_MEAN: find_summary_mean_column(table),
        COLUMN_ROLE_SEM: find_exact_header_column(table, ["sem"]),
        COLUMN_ROLE_SAMPLE_SIZE: find_sample_size_column(table),
    }


def column_mapping_display(table: pd.DataFrame, column: int) -> str:
    from openpyxl.utils import get_column_letter

    top, bottom = column_signature(table, column)
    parts = list(dict.fromkeys(part for part in [top, bottom] if part))
    header = " | ".join(parts) if parts else "Unnamed column"
    return f"{get_column_letter(int(column))}: {header}"


def extraction_preview(summary: pd.DataFrame) -> tuple[int, int, list[str]]:
    bar_count = len(summary.drop_duplicates(subset=summary_columns()))
    individual_count = int(
        summary[INDIVIDUAL_RQ_COLUMN]
        .map(lambda value: not is_blank(value))
        .sum()
    )
    issues = []
    bar_keys = ["Group", "Compound ID", REFERENCE_SOURCE_COLUMN]
    for key, subset in summary.groupby(bar_keys, sort=False, dropna=False):
        values = pd.to_numeric(
            subset[INDIVIDUAL_RQ_COLUMN],
            errors="coerce",
        ).dropna()
        provided = pd.to_numeric(
            subset[MEAN_RQ_COLUMN],
            errors="coerce",
        ).dropna()
        if values.empty or provided.empty:
            continue
        calculated_mean = float(values.mean())
        provided_mean = float(provided.iloc[0])
        tolerance = max(0.01, abs(provided_mean) * 0.02)
        if abs(calculated_mean - provided_mean) > tolerance:
            group, compound, reference = key
            label = f"{clean_text(group)} | {clean_text(compound)}"
            if not is_blank(reference):
                label += f" | {clean_text(reference)}"
            issues.append(
                f"{label}: provided {provided_mean:g}, calculated {calculated_mean:g}"
            )
    return bar_count, individual_count, issues


def choose_column_mapping_gui(root, cells: pd.DataFrame) -> Optional[ColumnMapping]:
    import tkinter as tk
    from tkinter import messagebox, ttk

    table, header_row, table_count = mapping_table_from_cells(cells)
    display_by_column = {
        column: column_mapping_display(table, column)
        for column in table.columns
    }
    signature_by_display = {
        display: column_signature(table, column)
        for column, display in display_by_column.items()
    }
    detected = detected_mapping_columns(table)

    select_column = "— Select a column —"
    not_used = "— Not used —"
    role_rows = [
        (COLUMN_ROLE_GROUP, "Group", True),
        (COLUMN_ROLE_COMPOUND, "Compound ID", True),
        (COLUMN_ROLE_SAMPLE_ID, "Sample ID", False),
        (COLUMN_ROLE_ANIMAL_ID, "Animal ID", False),
        (COLUMN_ROLE_INDIVIDUAL, "Individual plotted value", False),
        (COLUMN_ROLE_MEAN, "Bar mean", True),
        (COLUMN_ROLE_SEM, "Error bar (SEM)", True),
        (COLUMN_ROLE_SAMPLE_SIZE, "Sample size", False),
    ]

    window = tk.Toplevel(root)
    window.title("Review qPCR columns")
    window.resizable(True, False)
    window.columnconfigure(1, weight=1)

    ttk.Label(
        window,
        text="Review the detected columns before extraction.",
        font=("TkDefaultFont", 10, "bold"),
    ).grid(row=0, column=0, columnspan=2, padx=18, pady=(16, 4), sticky="w")
    ttk.Label(
        window,
        text=(
            f"Header rows: {header_row}–{header_row + 1}   |   "
            f"Compatible table blocks: {table_count}"
        ),
    ).grid(row=1, column=0, columnspan=2, padx=18, pady=(0, 4), sticky="w")

    genes = detect_gene_names(table)
    gene_text = "Not detected"
    if genes:
        target = genes[0]
        references = "; ".join(genes[1:]) or "Not detected"
        gene_text = f"Target: {target}   |   Reference: {references}"
    ttk.Label(window, text=gene_text).grid(
        row=2,
        column=0,
        columnspan=2,
        padx=18,
        pady=(0, 12),
        sticky="w",
    )

    variables: dict[str, tk.StringVar] = {}
    column_displays = list(display_by_column.values())
    for row_index, (role, label, required) in enumerate(role_rows, start=3):
        label_text = f"{label} *" if required else label
        ttk.Label(window, text=label_text).grid(
            row=row_index,
            column=0,
            padx=(18, 12),
            pady=4,
            sticky="w",
        )
        detected_column = detected[role]
        default = (
            display_by_column[detected_column]
            if detected_column in display_by_column
            else (select_column if required else not_used)
        )
        variable = tk.StringVar(value=default)
        variables[role] = variable
        values = [select_column] + column_displays
        if not required:
            values.insert(1, not_used)
        ttk.Combobox(
            window,
            textvariable=variable,
            values=values,
            state="readonly",
            width=72,
        ).grid(
            row=row_index,
            column=1,
            padx=(0, 18),
            pady=4,
            sticky="ew",
        )

    ttk.Label(
        window,
        text=(
            "* Required. Two-row headers are shown together so similarly named "
            "columns remain distinct."
        ),
    ).grid(
        row=11,
        column=0,
        columnspan=2,
        padx=18,
        pady=(10, 4),
        sticky="w",
    )
    ttk.Label(
        window,
        text="Selections are matched by header name across compatible table blocks.",
    ).grid(
        row=12,
        column=0,
        columnspan=2,
        padx=18,
        pady=(0, 10),
        sticky="w",
    )

    selected: dict[str, Optional[ColumnMapping]] = {"value": None}

    def use_mapping() -> None:
        missing = [
            label
            for role, label, required in role_rows
            if required and variables[role].get() not in signature_by_display
        ]
        if missing:
            messagebox.showerror(
                "Select required columns",
                f"Choose a column for: {', '.join(missing)}.",
                parent=window,
            )
            return

        mapping: ColumnMapping = {}
        for role, _label, required in role_rows:
            display = variables[role].get()
            mapping[role] = (
                signature_by_display[display]
                if display in signature_by_display
                else None
            )
            if required and mapping[role] is None:
                return
        selected["value"] = mapping
        window.destroy()

    def cancel() -> None:
        selected["value"] = None
        window.destroy()

    buttons = ttk.Frame(window)
    buttons.grid(row=13, column=0, columnspan=2, padx=18, pady=(6, 16), sticky="e")
    ttk.Button(buttons, text="Cancel", command=cancel).grid(row=0, column=0, padx=(0, 8))
    ttk.Button(buttons, text="Continue", command=use_mapping).grid(row=0, column=1)
    window.protocol("WM_DELETE_WINDOW", cancel)
    window.bind("<Return>", lambda _event: use_mapping())
    window.bind("<Escape>", lambda _event: cancel())
    window.grab_set()
    window.wait_window()
    return selected["value"]


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

        cells = read_excel_cells_with_merged_values(input_file, sheet_name)
        column_mapping = choose_column_mapping_gui(root, cells)
        if column_mapping is None:
            return 0

        summary = build_summary(input_file, sheet_name, column_mapping)
        if summary.empty:
            raise ValueError(
                "No rows were extracted. Check the selected column mapping and make "
                "sure the mean, SEM, and individual-value cells contain calculated values."
            )
        bar_count, individual_count, mean_issues = extraction_preview(summary)
        preview_lines = [
            f"Bars: {bar_count}",
            f"Individual values: {individual_count}",
        ]
        if individual_count == 0:
            preview_lines.extend(
                [
                    "",
                    "No individual values were detected; plots will contain bars only.",
                ]
            )
        elif mean_issues:
            preview_lines.extend(
                [
                    "",
                    f"Mean check: {len(mean_issues)} bar(s) differ from the calculated mean.",
                    *mean_issues[:3],
                ]
            )
            if len(mean_issues) > 3:
                preview_lines.append(f"...and {len(mean_issues) - 3} more.")
        else:
            preview_lines.extend(["", "Mean check: all available individual values agree."])
        preview_lines.extend(["", "Continue with extraction?"])
        if not messagebox.askokcancel(
            "Extraction preview",
            "\n".join(preview_lines),
            parent=root,
        ):
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

        save_summary(
            summary,
            input_file,
            append_sheet=append_sheet,
            output_sheet=default_sheet_name,
        )
        messagebox.showinfo(
            "Done",
            (
                f"Extracted {bar_count} bars and {individual_count} individual values."
                f"\n\nSaved to:\n{input_file}\nSheet: {append_sheet}"
            ),
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
