"""Table input/output helpers for oligo processing."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from tools_for_pharma.oligo.core import (
    DEFAULT_END,
    DEFAULT_START,
    antisense_region_to_sense,
    normalize_rna,
)
from tools_for_pharma.shared.excel_utils import (
    append_or_replace_sheet,
    list_excel_sheets,
    read_table,
    write_table,
)


def default_output_path(input_path: str) -> Path:
    """Return the default destination for processed table results."""
    return Path(input_path)


def output_column_names(start: int = DEFAULT_START, end: int = DEFAULT_END) -> tuple[str, str]:
    """Return clear output names for antisense region and sense complement."""
    return f"antisense_{start}-{end}_5to3", f"sense_{end}-{start}_5to3"


def save_processed_table(
    result,
    source: Path,
    output_path: Optional[str] = None,
    sheet_name: Optional[str] = None,
) -> Path:
    """Save processed results, updating the source table by default."""
    destination = Path(output_path) if output_path else source
    suffix = destination.suffix.lower()

    if suffix == ".csv":
        return write_table(result, str(destination))
    if suffix in {".xlsx", ".xlsm"}:
        if output_path:
            return write_table(result, str(destination), sheet_name=sheet_name)
        target_sheet = sheet_name or list_excel_sheets(destination)[0]
        append_or_replace_sheet(destination, result, target_sheet)
        return destination
    if suffix == ".xls":
        if not output_path:
            raise ValueError(
                "Cannot update .xls files in place. Save as .xlsx first, or pass --output."
            )
        return write_table(result, str(destination))

    return write_table(result, str(destination))


def process_table(
    input_path: str,
    column: str,
    output_path: Optional[str] = None,
    sheet_name: Optional[str] = None,
    start: int = DEFAULT_START,
    end: int = DEFAULT_END,
) -> Path:
    """Process one antisense column and save results to Excel or CSV."""
    source = Path(input_path)
    if not source.exists():
        raise ValueError(f"Input file does not exist: {source}")

    data = read_table(str(source), sheet_name=sheet_name)
    if column not in data.columns:
        available = ", ".join(str(name) for name in data.columns)
        raise ValueError(f"Column '{column}' not found. Available columns: {available}")

    antisense_column, sense_column = output_column_names(start=start, end=end)
    result = data.copy()
    result["normalized_antisense"] = ""
    result[antisense_column] = ""
    result[sense_column] = ""
    result["error"] = ""

    for index, value in result[column].items():
        try:
            normalized = normalize_rna(value)
            antisense_region, sense_5to3 = antisense_region_to_sense(
                normalized,
                start=start,
                end=end,
            )
            result.at[index, "normalized_antisense"] = normalized
            result.at[index, antisense_column] = antisense_region
            result.at[index, sense_column] = sense_5to3
        except ValueError as error:
            result.at[index, "error"] = str(error)

    return save_processed_table(result, source, output_path, sheet_name)
