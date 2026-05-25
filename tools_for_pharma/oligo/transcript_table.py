"""Excel/CSV table helpers for transcript coordinate matching."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd

from tools_for_pharma.oligo.table import read_table, write_table
from tools_for_pharma.oligo.transcript import (
    fasta_or_plain_text_to_sequence,
    read_text_sequence_file,
    transcript_region_to_oligo,
)


SS_MATCHED_COLUMN = "SS matched"
AS_MATCHED_COLUMN = "AS matched"
MATCHED_LENGTH_COLUMN = "matched length (nt)"
MATCH_ERROR_COLUMN = "match error"


def default_transcript_range_output_path(input_path: str) -> Path:
    """Return the default output path for a transcript coordinate table."""
    source = Path(input_path)
    return source.with_name(f"{source.stem}_oligo_matched.xlsx")


def find_column_case_insensitive(data: pd.DataFrame, column_name: str) -> str:
    """Return the actual DataFrame column matching a name case-insensitively."""
    target = column_name.strip().lower()
    for column in data.columns:
        if str(column).strip().lower() == target:
            return str(column)

    available = ", ".join(str(column) for column in data.columns)
    raise ValueError(f"Column '{column_name}' not found. Available columns: {available}")


def parse_position(value: object, column_name: str) -> int:
    """Parse a 1-based coordinate from an Excel/CSV cell."""
    text = str(value).strip().replace(",", "")
    if not text:
        raise ValueError(f"{column_name} is empty")

    try:
        number = float(text)
    except ValueError as error:
        raise ValueError(f"{column_name} is not a number: {value}") from error

    if not number.is_integer():
        raise ValueError(f"{column_name} must be a whole number: {value}")

    position = int(number)
    if position < 1:
        raise ValueError(f"{column_name} must be 1 or greater: {value}")
    return position


def ensure_result_columns(data: pd.DataFrame) -> pd.DataFrame:
    """Ensure matched-output columns exist, appending missing ones at the end."""
    result = data.copy()
    for column in [
        SS_MATCHED_COLUMN,
        AS_MATCHED_COLUMN,
        MATCHED_LENGTH_COLUMN,
        MATCH_ERROR_COLUMN,
    ]:
        if column not in result.columns:
            result[column] = ""
    result[MATCHED_LENGTH_COLUMN] = result[MATCHED_LENGTH_COLUMN].astype(object)
    return result


def is_empty_cell(value: object) -> bool:
    """Return True when a coordinate cell should be treated as blank."""
    return pd.isna(value) or str(value).strip() == ""


def process_transcript_range_table(
    transcript_path: str,
    table_path: str,
    output_path: Optional[str] = None,
    sheet_name: Optional[str] = None,
) -> Path:
    """Fill SS/AS matched columns using start/end coordinates from a table."""
    transcript_text = read_text_sequence_file(transcript_path)
    transcript_sequence = fasta_or_plain_text_to_sequence(transcript_text)

    data = read_table(table_path, sheet_name=sheet_name)
    start_column = find_column_case_insensitive(data, "start")
    end_column = find_column_case_insensitive(data, "end")
    result = ensure_result_columns(data)

    for index, row in result.iterrows():
        if is_empty_cell(row[start_column]) or is_empty_cell(row[end_column]):
            result.at[index, SS_MATCHED_COLUMN] = ""
            result.at[index, AS_MATCHED_COLUMN] = ""
            result.at[index, MATCHED_LENGTH_COLUMN] = ""
            result.at[index, MATCH_ERROR_COLUMN] = ""
            continue

        try:
            start = parse_position(row[start_column], "start")
            end = parse_position(row[end_column], "end")
            matched = transcript_region_to_oligo(
                transcript_sequence,
                start=start,
                end=end,
            )
            result.at[index, SS_MATCHED_COLUMN] = matched.ss_5to3
            result.at[index, AS_MATCHED_COLUMN] = matched.as_5to3
            result.at[index, MATCHED_LENGTH_COLUMN] = matched.selected_length
            result.at[index, MATCH_ERROR_COLUMN] = ""
        except ValueError as error:
            result.at[index, SS_MATCHED_COLUMN] = ""
            result.at[index, AS_MATCHED_COLUMN] = ""
            result.at[index, MATCHED_LENGTH_COLUMN] = ""
            result.at[index, MATCH_ERROR_COLUMN] = str(error)

    destination = (
        Path(output_path)
        if output_path
        else default_transcript_range_output_path(table_path)
    )
    return write_table(result, str(destination))
