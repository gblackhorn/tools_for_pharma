"""Plot fitted inhibition curves and mark useful IC values."""

from __future__ import annotations

import argparse
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
import re
import sys
import warnings

import numpy as np
import pandas as pd

from tools_for_pharma.plotting.style import (
    DEFAULT_THEME,
    apply_main_axis_style,
    get_pyplot,
    save_figure_pair,
)
from tools_for_pharma.shared.excel_utils import list_excel_sheets
from tools_for_pharma.shared.text_utils import clean_text, sanitize_filename


INPUT_FILE = Path("inhibition_curve.xlsx")
DEFAULT_IC_MARKERS = [50.0, 75.0, 90.0]
DEFAULT_IC_MARKER_TEXT = "IC50,IC75,IC90"
FIT_METHOD_4PL = "4pl"
FIT_METHOD_INTERPOLATION = "interpolation"
FIT_METHOD_CHOICES = [FIT_METHOD_4PL, FIT_METHOD_INTERPOLATION]
PLOT_MODE_COMBINED = "combined"
PLOT_MODE_SINGLE = "single"
PLOT_MODE_BOTH = "both"
PLOT_MODE_CHOICES = [PLOT_MODE_COMBINED, PLOT_MODE_SINGLE, PLOT_MODE_BOTH]
WIDE_SERIES_COLUMN = "Compound"
WIDE_RESPONSE_COLUMN = "Inhibition rate (%)"
ERROR_COLUMN = "Error"

NUMBER_PATTERN = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
MEAN_ERROR_PATTERN = re.compile(
    rf"^\s*(?P<mean>{NUMBER_PATTERN})\s*(?:\u00b1|\+/-|\+-)\s*"
    rf"(?P<error>{NUMBER_PATTERN})\s*$"
)

FONT_FAMILY = DEFAULT_THEME.font_family
FONT_SIZE = DEFAULT_THEME.body_font_size
TITLE_FONT_SIZE = DEFAULT_THEME.title_font_size
AXIS_LABEL_FONT_SIZE = DEFAULT_THEME.axis_label_font_size
TICK_LABEL_FONT_SIZE = DEFAULT_THEME.tick_font_size
LEGEND_FONT_SIZE = DEFAULT_THEME.legend_font_size
FIGURE_HEIGHT = 10.5
FIGURE_WIDTH = 15.0
SUMMARY_FIGURE_WIDTH = 19.0
GRID_COLOR = DEFAULT_THEME.grid_color
TEXT_COLOR = DEFAULT_THEME.text_color
AXIS_COLOR = DEFAULT_THEME.axis_color
MARKER_COLOR = "#404040"
CURVE_PALETTE = DEFAULT_THEME.palette


@dataclass(frozen=True)
class ICMarker:
    """Interpolated concentration for a target inhibition level."""

    label: str
    target: float
    x: float | None
    y: float | None
    slope: float | None = None


@dataclass(frozen=True)
class CurveDescriptor:
    """Compact numeric summary for an interpolated curve."""

    series_label: str
    fit_method: str
    min_response: float
    max_response: float
    auc: float
    markers: list[ICMarker]
    hill_slope: float | None = None
    r_squared: float | None = None
    fallback_reason: str | None = None


@dataclass(frozen=True)
class CurveFit:
    """Fitted or interpolated curve values and model parameters."""

    method: str
    x_grid: np.ndarray
    y_grid: np.ndarray
    bottom: float | None = None
    top: float | None = None
    log_ic50: float | None = None
    hill_slope: float | None = None
    r_squared: float | None = None
    fallback_reason: str | None = None


def default_plot_dir(input_file: Path) -> Path:
    return input_file.with_name(input_file.stem)


def parse_response_value(value: object) -> tuple[float | None, float | None]:
    text = clean_text(value)
    if not text:
        return None, None

    match = MEAN_ERROR_PATTERN.match(text)
    if match:
        return float(match.group("mean")), float(match.group("error"))

    numeric_value = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if pd.isna(numeric_value):
        return None, None
    return float(numeric_value), None


def response_mean(value: object) -> float | None:
    mean, _error = parse_response_value(value)
    return mean


def row_is_numeric(values: pd.Series) -> bool:
    non_blank = [value for value in values.tolist() if clean_text(value)]
    if not non_blank:
        return False
    return all(response_mean(value) is not None for value in non_blank)


def normalize_raw_table(raw: pd.DataFrame) -> pd.DataFrame:
    raw = raw.dropna(how="all").dropna(axis=1, how="all")
    if raw.empty:
        raise ValueError("Input table is empty.")

    first_row = raw.iloc[0]
    if row_is_numeric(first_row):
        table = raw.copy()
        table.columns = ["Concentration"] + [
            f"Compound {index}"
            for index in range(1, len(table.columns))
        ]
        return table.reset_index(drop=True)

    columns = [
        clean_text(value) or f"Column {index + 1}"
        for index, value in enumerate(first_row)
    ]
    table = raw.iloc[1:].copy()
    table.columns = columns
    return table.reset_index(drop=True)


def read_raw_table(input_file: Path, sheet_name: str | int | None) -> pd.DataFrame:
    suffix = input_file.suffix.lower()
    if suffix in {".xlsx", ".xlsm", ".xls"}:
        raw = pd.read_excel(input_file, sheet_name=sheet_name or 0, header=None)
        return normalize_raw_table(raw)
    if suffix in {".csv", ".txt"}:
        raw = pd.read_csv(input_file, header=None)
        return normalize_raw_table(raw)
    raise ValueError("Input file must be an Excel workbook or CSV/text table.")


def parse_ic_markers(values: str | list[float] | tuple[float, ...] | None) -> list[float]:
    if values is None:
        return DEFAULT_IC_MARKERS.copy()
    if isinstance(values, str):
        markers = [
            float(clean_text(value).upper().removeprefix("IC"))
            for value in values.split(",")
            if clean_text(value)
        ]
    else:
        markers = [float(value) for value in values]
    if not markers:
        raise ValueError("At least one IC marker must be provided.")
    invalid = [marker for marker in markers if marker < 0 or marker > 100]
    if invalid:
        raise ValueError("IC markers must be inhibition percentages from 0 to 100.")
    return markers


def marker_label(target: float) -> str:
    return f"IC{target:g}"


def detect_column(
    table: pd.DataFrame,
    requested_column: str | None,
    label: str,
    name_keywords: tuple[str, ...],
    fallback_numeric_index: int,
) -> str:
    if requested_column:
        if requested_column not in table.columns:
            raise ValueError(f"Input table is missing {label} column: {requested_column}.")
        return requested_column

    for column in table.columns:
        normalized = clean_text(column).lower()
        if any(keyword in normalized for keyword in name_keywords):
            return str(column)

    numeric_columns = []
    for column in table.columns:
        values = table[column].map(response_mean)
        if values.notna().any():
            numeric_columns.append(str(column))

    if len(numeric_columns) > fallback_numeric_index:
        return numeric_columns[fallback_numeric_index]
    raise ValueError(f"Could not detect the {label} column.")


def detect_series_column(
    table: pd.DataFrame,
    requested_column: str | None,
    x_column: str,
    y_column: str,
) -> str | None:
    if requested_column:
        if requested_column not in table.columns:
            raise ValueError(f"Input table is missing series column: {requested_column}.")
        return requested_column

    for column in table.columns:
        normalized = clean_text(column).lower()
        if column in {x_column, y_column}:
            continue
        if any(keyword in normalized for keyword in ("compound", "group", "sample", "series")):
            return str(column)
    return None


def is_wide_curve_table(
    table: pd.DataFrame,
    x_column: str | None = None,
    y_column: str | None = None,
    series_column: str | None = None,
) -> bool:
    if y_column or series_column or len(table.columns) < 2:
        return False
    candidate_x_column = x_column or str(table.columns[0])
    if candidate_x_column not in table.columns:
        raise ValueError(f"Input table is missing concentration/dose column: {candidate_x_column}.")

    x_values = pd.to_numeric(table[candidate_x_column], errors="coerce")
    x_name = clean_text(candidate_x_column).lower()
    first_column_is_x = x_values.notna().any() or any(
        keyword in x_name
        for keyword in ("concentration", "conc", "dose", "amount")
    )
    if not first_column_is_x:
        return False

    response_columns = [
        column
        for column in table.columns
        if column != candidate_x_column
        and table[column].map(response_mean).notna().any()
    ]
    return bool(response_columns)


def prepare_wide_curve_data(
    table: pd.DataFrame,
    x_column: str | None = None,
) -> tuple[pd.DataFrame, str, str, str]:
    x_column = x_column or str(table.columns[0])
    if x_column not in table.columns:
        raise ValueError(f"Input table is missing concentration/dose column: {x_column}.")

    response_columns = [
        column
        for column in table.columns
        if column != x_column
        and table[column].map(response_mean).notna().any()
    ]
    if not response_columns:
        raise ValueError("Wide input must contain at least one numeric inhibition-rate column.")

    wide = table[[x_column, *response_columns]].copy()
    wide[x_column] = pd.to_numeric(wide[x_column], errors="coerce")
    for column in response_columns:
        parsed = wide[column].map(parse_response_value)
        wide[column] = parsed.map(lambda item: item[0])
        wide[f"{column} {ERROR_COLUMN}"] = parsed.map(lambda item: item[1])

    mean_summary = wide.melt(
        id_vars=[x_column],
        value_vars=response_columns,
        var_name=WIDE_SERIES_COLUMN,
        value_name=WIDE_RESPONSE_COLUMN,
    )
    error_summary = wide.melt(
        id_vars=[x_column],
        value_vars=[f"{column} {ERROR_COLUMN}" for column in response_columns],
        var_name=WIDE_SERIES_COLUMN,
        value_name=ERROR_COLUMN,
    )
    error_summary[WIDE_SERIES_COLUMN] = error_summary[WIDE_SERIES_COLUMN].str.removesuffix(
        f" {ERROR_COLUMN}"
    )
    summary = mean_summary.copy()
    summary[ERROR_COLUMN] = error_summary[ERROR_COLUMN]
    summary = summary.dropna(subset=[x_column, WIDE_RESPONSE_COLUMN])
    summary = summary[summary[x_column] > 0]
    summary[WIDE_SERIES_COLUMN] = summary[WIDE_SERIES_COLUMN].map(clean_text)
    if summary.empty:
        raise ValueError("No plottable curve rows were found. Concentrations must be positive numbers.")
    return summary, WIDE_SERIES_COLUMN, x_column, WIDE_RESPONSE_COLUMN


def prepare_curve_data(
    table: pd.DataFrame,
    x_column: str | None = None,
    y_column: str | None = None,
    series_column: str | None = None,
) -> tuple[pd.DataFrame, str | None, str, str]:
    if is_wide_curve_table(table, x_column, y_column, series_column):
        return prepare_wide_curve_data(table, x_column)

    x_column = detect_column(
        table,
        x_column,
        "concentration/dose",
        ("concentration", "conc", "dose", "amount"),
        0,
    )
    y_column = detect_column(
        table,
        y_column,
        "inhibition/response",
        ("inhibition", "inhibit", "response", "rate", "effect"),
        1,
    )
    series_column = detect_series_column(table, series_column, x_column, y_column)

    columns = [x_column, y_column] + ([series_column] if series_column else [])
    summary = table[columns].copy()
    summary[x_column] = pd.to_numeric(summary[x_column], errors="coerce")
    parsed_response = summary[y_column].map(parse_response_value)
    summary[y_column] = parsed_response.map(lambda item: item[0])
    summary[ERROR_COLUMN] = parsed_response.map(lambda item: item[1])
    summary = summary.dropna(subset=[x_column, y_column])
    summary = summary[summary[x_column] > 0]
    if series_column:
        summary[series_column] = summary[series_column].map(clean_text)
        summary = summary[summary[series_column] != ""]
    if summary.empty:
        raise ValueError("No plottable curve rows were found. Concentrations must be positive numbers.")
    return summary, series_column, x_column, y_column


def transform_x(values: pd.Series | np.ndarray, log_x: bool) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    return np.log10(array) if log_x else array


def grouped_curve_points(
    subset: pd.DataFrame,
    x_column: str,
    y_column: str,
) -> pd.DataFrame:
    grouped = subset.sort_values(x_column).groupby(x_column, as_index=False)[y_column].mean()
    if len(grouped) < 2:
        raise ValueError("At least two concentrations are required to draw a curve.")
    return grouped


def interpolate_grid(
    subset: pd.DataFrame,
    x_column: str,
    y_column: str,
    log_x: bool = True,
    points: int = 200,
) -> tuple[np.ndarray, np.ndarray]:
    grouped = grouped_curve_points(subset, x_column, y_column)
    x_values = grouped[x_column].to_numpy(dtype=float)
    y_values = grouped[y_column].to_numpy(dtype=float)
    x_space = transform_x(x_values, log_x)
    grid_space = np.linspace(x_space.min(), x_space.max(), points)
    y_grid = np.interp(grid_space, x_space, y_values)
    x_grid = np.power(10, grid_space) if log_x else grid_space
    return x_grid, y_grid


def four_parameter_logistic(
    x_values: np.ndarray,
    bottom: float,
    top: float,
    log_ic50: float,
    hill_slope: float,
) -> np.ndarray:
    log_x = np.log10(np.asarray(x_values, dtype=float))
    exponent = np.clip((log_ic50 - log_x) * hill_slope, -60, 60)
    return bottom + (top - bottom) / (1 + np.power(10.0, exponent))


def four_pl_sse(params: np.ndarray, x_values: np.ndarray, y_values: np.ndarray) -> float:
    bottom, top, log_ic50, hill_slope = params
    if not np.isfinite(params).all() or hill_slope <= 0:
        return float("inf")
    predicted = four_parameter_logistic(x_values, bottom, top, log_ic50, hill_slope)
    return float(np.sum(np.square(y_values - predicted)))


def estimate_log_ic50(x_values: np.ndarray, y_values: np.ndarray) -> float:
    target = 50.0
    order = np.argsort(x_values)
    x_ordered = x_values[order]
    y_ordered = y_values[order]
    log_x = np.log10(x_ordered)
    for index in range(len(x_ordered) - 1):
        y1 = y_ordered[index]
        y2 = y_ordered[index + 1]
        if y1 == y2:
            continue
        if min(y1, y2) <= target <= max(y1, y2):
            fraction = (target - y1) / (y2 - y1)
            return float(log_x[index] + fraction * (log_x[index + 1] - log_x[index]))
    return float(np.median(log_x))


def optimize_four_pl_candidate(
    initial: np.ndarray,
    x_values: np.ndarray,
    y_values: np.ndarray,
) -> tuple[np.ndarray, float]:
    params = initial.astype(float)
    best_score = four_pl_sse(params, x_values, y_values)
    y_range = max(float(np.nanmax(y_values) - np.nanmin(y_values)), 1.0)
    log_range = max(float(np.log10(np.nanmax(x_values)) - np.log10(np.nanmin(x_values))), 0.1)
    steps = np.array([y_range / 2, y_range / 2, log_range / 3, 0.5], dtype=float)

    for _iteration in range(160):
        improved = False
        for param_index in range(len(params)):
            for direction in (1.0, -1.0):
                candidate = params.copy()
                candidate[param_index] += direction * steps[param_index]
                score = four_pl_sse(candidate, x_values, y_values)
                if score < best_score:
                    params = candidate
                    best_score = score
                    improved = True
        if not improved:
            steps *= 0.62
            if float(np.max(steps)) < 1e-5:
                break
    return params, best_score


def fit_four_parameter_logistic(
    subset: pd.DataFrame,
    x_column: str,
    y_column: str,
    log_x: bool = True,
    points: int = 240,
) -> CurveFit:
    grouped = grouped_curve_points(subset, x_column, y_column)
    if len(grouped) < 4:
        x_grid, y_grid = interpolate_grid(subset, x_column, y_column, log_x, points)
        return CurveFit(
            method=FIT_METHOD_INTERPOLATION,
            x_grid=x_grid,
            y_grid=y_grid,
            fallback_reason="4PL requires at least four concentration levels",
        )

    x_values = grouped[x_column].to_numpy(dtype=float)
    y_values = grouped[y_column].to_numpy(dtype=float)
    y_min = float(np.nanmin(y_values))
    y_max = float(np.nanmax(y_values))
    y_range = max(y_max - y_min, 1.0)
    log_ic50 = estimate_log_ic50(x_values, y_values)
    initial_params = np.array(
        [
            max(0.0, y_min - 0.15 * y_range),
            min(100.0, y_max + 0.15 * y_range),
            log_ic50,
            1.0,
        ],
        dtype=float,
    )

    try:
        from scipy.optimize import curve_fit
        from scipy.optimize import OptimizeWarning

        log_min = float(np.log10(np.nanmin(x_values)))
        log_max = float(np.log10(np.nanmax(x_values)))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", OptimizeWarning)
            params, _covariance = curve_fit(
                four_parameter_logistic,
                x_values,
                y_values,
                p0=initial_params,
                bounds=(
                    [-100.0, -100.0, log_min - 3.0, 0.01],
                    [200.0, 200.0, log_max + 3.0, 10.0],
                ),
                maxfev=20000,
            )
        best_params = np.asarray(params, dtype=float)
        best_score = four_pl_sse(best_params, x_values, y_values)
    except Exception:
        best_params = None
        best_score = float("inf")

    candidate_values = []
    for bottom in (y_min, max(0.0, y_min - 0.15 * y_range), 0.0):
        for top in (y_max, min(100.0, y_max + 0.15 * y_range), 100.0):
            if top <= bottom:
                continue
            for hill in (0.7, 1.0, 1.5, 2.0):
                candidate_values.append(np.array([bottom, top, log_ic50, hill], dtype=float))

    for candidate in candidate_values:
        params, score = optimize_four_pl_candidate(candidate, x_values, y_values)
        if score < best_score:
            best_params = params
            best_score = score

    if best_params is None or not np.isfinite(best_score):
        x_grid, y_grid = interpolate_grid(subset, x_column, y_column, log_x, points)
        return CurveFit(
            method=FIT_METHOD_INTERPOLATION,
            x_grid=x_grid,
            y_grid=y_grid,
            fallback_reason="4PL fit did not converge",
        )

    x_space = transform_x(x_values, log_x)
    grid_space = np.linspace(float(x_space.min()), float(x_space.max()), points)
    x_grid = np.power(10, grid_space) if log_x else grid_space
    y_grid = four_parameter_logistic(x_grid, *best_params)
    predicted = four_parameter_logistic(x_values, *best_params)
    total_sse = float(np.sum(np.square(y_values - np.mean(y_values))))
    r_squared = None if total_sse == 0 else 1 - float(np.sum(np.square(y_values - predicted))) / total_sse
    return CurveFit(
        method=FIT_METHOD_4PL,
        x_grid=x_grid,
        y_grid=y_grid,
        bottom=float(best_params[0]),
        top=float(best_params[1]),
        log_ic50=float(best_params[2]),
        hill_slope=float(best_params[3]),
        r_squared=r_squared,
    )


def fit_curve(
    subset: pd.DataFrame,
    x_column: str,
    y_column: str,
    fit_method: str = FIT_METHOD_4PL,
    log_x: bool = True,
    points: int = 240,
) -> CurveFit:
    if fit_method == FIT_METHOD_INTERPOLATION:
        x_grid, y_grid = interpolate_grid(subset, x_column, y_column, log_x, points)
        return CurveFit(method=FIT_METHOD_INTERPOLATION, x_grid=x_grid, y_grid=y_grid)
    if fit_method == FIT_METHOD_4PL:
        return fit_four_parameter_logistic(subset, x_column, y_column, log_x, points)
    raise ValueError(f"Unknown fit method: {fit_method}")


def interpolate_target(
    subset: pd.DataFrame,
    x_column: str,
    y_column: str,
    target: float,
    log_x: bool = True,
) -> tuple[float | None, float | None]:
    ordered = subset.sort_values(x_column)
    grouped = ordered.groupby(x_column, as_index=False)[y_column].mean()
    x_values = grouped[x_column].to_numpy(dtype=float)
    y_values = grouped[y_column].to_numpy(dtype=float)
    x_space = transform_x(x_values, log_x)

    for index in range(len(grouped) - 1):
        y1 = y_values[index]
        y2 = y_values[index + 1]
        if y1 == y2:
            if target == y1:
                value = x_space[index]
                return float(10**value if log_x else value), 0.0
            continue
        low = min(y1, y2)
        high = max(y1, y2)
        if low <= target <= high:
            fraction = (target - y1) / (y2 - y1)
            value = x_space[index] + fraction * (x_space[index + 1] - x_space[index])
            slope = (y2 - y1) / (x_space[index + 1] - x_space[index])
            return float(10**value if log_x else value), float(slope)
    return None, None


def interpolate_target_x(
    subset: pd.DataFrame,
    x_column: str,
    y_column: str,
    target: float,
    log_x: bool = True,
) -> float | None:
    x_value, _slope = interpolate_target(subset, x_column, y_column, target, log_x)
    return x_value


def calculate_ic_markers(
    subset: pd.DataFrame,
    x_column: str,
    y_column: str,
    targets: list[float] | tuple[float, ...] = DEFAULT_IC_MARKERS,
    log_x: bool = True,
) -> list[ICMarker]:
    markers = []
    for target in targets:
        x_value, slope = interpolate_target(subset, x_column, y_column, float(target), log_x)
        markers.append(
            ICMarker(
                label=marker_label(float(target)),
                target=float(target),
                x=x_value,
                y=float(target) if x_value is not None else None,
                slope=slope,
            )
        )
    return markers


def calculate_fit_ic_markers(
    subset: pd.DataFrame,
    x_column: str,
    y_column: str,
    curve_fit: CurveFit,
    targets: list[float] | tuple[float, ...] = DEFAULT_IC_MARKERS,
    log_x: bool = True,
) -> list[ICMarker]:
    if (
        curve_fit.method != FIT_METHOD_4PL
        or curve_fit.bottom is None
        or curve_fit.top is None
        or curve_fit.log_ic50 is None
        or curve_fit.hill_slope is None
    ):
        return calculate_ic_markers(subset, x_column, y_column, targets, log_x)

    bottom = curve_fit.bottom
    top = curve_fit.top
    response_range = top - bottom
    markers = []
    for target in targets:
        target = float(target)
        if response_range == 0 or not min(bottom, top) <= target <= max(bottom, top):
            markers.append(ICMarker(marker_label(target), target, None, None, None))
            continue

        ratio = response_range / (target - bottom) - 1
        if ratio <= 0:
            markers.append(ICMarker(marker_label(target), target, None, None, None))
            continue

        log_x = curve_fit.log_ic50 - np.log10(ratio) / curve_fit.hill_slope
        x_value = float(10**log_x)
        local_term = np.power(10.0, (curve_fit.log_ic50 - log_x) * curve_fit.hill_slope)
        slope = (
            response_range
            * np.log(10.0)
            * curve_fit.hill_slope
            * local_term
            / np.square(1 + local_term)
        )
        markers.append(
            ICMarker(
                label=marker_label(target),
                target=target,
                x=x_value,
                y=target,
                slope=float(slope),
            )
        )
    return markers


def concentration_label(value: float) -> str:
    return f"{value:.3g}"


def metric_label(value: float) -> str:
    return f"{value:.3g}"


def calculate_curve_descriptor(
    subset: pd.DataFrame,
    series_label: object,
    x_column: str,
    y_column: str,
    targets: list[float] | tuple[float, ...] = DEFAULT_IC_MARKERS,
    log_x: bool = True,
    fit_method: str = FIT_METHOD_INTERPOLATION,
    curve_fit: CurveFit | None = None,
) -> CurveDescriptor:
    curve_fit = curve_fit or fit_curve(subset, x_column, y_column, fit_method, log_x)
    x_space = transform_x(curve_fit.x_grid, log_x)
    return CurveDescriptor(
        series_label=clean_text(series_label),
        fit_method=curve_fit.method,
        min_response=float(subset[y_column].min()),
        max_response=float(subset[y_column].max()),
        auc=float(np.trapezoid(curve_fit.y_grid, x_space)),
        markers=calculate_fit_ic_markers(subset, x_column, y_column, curve_fit, targets, log_x),
        hill_slope=curve_fit.hill_slope,
        r_squared=curve_fit.r_squared,
        fallback_reason=curve_fit.fallback_reason,
    )


def save_plot(axis, output_path: Path) -> list[Path]:
    return save_figure_pair(axis.figure, output_path)


def plot_curve_summary(
    summary: pd.DataFrame,
    output_path: Path,
    title: str,
    x_column: str,
    y_column: str,
    series_column: str | None = None,
    ic_markers: list[float] | tuple[float, ...] = DEFAULT_IC_MARKERS,
    log_x: bool = True,
    fit_method: str = FIT_METHOD_4PL,
    x_label: str | None = None,
    y_label: str | None = None,
    summary_renderer: Callable[[object, list[CurveDescriptor], bool, list[str]], None] | None = None,
) -> list[Path]:
    plt = get_pyplot()
    figure_width = SUMMARY_FIGURE_WIDTH if summary_renderer is not None else FIGURE_WIDTH
    figure, axis = plt.subplots(figsize=(figure_width, FIGURE_HEIGHT))

    if series_column:
        series_items = list(summary.groupby(series_column, sort=False))
    else:
        series_items = [(title, summary)]

    all_marker_rows = []
    descriptors = []
    descriptor_colors = []
    for series_index, (series_label, subset) in enumerate(series_items):
        color = CURVE_PALETTE[series_index % len(CURVE_PALETTE)]
        descriptor_colors.append(color)
        fitted_curve = fit_curve(subset, x_column, y_column, fit_method, log_x)
        label = clean_text(series_label)
        descriptors.append(
            calculate_curve_descriptor(
                subset,
                label,
                x_column,
                y_column,
                ic_markers,
                log_x,
                fit_method,
                fitted_curve,
            )
        )
        axis.plot(fitted_curve.x_grid, fitted_curve.y_grid, color=color, linewidth=1.8, label=label)
        has_error = ERROR_COLUMN in subset.columns and subset[ERROR_COLUMN].notna().any()
        if has_error:
            axis.errorbar(
                subset[x_column],
                subset[y_column],
                yerr=subset[ERROR_COLUMN],
                fmt="o",
                color=color,
                ecolor=color,
                elinewidth=0.8,
                capsize=2.5,
                markersize=4.8,
                markeredgecolor="white",
                markeredgewidth=0.7,
                linewidth=0,
                zorder=3,
            )
        else:
            axis.scatter(
                subset[x_column],
                subset[y_column],
                color=color,
                edgecolor="white",
                linewidth=0.7,
                s=32,
                zorder=3,
            )

        for marker in calculate_fit_ic_markers(
            subset,
            x_column,
            y_column,
            fitted_curve,
            ic_markers,
            log_x,
        ):
            if marker.x is None or marker.y is None:
                continue
            all_marker_rows.append((marker, color, label))
            axis.scatter(
                [marker.x],
                [marker.y],
                color=color,
                edgecolor=MARKER_COLOR,
                linewidth=0.7,
                s=42,
                zorder=4,
            )
            axis.axhline(marker.y, color=color, linewidth=0.6, linestyle=":", alpha=0.45)
            axis.axvline(marker.x, color=color, linewidth=0.6, linestyle=":", alpha=0.45)

    if len(series_items) == 1:
        for marker_index, (marker, color, _series_label) in enumerate(all_marker_rows):
            axis.annotate(
                f"{marker.label}={concentration_label(marker.x)}",
                xy=(marker.x, marker.y),
                xytext=(6, 8 + marker_index % 3 * 10),
                textcoords="offset points",
                color=color,
                fontsize=DEFAULT_THEME.annotation_font_size,
            )

    if log_x:
        axis.set_xscale("log")
    axis.set_title(title, pad=12, color=TEXT_COLOR)
    axis.set_xlabel(x_label or x_column)
    axis.set_ylabel(y_label or y_column)
    axis.yaxis.grid(
        True,
        color=GRID_COLOR,
        linewidth=DEFAULT_THEME.primary_grid_line_width,
    )
    axis.xaxis.grid(
        True,
        color=GRID_COLOR,
        linewidth=DEFAULT_THEME.secondary_grid_line_width,
        alpha=DEFAULT_THEME.secondary_grid_alpha,
    )
    axis.set_axisbelow(True)
    apply_main_axis_style(axis)
    if series_column:
        axis.legend(frameon=False, loc="best")
    if summary_renderer is not None:
        summary_renderer(axis, descriptors, log_x, descriptor_colors)
        axis.figure.tight_layout(rect=(0, 0, 0.69, 1))
    else:
        axis.figure.tight_layout()
    return save_plot(axis, output_path)


def grouped_series_items(
    summary: pd.DataFrame,
    series_column: str | None,
    title: str,
) -> list[tuple[str, pd.DataFrame]]:
    if series_column:
        return [
            (clean_text(series_label), subset.copy())
            for series_label, subset in summary.groupby(series_column, sort=False)
        ]
    return [(title, summary)]


def create_curve_plot(
    input_file: Path,
    output_dir: Path | None = None,
    sheet_name: str | int | None = None,
    x_column: str | None = None,
    y_column: str | None = None,
    series_column: str | None = None,
    title: str | None = None,
    ic_markers: str | list[float] | tuple[float, ...] | None = None,
    log_x: bool = True,
    fit_method: str = FIT_METHOD_4PL,
    plot_mode: str = PLOT_MODE_COMBINED,
    summary_renderer: Callable[[object, list[CurveDescriptor], bool, list[str]], None] | None = None,
    summary_max_curves: int = 3,
) -> list[Path]:
    if plot_mode not in PLOT_MODE_CHOICES:
        raise ValueError(f"Unknown plot mode: {plot_mode}. Choose from {', '.join(PLOT_MODE_CHOICES)}.")

    table = read_raw_table(input_file, sheet_name)
    summary, detected_series_column, detected_x_column, detected_y_column = prepare_curve_data(
        table,
        x_column,
        y_column,
        series_column,
    )
    plot_dir = output_dir if output_dir else default_plot_dir(input_file)
    plot_title = title or input_file.stem
    parsed_ic_markers = parse_ic_markers(ic_markers)
    series_items = grouped_series_items(summary, detected_series_column, plot_title)
    output_paths: list[Path] = []

    should_create_combined = plot_mode in {PLOT_MODE_COMBINED, PLOT_MODE_BOTH}
    should_create_single = plot_mode in {PLOT_MODE_SINGLE, PLOT_MODE_BOTH}

    if should_create_combined or (plot_mode == PLOT_MODE_SINGLE and len(series_items) == 1):
        output_name = sanitize_filename(plot_title) + "_dose_response_curve.png"
        combined_renderer = (
            summary_renderer
            if summary_renderer is not None and len(series_items) <= summary_max_curves
            else None
        )
        output_paths.extend(
            plot_curve_summary(
                summary,
                plot_dir / output_name,
                plot_title,
                detected_x_column,
                detected_y_column,
                detected_series_column,
                parsed_ic_markers,
                log_x,
                fit_method,
                summary_renderer=combined_renderer,
            )
        )

    if should_create_single and len(series_items) > 1:
        for series_label, subset in series_items:
            output_name = (
                f"{sanitize_filename(plot_title)}_{sanitize_filename(series_label)}"
                "_dose_response_curve.png"
            )
            output_paths.extend(
                plot_curve_summary(
                    subset,
                    plot_dir / output_name,
                    f"{plot_title} - {series_label}",
                    detected_x_column,
                    detected_y_column,
                    detected_series_column,
                    parsed_ic_markers,
                    log_x,
                    fit_method,
                    summary_renderer=summary_renderer,
                )
            )

    return output_paths


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
        PLOT_MODE_COMBINED: "Combined overview",
        PLOT_MODE_SINGLE: "One image per compound",
        PLOT_MODE_BOTH: "Combined plus one per compound",
    }
    selected = {"value": PLOT_MODE_COMBINED}
    window = tk.Toplevel(root)
    window.title("Select plot mode")
    window.resizable(False, False)
    window.columnconfigure(1, weight=1)

    ttk.Label(window, text="Plot mode").grid(
        row=0, column=0, padx=16, pady=(16, 8), sticky="w"
    )
    mode_var = tk.StringVar(value=labels[PLOT_MODE_COMBINED])
    mode_box = ttk.Combobox(
        window,
        textvariable=mode_var,
        values=[labels[mode] for mode in PLOT_MODE_CHOICES],
        state="readonly",
        width=34,
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


def run_gui(
    summary_renderer: Callable[[object, list[CurveDescriptor], bool, list[str]], None] | None = None,
) -> int:
    import tkinter as tk
    from tkinter import filedialog, messagebox, simpledialog

    root = tk.Tk()
    root.withdraw()

    try:
        input_path = filedialog.askopenfilename(
            title="Select inhibition curve data file",
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
        if sheet_name is None and is_excel and len(list_excel_sheets(input_file)) > 1:
            return 0

        marker_text = simpledialog.askstring(
            "IC markers",
            "IC values to mark, separated by commas",
            initialvalue=DEFAULT_IC_MARKER_TEXT,
            parent=root,
        )
        if marker_text is None:
            return 0

        plot_mode = choose_plot_mode_gui(root)
        if plot_mode is None:
            return 0

        plot_paths = create_curve_plot(
            input_file,
            sheet_name=sheet_name,
            ic_markers=marker_text,
            plot_mode=plot_mode,
            summary_renderer=summary_renderer,
        )
        message = f"Created {len(plot_paths)} curve plot files."
        if plot_paths:
            message += f"\n\nSaved to:\n{default_plot_dir(input_file)}"
        messagebox.showinfo("Done", message)
        return 0
    except Exception as error:
        messagebox.showerror("Curve plot failed", str(error))
        return 1
    finally:
        root.destroy()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Create a fitted inhibition curve and mark IC values."
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
        "--x-column",
        help="Concentration or dose column. Defaults to the first column for wide tables.",
    )
    parser.add_argument(
        "--y-column",
        help="Inhibition or response column for long-format tables.",
    )
    parser.add_argument(
        "--series-column",
        help="Compound/group column for long-format tables.",
    )
    parser.add_argument("--title", help="Plot title. Defaults to the input filename.")
    parser.add_argument(
        "--ic",
        default=DEFAULT_IC_MARKER_TEXT,
        help=f"Comma-separated IC values to mark. Defaults to {DEFAULT_IC_MARKER_TEXT}.",
    )
    parser.add_argument(
        "--linear-x",
        action="store_true",
        help="Use a linear x-axis instead of the default log concentration axis.",
    )
    parser.add_argument(
        "--fit-method",
        choices=FIT_METHOD_CHOICES,
        default=FIT_METHOD_4PL,
        help="Curve method. Defaults to 4pl, the common biopharma dose-response fit.",
    )
    parser.add_argument(
        "--plot-mode",
        choices=PLOT_MODE_CHOICES,
        default=PLOT_MODE_COMBINED,
        help=(
            "combined = one overview plot; single = one image per compound; "
            "both = create both. Defaults to combined."
        ),
    )
    parser.add_argument(
        "--gui",
        action="store_true",
        help="Choose the data Excel/CSV file with dialogs.",
    )
    return parser


def main(
    summary_renderer: Callable[[object, list[CurveDescriptor], bool, list[str]], None] | None = None,
) -> int:
    args = build_parser().parse_args()
    if args.gui:
        return run_gui(summary_renderer)

    try:
        plot_paths = create_curve_plot(
            args.input,
            args.output_dir,
            args.sheet,
            args.x_column,
            args.y_column,
            args.series_column,
            args.title,
            args.ic,
            not args.linear_x,
            args.fit_method,
            args.plot_mode,
            summary_renderer,
        )
    except Exception as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1

    if not plot_paths:
        print("No curve plots were created.")
        return 0

    print("Curve plot:")
    for plot_path in plot_paths:
        print(f"  {plot_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
