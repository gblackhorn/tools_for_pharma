"""IC50-specific summary panel rendering for fitted dose-response plots."""

from __future__ import annotations

from tools_for_pharma.plotting.curve import (
    CurveDescriptor,
    FIT_METHOD_4PL,
    concentration_label,
    metric_label,
)
from tools_for_pharma.plotting.style import DEFAULT_THEME, PlotTheme


def descriptor_lines(descriptors: list[CurveDescriptor], log_x: bool) -> list[str]:
    auc_unit = "AUC(log dose)" if log_x else "AUC"
    lines = ["Curve summary"]
    for descriptor in descriptors:
        slope_unit = (
            "per log10 dose"
            if descriptor.fit_method == FIT_METHOD_4PL or log_x
            else "per dose"
        )
        lines.append(descriptor.series_label)
        if descriptor.fit_method == FIT_METHOD_4PL:
            fit_parts = ["  fit 4PL"]
            if descriptor.hill_slope is not None:
                fit_parts.append(f"Hill {metric_label(descriptor.hill_slope)}")
            if descriptor.r_squared is not None:
                fit_parts.append(f"R2 {metric_label(descriptor.r_squared)}")
            lines.append("; ".join(fit_parts))
        elif descriptor.fallback_reason:
            lines.append(f"  fit interpolation; {descriptor.fallback_reason}")
        else:
            lines.append("  fit interpolation")
        lines.append(
            f"  range {metric_label(descriptor.min_response)}-"
            f"{metric_label(descriptor.max_response)}; "
            f"{auc_unit} {metric_label(descriptor.auc)}"
        )
        for marker in descriptor.markers:
            if marker.x is None:
                lines.append(f"  {marker.label}: outside range")
                continue
            slope = (
                "NA"
                if marker.slope is None
                else f"{metric_label(marker.slope)} {slope_unit}"
            )
            lines.append(f"  {marker.label}: {concentration_label(marker.x)}, slope {slope}")
    return lines


def summary_rows(descriptor: CurveDescriptor, log_x: bool) -> list[tuple[str, str, bool]]:
    slope_unit = (
        "/log10 dose"
        if descriptor.fit_method == FIT_METHOD_4PL or log_x
        else "/dose"
    )
    rows: list[tuple[str, str, bool]] = []
    if descriptor.fit_method == FIT_METHOD_4PL:
        rows.append(
            (
                "Hill",
                metric_label(descriptor.hill_slope) if descriptor.hill_slope is not None else "NA",
                True,
            )
        )
        rows.append(
            (
                "R2",
                metric_label(descriptor.r_squared) if descriptor.r_squared is not None else "NA",
                False,
            )
        )
    elif descriptor.fallback_reason:
        rows.append(("Fit", f"interpolation ({descriptor.fallback_reason})", False))
    else:
        rows.append(("Fit", "interpolation", False))

    for marker in descriptor.markers:
        if marker.x is None:
            rows.append((marker.label, "outside range", True))
        else:
            slope = "NA" if marker.slope is None else f"{metric_label(marker.slope)} {slope_unit}"
            rows.append((marker.label, concentration_label(marker.x), True))
            rows.append(("slope", slope, False))

    rows.append(
        (
            "Range",
            f"{metric_label(descriptor.min_response)}-{metric_label(descriptor.max_response)}%",
            False,
        )
    )
    rows.append(("AUC", metric_label(descriptor.auc), False))
    return rows


def add_summary_panel(
    axis,
    descriptors: list[CurveDescriptor],
    log_x: bool,
    descriptor_colors: list[str],
    *,
    theme: PlotTheme = DEFAULT_THEME,
) -> None:
    if not descriptors:
        return
    axis.figure.subplots_adjust(right=0.69)
    panel = axis.inset_axes(
        [1.05, 0.04, 0.58, 0.9],
        transform=axis.transAxes,
    )
    panel.set_facecolor("white")
    panel.set_xticks([])
    panel.set_yticks([])
    for spine in panel.spines.values():
        spine.set_visible(True)
        spine.set_edgecolor(theme.grid_color)
        spine.set_linewidth(theme.axis_line_width)
    panel.set_xlim(0, 1)
    panel.set_ylim(0, 1)

    panel.text(
        0.05,
        0.96,
        "Curve summary",
        ha="left",
        va="top",
        color=theme.text_color,
        fontsize=theme.panel_title_font_size,
        fontweight="bold",
    )
    y = 0.88
    for index, descriptor in enumerate(descriptors):
        color = descriptor_colors[index % len(descriptor_colors)]
        panel.text(
            0.05,
            y,
            descriptor.series_label,
            ha="left",
            va="top",
            color=color,
            fontsize=theme.panel_series_font_size,
            fontweight="bold",
        )
        y -= 0.048
        for metric, value, is_key in summary_rows(descriptor, log_x):
            if y < 0.05:
                panel.text(
                    0.05,
                    y,
                    "...",
                    ha="left",
                    va="top",
                    color=theme.text_color,
                    fontsize=theme.panel_compact_font_size,
                )
                return
            row_weight = "bold" if is_key else "normal"
            is_slope_row = metric == "slope"
            metric_x = 0.10 if is_slope_row else 0.05
            value_x = 0.42
            panel.text(
                metric_x,
                y,
                metric,
                ha="left",
                va="top",
                color=theme.text_color,
                fontsize=(
                    theme.panel_compact_font_size
                    if is_slope_row
                    else theme.panel_metric_font_size
                ),
                fontweight=row_weight,
                family="monospace",
            )
            panel.text(
                value_x,
                y,
                value,
                ha="left",
                va="top",
                color=theme.text_color,
                fontsize=(
                    theme.panel_compact_font_size
                    if is_slope_row
                    else theme.panel_metric_font_size
                ),
                fontweight=row_weight,
            )
            y -= 0.031 if is_slope_row else 0.035
        y -= 0.048
