"""Shared visual theme and reusable styling helpers for Matplotlib plots."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class PlotTheme:
    """Visual defaults shared by the repository's scientific plots."""

    font_family: str = "Arial"
    font_fallbacks: tuple[str, ...] = ("DejaVu Sans", "Calibri")
    body_font_size: float = 20
    title_font_size: float = 30
    axis_label_font_size: float = 28
    tick_font_size: float = 22
    legend_font_size: float = 22
    annotation_font_size: float = 18
    panel_title_font_size: float = 16
    panel_series_font_size: float = 14
    panel_metric_font_size: float = 11.5
    panel_compact_font_size: float = 11

    text_color: str = "#222222"
    axis_color: str = "#444444"
    grid_color: str = "#E4E4E4"
    figure_facecolor: str = "white"
    axes_facecolor: str = "white"
    palette: tuple[str, ...] = (
        "#4E79A7",
        "#59A14F",
        "#F28E2B",
        "#B07AA1",
        "#E15759",
        "#76B7B2",
    )

    axis_line_width: float = 0.8
    tick_length: float = 3
    tick_width: float = 0.8
    spine_offset_points: float = 7
    primary_grid_line_width: float = 0.55
    secondary_grid_line_width: float = 0.45
    secondary_grid_alpha: float = 0.6

    png_dpi: int = 300
    svg_fonttype: str = "none"

    def rc_params(self) -> dict[str, object]:
        """Return Matplotlib rcParams for typography, colors, and export behavior."""

        return {
            "font.family": "sans-serif",
            "font.sans-serif": [self.font_family, *self.font_fallbacks],
            "font.size": self.body_font_size,
            "axes.titlesize": self.title_font_size,
            "axes.labelsize": self.axis_label_font_size,
            "axes.titleweight": "bold",
            "axes.edgecolor": self.axis_color,
            "axes.labelcolor": self.text_color,
            "xtick.color": self.text_color,
            "ytick.color": self.text_color,
            "xtick.labelsize": self.tick_font_size,
            "ytick.labelsize": self.tick_font_size,
            "legend.fontsize": self.legend_font_size,
            "figure.facecolor": self.figure_facecolor,
            "axes.facecolor": self.axes_facecolor,
            "savefig.facecolor": self.figure_facecolor,
            "svg.fonttype": self.svg_fonttype,
        }


DEFAULT_THEME = PlotTheme()


def get_pyplot(theme: PlotTheme = DEFAULT_THEME):
    """Return pyplot configured with the selected repository theme."""

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(theme.rc_params())
    return plt


def apply_main_axis_style(
    axis,
    theme: PlotTheme = DEFAULT_THEME,
    *,
    detach_spines: bool = True,
) -> None:
    """Apply the shared open-axis treatment to a plot's main data axes."""

    axis.tick_params(
        axis="both",
        labelsize=theme.tick_font_size,
        length=theme.tick_length,
        width=theme.tick_width,
    )
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_linewidth(theme.axis_line_width)
    axis.spines["bottom"].set_linewidth(theme.axis_line_width)
    offset = theme.spine_offset_points if detach_spines else 0
    axis.spines["left"].set_position(("outward", offset))
    axis.spines["bottom"].set_position(("outward", offset))


def save_figure_pair(
    figure,
    output_path: Path,
    theme: PlotTheme = DEFAULT_THEME,
) -> list[Path]:
    """Save a figure as PNG and editable-text SVG, then close it."""

    plt = get_pyplot(theme)
    png_path = output_path.with_suffix(".png")
    svg_path = output_path.with_suffix(".svg")
    png_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(png_path, dpi=theme.png_dpi, bbox_inches="tight")
    figure.savefig(svg_path, bbox_inches="tight")
    plt.close(figure)
    return [png_path, svg_path]
