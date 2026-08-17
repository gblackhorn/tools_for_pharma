from __future__ import annotations

from dataclasses import replace

from tools_for_pharma.plotting.curve import CurveDescriptor, FIT_METHOD_4PL
from tools_for_pharma.plotting.ic50_summary import add_summary_panel
from tools_for_pharma.plotting.style import (
    DEFAULT_THEME,
    apply_main_axis_style,
    get_pyplot,
)


def test_default_theme_matches_large_readable_plot_typography() -> None:
    assert DEFAULT_THEME.body_font_size == 20
    assert DEFAULT_THEME.title_font_size == 30
    assert DEFAULT_THEME.axis_label_font_size == 28
    assert DEFAULT_THEME.tick_font_size == 22
    assert DEFAULT_THEME.legend_font_size == 22


def test_main_axis_style_supports_local_overrides() -> None:
    local_theme = replace(
        DEFAULT_THEME,
        tick_font_size=17,
        spine_offset_points=11,
    )
    plt = get_pyplot(local_theme)
    figure, axis = plt.subplots()

    apply_main_axis_style(axis, local_theme)

    assert not axis.spines["top"].get_visible()
    assert not axis.spines["right"].get_visible()
    assert axis.spines["left"].get_position() == ("outward", 11)
    assert axis.spines["bottom"].get_position() == ("outward", 11)
    assert {label.get_fontsize() for label in axis.get_xticklabels()} == {17}
    plt.close(figure)


def test_main_axis_style_can_leave_spines_attached() -> None:
    plt = get_pyplot()
    figure, axis = plt.subplots()

    apply_main_axis_style(axis, detach_spines=False)

    assert axis.spines["left"].get_position() == ("outward", 0)
    assert axis.spines["bottom"].get_position() == ("outward", 0)
    plt.close(figure)


def test_ic50_panel_uses_separate_overrideable_typography() -> None:
    panel_theme = replace(
        DEFAULT_THEME,
        panel_title_font_size=15,
        panel_series_font_size=13,
        panel_metric_font_size=10.5,
        panel_compact_font_size=10,
    )
    descriptor = CurveDescriptor(
        series_label="Compound A",
        fit_method=FIT_METHOD_4PL,
        min_response=5,
        max_response=95,
        auc=120,
        markers=[],
        hill_slope=1.2,
        r_squared=0.99,
    )
    plt = get_pyplot(panel_theme)
    figure, axis = plt.subplots()

    add_summary_panel(
        axis,
        [descriptor],
        log_x=True,
        descriptor_colors=["#4E79A7"],
        theme=panel_theme,
    )

    panel = axis.child_axes[0]
    text_by_value = {text.get_text(): text for text in panel.texts}
    assert text_by_value["Curve summary"].get_fontsize() == 15
    assert text_by_value["Compound A"].get_fontsize() == 13
    assert panel_theme.panel_metric_font_size < panel_theme.tick_font_size
    plt.close(figure)
