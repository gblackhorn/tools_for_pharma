from __future__ import annotations

from io import StringIO

import pytest

from tools_for_pharma.plotting.bar import get_pyplot as get_bar_pyplot
from tools_for_pharma.plotting.curve import get_pyplot as get_curve_pyplot
from tools_for_pharma.qpcr.plot import get_pyplot as get_qpcr_pyplot
from tools_for_pharma.qpcr.ref_qc_plot import get_pyplot as get_ref_qc_pyplot


@pytest.mark.parametrize(
    "get_pyplot",
    [get_bar_pyplot, get_curve_pyplot, get_qpcr_pyplot, get_ref_qc_pyplot],
)
def test_svg_output_keeps_text_editable(get_pyplot) -> None:
    plt = get_pyplot()
    figure, axis = plt.subplots()
    axis.set_title("Editable title")
    axis.set_xlabel("Dose")

    output = StringIO()
    figure.savefig(output, format="svg")
    plt.close(figure)
    svg = output.getvalue()

    assert "<text" in svg
    assert "Editable title" in svg
    assert "Arial" in svg
