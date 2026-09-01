"""Tests for reusable equal-length sequence comparison."""

from __future__ import annotations

import pytest

from tools_for_pharma.sequence.comparison import (
    hamming_distance,
    mismatch_positions_1based,
)


def test_equal_length_comparison_reports_one_based_positions() -> None:
    assert mismatch_positions_1based("AUGC", "AUAC") == (3,)
    assert hamming_distance("AUGC", "AUAC") == 1
    assert hamming_distance("AUGC", "AUGC") == 0


def test_equal_length_comparison_rejects_implicit_alignment() -> None:
    with pytest.raises(ValueError, match="same length; received 4 and 3"):
        mismatch_positions_1based("AUGC", "AUG")
