"""Core oligo sequence utilities."""

from __future__ import annotations

from typing import Tuple

from tools_for_pharma.sequence.nucleotides import (
    complement_rna,
    normalize_rna,
    reverse_complement_rna,
    subsequence_1based,
)


DEFAULT_START = 2
DEFAULT_END = 18
def get_subsequence(
    sequence: str,
    start: int = DEFAULT_START,
    end: int = DEFAULT_END,
) -> str:
    """Return a 1-based inclusive subsequence from normalized RNA input."""
    normalized = normalize_rna(sequence)
    return subsequence_1based(normalized, start=start, end=end)


def get_complementary_sequence(sequence: str, reverse: bool = True) -> str:
    """Return the complementary RNA sequence.

    Args:
        sequence: RNA or DNA sequence.
        reverse: Use True to return reverse-complement in 5'->3' orientation.
    """
    if reverse:
        return reverse_complement_rna(sequence)
    return complement_rna(sequence)


def antisense_region_to_sense(
    antisense_5to3: str,
    start: int = DEFAULT_START,
    end: int = DEFAULT_END,
) -> Tuple[str, str]:
    """Return antisense region and complementary sense sequence, both 5'->3'."""
    antisense_region = get_subsequence(antisense_5to3, start=start, end=end)
    sense_5to3 = get_complementary_sequence(antisense_region, reverse=True)
    return antisense_region, sense_5to3


def antisense_2_18_to_sense(antisense_5to3: str) -> Tuple[str, str]:
    """Backward-compatible helper for the default 2-18 antisense region."""
    return antisense_region_to_sense(
        antisense_5to3,
        start=DEFAULT_START,
        end=DEFAULT_END,
    )
