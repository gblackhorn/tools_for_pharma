"""Core oligo sequence utilities."""

from __future__ import annotations

import re
from typing import Tuple


DEFAULT_START = 2
DEFAULT_END = 18
_RNA_BASES = {"A", "U", "C", "G"}
_RNA_COMPLEMENT = str.maketrans({"A": "U", "U": "A", "C": "G", "G": "C"})


def normalize_rna(sequence: str) -> str:
    """Normalize a nucleotide sequence to uppercase RNA letters."""
    cleaned = re.sub(r"[^A-Za-z]", "", str(sequence)).upper().replace("T", "U")
    if not cleaned:
        raise ValueError("Sequence is empty after cleanup.")

    invalid = sorted(set(cleaned) - _RNA_BASES)
    if invalid:
        raise ValueError(f"Sequence contains invalid bases: {', '.join(invalid)}")

    return cleaned


def get_subsequence(
    sequence: str,
    start: int = DEFAULT_START,
    end: int = DEFAULT_END,
) -> str:
    """Return a 1-based inclusive subsequence from normalized RNA input."""
    if start < 1:
        raise ValueError("Start position must be 1 or greater.")
    if end < start:
        raise ValueError("End position must be greater than or equal to start.")

    normalized = normalize_rna(sequence)
    if len(normalized) < end:
        raise ValueError(f"Sequence must be at least {end} nt long.")

    return normalized[start - 1 : end]


def get_complementary_sequence(sequence: str, reverse: bool = True) -> str:
    """Return the complementary RNA sequence.

    Args:
        sequence: RNA or DNA sequence.
        reverse: Use True to return reverse-complement in 5'->3' orientation.
    """
    complement = normalize_rna(sequence).translate(_RNA_COMPLEMENT)
    return complement[::-1] if reverse else complement


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
