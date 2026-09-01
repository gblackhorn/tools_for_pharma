"""Domain-independent positional comparison for equal-length sequences."""

from __future__ import annotations


def mismatch_positions_1based(left: str, right: str) -> tuple[int, ...]:
    """Return 1-based positions that differ between equal-length sequences."""
    if len(left) != len(right):
        raise ValueError(
            "Equal-length sequence comparison requires sequences of the same "
            f"length; received {len(left)} and {len(right)}."
        )
    return tuple(
        index
        for index, (left_base, right_base) in enumerate(
            zip(left, right),
            start=1,
        )
        if left_base != right_base
    )


def hamming_distance(left: str, right: str) -> int:
    """Return the number of differing positions in equal-length sequences."""
    return len(mismatch_positions_1based(left, right))
