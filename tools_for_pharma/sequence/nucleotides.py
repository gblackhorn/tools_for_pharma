"""Domain-independent nucleotide normalization and sequence operations."""

from __future__ import annotations

import re
from typing import Literal


Alphabet = Literal["DNA", "RNA"]
CleanupMode = Literal["letters", "whitespace"]

_DNA_BASES = frozenset("ACGT")
_RNA_BASES = frozenset("ACGU")
_IUPAC_AMBIGUITY_CODES = frozenset("RYSWKMBDHVN")
_DNA_COMPLEMENT = str.maketrans(
    {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
        "R": "Y",
        "Y": "R",
        "S": "S",
        "W": "W",
        "K": "M",
        "M": "K",
        "B": "V",
        "V": "B",
        "D": "H",
        "H": "D",
        "N": "N",
    }
)
_RNA_COMPLEMENT = str.maketrans(
    {
        "A": "U",
        "U": "A",
        "C": "G",
        "G": "C",
        "R": "Y",
        "Y": "R",
        "S": "S",
        "W": "W",
        "K": "M",
        "M": "K",
        "B": "V",
        "V": "B",
        "D": "H",
        "H": "D",
        "N": "N",
    }
)


class SequenceNormalizationError(ValueError):
    """Structured normalization failure usable by compatibility adapters."""

    def __init__(
        self,
        reason: Literal["empty", "invalid"],
        invalid_bases: tuple[str, ...] = (),
    ) -> None:
        self.reason = reason
        self.invalid_bases = invalid_bases
        if reason == "empty":
            message = "Sequence is empty after cleanup."
        else:
            message = f"Sequence contains invalid bases: {', '.join(invalid_bases)}"
        super().__init__(message)


def _normalized_ambiguity_codes(
    codes: str,
    alphabet: Alphabet,
) -> frozenset[str]:
    normalized = str(codes).upper()
    normalized = (
        normalized.replace("T", "U")
        if alphabet == "RNA"
        else normalized.replace("U", "T")
    )
    allowed = frozenset(normalized)
    unsupported = allowed - _IUPAC_AMBIGUITY_CODES
    if unsupported:
        raise ValueError(
            "Unsupported ambiguity codes: " + ", ".join(sorted(unsupported))
        )
    return allowed


def normalize_nucleotides(
    sequence: object,
    *,
    alphabet: Alphabet,
    allowed_ambiguity_codes: str = "",
    cleanup: CleanupMode = "letters",
) -> str:
    """Normalize and validate a DNA or RNA sequence.

    ``cleanup="letters"`` removes non-letter formatting before validation,
    matching the existing oligo input behavior. ``cleanup="whitespace"``
    removes whitespace only, so punctuation remains visible to validation.
    """
    if alphabet not in {"DNA", "RNA"}:
        raise ValueError("alphabet must be 'DNA' or 'RNA'")
    if cleanup == "letters":
        normalized = re.sub(r"[^A-Za-z]", "", str(sequence)).upper()
    elif cleanup == "whitespace":
        normalized = re.sub(r"\s+", "", str(sequence)).upper()
    else:
        raise ValueError("cleanup must be 'letters' or 'whitespace'")

    if alphabet == "RNA":
        normalized = normalized.replace("T", "U")
        canonical_bases = _RNA_BASES
    else:
        normalized = normalized.replace("U", "T")
        canonical_bases = _DNA_BASES

    if not normalized:
        raise SequenceNormalizationError("empty")

    allowed_bases = canonical_bases | _normalized_ambiguity_codes(
        allowed_ambiguity_codes,
        alphabet,
    )
    invalid_bases = tuple(sorted(set(normalized) - allowed_bases))
    if invalid_bases:
        raise SequenceNormalizationError("invalid", invalid_bases)
    return normalized


def normalize_rna(
    sequence: object,
    *,
    allowed_ambiguity_codes: str = "",
    cleanup: CleanupMode = "letters",
) -> str:
    """Return an uppercase RNA sequence, converting thymine to uracil."""
    return normalize_nucleotides(
        sequence,
        alphabet="RNA",
        allowed_ambiguity_codes=allowed_ambiguity_codes,
        cleanup=cleanup,
    )


def normalize_dna(
    sequence: object,
    *,
    allowed_ambiguity_codes: str = "",
    cleanup: CleanupMode = "letters",
) -> str:
    """Return an uppercase DNA sequence, converting uracil to thymine."""
    return normalize_nucleotides(
        sequence,
        alphabet="DNA",
        allowed_ambiguity_codes=allowed_ambiguity_codes,
        cleanup=cleanup,
    )


def complement_rna(
    sequence: object,
    *,
    allowed_ambiguity_codes: str = "",
) -> str:
    """Return the RNA complement in the input sequence's direction."""
    normalized = normalize_rna(
        sequence,
        allowed_ambiguity_codes=allowed_ambiguity_codes,
    )
    return normalized.translate(_RNA_COMPLEMENT)


def reverse_complement_rna(
    sequence: object,
    *,
    allowed_ambiguity_codes: str = "",
) -> str:
    """Return the RNA reverse complement in 5' to 3' orientation."""
    return complement_rna(
        sequence,
        allowed_ambiguity_codes=allowed_ambiguity_codes,
    )[::-1]


def complement_dna(
    sequence: object,
    *,
    allowed_ambiguity_codes: str = "",
) -> str:
    """Return the DNA complement in the input sequence's direction."""
    normalized = normalize_dna(
        sequence,
        allowed_ambiguity_codes=allowed_ambiguity_codes,
    )
    return normalized.translate(_DNA_COMPLEMENT)


def reverse_complement_dna(
    sequence: object,
    *,
    allowed_ambiguity_codes: str = "",
) -> str:
    """Return the DNA reverse complement in 5' to 3' orientation."""
    return complement_dna(
        sequence,
        allowed_ambiguity_codes=allowed_ambiguity_codes,
    )[::-1]


def subsequence_1based(sequence: str, start: int, end: int) -> str:
    """Return a 1-based inclusive subsequence without changing its alphabet."""
    if start < 1:
        raise ValueError("Start position must be 1 or greater.")
    if end < start:
        raise ValueError("End position must be greater than or equal to start.")
    if len(sequence) < end:
        raise ValueError(f"Sequence must be at least {end} nt long.")
    return sequence[start - 1 : end]
