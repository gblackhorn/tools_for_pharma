"""Pure AS/SS transcript-window scanning and best-match selection."""

from __future__ import annotations

from typing import Iterable

from tools_for_pharma.oligo.transcript_scan.models import (
    AntisenseQuery,
    AntisenseRegion,
    ComparisonResult,
    TranscriptMatch,
)
from tools_for_pharma.oligo.transcript_scan.queries import (
    antisense_region_sequence,
    normalize_sequence_type,
)
from tools_for_pharma.sequence.comparison import mismatch_positions_1based
from tools_for_pharma.sequence.nucleotides import (
    normalize_rna,
    reverse_complement_rna,
)


DEFAULT_MAX_MISMATCHES = 3


def mismatch_positions(query: str, target: str) -> tuple[int, ...]:
    """Return 1-based mismatch positions between equal-length RNA strings."""
    return mismatch_positions_1based(query, target)


def scan_antisense_against_transcript(
    antisense_5to3: str,
    transcript_sequence: str,
    transcript_name: str = "target_transcript",
    antisense_name: str = "antisense_query",
    scan_region: AntisenseRegion | None = None,
    max_mismatches: int | None = DEFAULT_MAX_MISMATCHES,
    sequence_type: str = "AS",
) -> list[TranscriptMatch]:
    """Find AS reverse-complement or SS direct windows in a transcript."""
    normalized_type = normalize_sequence_type(sequence_type)
    antisense = normalize_rna(antisense_5to3)
    transcript = normalize_rna(transcript_sequence)
    region = scan_region or AntisenseRegion("full")
    region_sequence, region_start, region_end = antisense_region_sequence(
        antisense,
        region,
    )
    target = (
        reverse_complement_rna(region_sequence)
        if normalized_type == "AS"
        else region_sequence
    )
    if len(transcript) < len(target):
        return []

    matches = []
    for start_index in range(0, len(transcript) - len(target) + 1):
        window = transcript[start_index : start_index + len(target)]
        mismatches = mismatch_positions(target, window)
        if max_mismatches is None or len(mismatches) <= max_mismatches:
            transcript_match_as = (
                reverse_complement_rna(window)
                if normalized_type == "AS"
                else window
            )
            matches.append(
                TranscriptMatch(
                    transcript_name=transcript_name,
                    antisense_name=antisense_name,
                    scan_region=region.name,
                    as_region_start=region_start,
                    as_region_end=region_end,
                    antisense_5to3=antisense,
                    antisense_region_5to3=region_sequence,
                    target_5to3=target,
                    transcript_start=start_index + 1,
                    transcript_end=start_index + len(target),
                    mismatches=len(mismatches),
                    transcript_window_5to3=window,
                    transcript_match_as_5to3=transcript_match_as,
                    mismatch_positions_1based=mismatches,
                    as_mismatch_positions_1based=mismatch_positions(
                        region_sequence,
                        transcript_match_as,
                    ),
                    sequence_type=normalized_type,
                )
            )
    return sorted(
        matches,
        key=lambda item: (item.mismatches, item.transcript_start),
    )


def scan_sense_against_transcript(
    sense_5to3: str,
    transcript_sequence: str,
    transcript_name: str = "target_transcript",
    sense_name: str = "ss_query",
    scan_region: AntisenseRegion | None = None,
    max_mismatches: int | None = DEFAULT_MAX_MISMATCHES,
) -> list[TranscriptMatch]:
    """Find direct sense target windows in a transcript."""
    return scan_antisense_against_transcript(
        antisense_5to3=sense_5to3,
        transcript_sequence=transcript_sequence,
        transcript_name=transcript_name,
        antisense_name=sense_name,
        scan_region=scan_region,
        max_mismatches=max_mismatches,
        sequence_type="SS",
    )


def closest_transcript_matches(
    matches: Iterable[TranscriptMatch],
    limit: int,
) -> list[TranscriptMatch]:
    """Return the closest windows with stable positional tie breaking."""
    if limit < 1:
        raise ValueError("--closest must be 1 or greater.")
    return sorted(
        matches,
        key=lambda match: (
            match.mismatches,
            match.transcript_start,
            match.antisense_name,
            match.scan_region,
        ),
    )[:limit]


def comparison_result_for_region(
    *,
    input_order: int,
    query: AntisenseQuery,
    target_accession: str,
    scan_region: AntisenseRegion,
    passing_matches: Iterable[TranscriptMatch] = (),
    all_matches: Iterable[TranscriptMatch] = (),
    target_error: str = "",
) -> ComparisonResult:
    """Build one user-facing best comparison for a query/target/region."""
    query_region, region_start, region_end = antisense_region_sequence(
        query.sequence_5to3,
        scan_region,
    )
    passing = sorted(
        passing_matches,
        key=lambda match: (match.mismatches, match.transcript_start),
    )
    candidates = sorted(
        all_matches,
        key=lambda match: (match.mismatches, match.transcript_start),
    )
    best = passing[0] if passing else candidates[0] if candidates else None

    if target_error:
        result = "target_error"
    elif passing:
        result = "exact_match" if best and best.mismatches == 0 else "match"
    else:
        result = "no_match"

    mismatch_positions_in_query: tuple[int, ...] = ()
    differences = ""
    if best is not None:
        mismatch_positions_in_query = tuple(
            region_start + position - 1
            for position in best.as_mismatch_positions_1based
        )
        difference_items = []
        for relative_position, (expected_base, observed_base) in enumerate(
            zip(query_region, best.transcript_match_as_5to3),
            start=1,
        ):
            if expected_base != observed_base:
                full_position = region_start + relative_position - 1
                difference_items.append(
                    f"{full_position}:{expected_base}>{observed_base}"
                )
        differences = "; ".join(difference_items) or "None"

    return ComparisonResult(
        input_order=input_order,
        query_name=query.name,
        target_accession=target_accession,
        scan_region=scan_region.name,
        region_start=region_start,
        region_end=region_end,
        result=result,
        sites_within_threshold=len(passing),
        best_mismatches=best.mismatches if best is not None else None,
        mismatch_positions_in_query_1based=mismatch_positions_in_query,
        best_transcript_start=(
            best.transcript_start if best is not None else None
        ),
        best_transcript_end=best.transcript_end if best is not None else None,
        query_region_5to3=query_region,
        best_match_in_query_orientation_5to3=(
            best.transcript_match_as_5to3 if best is not None else ""
        ),
        differences=differences,
    )
