"""Dependency-light data contracts for local oligo/transcript scanning."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class AntisenseQuery:
    """One named AS or SS input sequence and its source-table fields."""

    name: str
    sequence_5to3: str
    target_accession: str = ""
    target_gene: str = ""
    species: str = ""
    notes: str = ""
    sequence_type: str = "AS"
    blast_query_id: str = field(default="", compare=False)
    source_fields: dict[str, object] = field(
        default_factory=dict,
        compare=False,
        repr=False,
    )


@dataclass(frozen=True)
class AntisenseRegion:
    """A 1-based inclusive AS or SS subregion to scan."""

    name: str
    start: int | None = None
    end: int | None = None


@dataclass(frozen=True)
class TranscriptMatch:
    """One local AS/SS-versus-transcript window match."""

    transcript_name: str
    antisense_name: str
    scan_region: str
    as_region_start: int
    as_region_end: int
    antisense_5to3: str
    antisense_region_5to3: str
    target_5to3: str
    transcript_start: int
    transcript_end: int
    mismatches: int
    transcript_window_5to3: str
    transcript_match_as_5to3: str
    mismatch_positions_1based: tuple[int, ...]
    as_mismatch_positions_1based: tuple[int, ...]
    sequence_type: str = "AS"


@dataclass(frozen=True)
class TranscriptTargetResult:
    """Retrieval and validation status for one transcript target."""

    requested_accession: str
    retrieved_accession: str = ""
    transcript_name: str = ""
    sequence_5to3: str = field(default="", repr=False)
    sequence_length_nt: int = 0
    cache_path: str = ""
    cache_status: str = ""
    exact_version_match: bool = False
    sequence_sha256: str = ""
    retrieved_at_utc: str = ""
    status: str = "error"
    error: str = ""


@dataclass(frozen=True)
class QueryTargetSummary:
    """One status row for a guide-versus-transcript comparison."""

    query_name: str
    sequence_type: str
    requested_accession: str
    retrieved_accession: str
    target_status: str
    scan_status: str
    scan_regions: str
    match_count: int
    exact_match_count: int
    best_mismatches: int | None
    error: str = ""


@dataclass(frozen=True)
class ComparisonResult:
    """One compact best result for a query, transcript, and scan region."""

    input_order: int
    query_name: str
    target_accession: str
    scan_region: str
    region_start: int
    region_end: int
    result: str
    sites_within_threshold: int
    best_mismatches: int | None
    mismatch_positions_in_query_1based: tuple[int, ...] = ()
    best_transcript_start: int | None = None
    best_transcript_end: int | None = None
    query_region_5to3: str = ""
    best_match_in_query_orientation_5to3: str = ""
    differences: str = ""


@dataclass(frozen=True)
class PrivatePanelScanResult:
    """Complete result of a private local transcript-panel scan."""

    targets: tuple[TranscriptTargetResult, ...]
    matches: tuple[TranscriptMatch, ...]
    summaries: tuple[QueryTargetSummary, ...]
    closest_matches: tuple[TranscriptMatch, ...] = ()
    comparison_results: tuple[ComparisonResult, ...] = ()
