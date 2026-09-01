"""Explicit local transcript-scan workflows.

This module contains no remote-BLAST client.  Public transcript accessions may
be retrieved through the target service, but private guide sequences are used
only by the local scanner.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable

from tools_for_pharma.oligo.ncbi_transport import NcbiHttpClient
from tools_for_pharma.oligo.transcript_scan.models import (
    AntisenseQuery,
    AntisenseRegion,
    ComparisonResult,
    PrivatePanelScanResult,
    QueryTargetSummary,
    TranscriptMatch,
    TranscriptTargetResult,
)
from tools_for_pharma.oligo.transcript_scan.queries import normalize_sequence_type
from tools_for_pharma.oligo.transcript_scan.scanner import (
    DEFAULT_MAX_MISMATCHES,
    closest_transcript_matches,
    comparison_result_for_region,
    scan_antisense_against_transcript,
)
from tools_for_pharma.oligo.transcript_scan.targets import (
    AccessionTargetSource,
    TranscriptTargetSource,
    local_transcript_target,
    read_transcript_input,
    retrieve_transcript_targets,
)


ProgressCallback = Callable[[int, int, str, str], None]
CancelCheck = Callable[[], bool]


@dataclass(frozen=True)
class LocalScanConfig:
    """Inputs needed to prepare transcript targets for a local table scan."""

    target_accessions: tuple[str, ...] = ()
    use_query_target_accession: bool = False
    target_sequence: str | None = None
    target_file: Path | None = None
    email: str | None = None
    tool: str = "tools_for_pharma_oligo"
    cache_dir: Path | None = None
    max_mismatches: int = DEFAULT_MAX_MISMATCHES


@dataclass(frozen=True)
class SingleSequenceScanConfig:
    """One local guide and one explicitly selected transcript source."""

    target_source: TranscriptTargetSource
    email: str | None = None
    tool: str = "tools_for_pharma_oligo"
    request_seconds: float = 15
    max_mismatches: int = DEFAULT_MAX_MISMATCHES
    closest: int | None = None


@dataclass(frozen=True)
class PrivatePanelWorkflowConfig:
    """Retrieval and local-scan settings for a transcript panel."""

    accessions: tuple[str, ...]
    cache_dir: Path
    email: str | None = None
    tool: str = "tools_for_pharma_oligo"
    offline: bool = False
    refresh: bool = False
    request_seconds: float = 15
    max_mismatches: int = DEFAULT_MAX_MISMATCHES
    closest: int | None = None
    download_targets_only: bool = False


def _validate_query_targets(queries: list[AntisenseQuery]) -> None:
    missing_accessions = [query.name for query in queries if not query.target_accession]
    if missing_accessions:
        names = ", ".join(missing_accessions)
        raise ValueError(
            "Missing target accession for query row(s): "
            f"{names}. Fill --target-accession-column for every query."
        )


def _shared_transcript(config: LocalScanConfig) -> tuple[str, str]:
    if len(config.target_accessions) > 1:
        raise ValueError("Multiple --target-accession values require private panel mode.")
    return read_transcript_input(
        transcript_sequence=config.target_sequence,
        transcript_file=config.target_file,
        accession=config.target_accessions[0] if config.target_accessions else None,
        email=config.email,
        tool=config.tool,
        cache_dir=config.cache_dir,
    )


def run_local_scan(
    config: LocalScanConfig,
    queries: list[AntisenseQuery],
    scan_regions: list[AntisenseRegion],
    max_mismatches: int | None = DEFAULT_MAX_MISMATCHES,
) -> list[TranscriptMatch]:
    """Scan one or more private guides against prepared transcript targets."""
    if max_mismatches == DEFAULT_MAX_MISMATCHES:
        max_mismatches = config.max_mismatches
    if config.use_query_target_accession:
        _validate_query_targets(queries)
        shared_transcript = None
    else:
        shared_transcript = _shared_transcript(config)

    matches: list[TranscriptMatch] = []
    for query in queries:
        if config.use_query_target_accession:
            transcript_name, transcript = read_transcript_input(
                accession=query.target_accession,
                email=config.email,
                tool=config.tool,
                cache_dir=config.cache_dir,
            )
        else:
            assert shared_transcript is not None
            transcript_name, transcript = shared_transcript
        for scan_region in scan_regions:
            matches.extend(
                scan_antisense_against_transcript(
                    query.sequence_5to3,
                    transcript,
                    transcript_name=transcript_name,
                    antisense_name=query.name,
                    scan_region=scan_region,
                    max_mismatches=max_mismatches,
                    sequence_type=query.sequence_type,
                )
            )
    return matches


def run_local_scan_with_comparison(
    config: LocalScanConfig,
    queries: list[AntisenseQuery],
    scan_regions: list[AntisenseRegion],
) -> tuple[list[TranscriptMatch], list[ComparisonResult]]:
    """Run a table scan and keep one compact comparison per selected region."""
    if config.use_query_target_accession:
        _validate_query_targets(queries)
        shared_transcript = None
        shared_target_label = ""
    else:
        shared_transcript = _shared_transcript(config)
        shared_target_label = (
            config.target_accessions[0]
            if config.target_accessions
            else shared_transcript[0]
        )

    matches: list[TranscriptMatch] = []
    comparison_results: list[ComparisonResult] = []
    for input_order, query in enumerate(queries, start=1):
        if config.use_query_target_accession:
            transcript_name, transcript = read_transcript_input(
                accession=query.target_accession,
                email=config.email,
                tool=config.tool,
                cache_dir=config.cache_dir,
            )
            target_label = query.target_accession
        else:
            assert shared_transcript is not None
            transcript_name, transcript = shared_transcript
            target_label = shared_target_label

        for scan_region in scan_regions:
            all_region_matches = scan_antisense_against_transcript(
                query.sequence_5to3,
                transcript,
                transcript_name=transcript_name,
                antisense_name=query.name,
                scan_region=scan_region,
                max_mismatches=None,
                sequence_type=query.sequence_type,
            )
            passing_matches = [
                match
                for match in all_region_matches
                if match.mismatches <= config.max_mismatches
            ]
            matches.extend(passing_matches)
            comparison_results.append(
                comparison_result_for_region(
                    input_order=input_order,
                    query=query,
                    target_accession=target_label,
                    scan_region=scan_region,
                    passing_matches=passing_matches,
                    all_matches=all_region_matches,
                )
            )
    return matches, comparison_results


def run_private_panel_scan(
    queries: list[AntisenseQuery],
    targets: list[TranscriptTargetResult],
    scan_regions: list[AntisenseRegion],
    max_mismatches: int,
    closest: int | None = None,
) -> PrivatePanelScanResult:
    """Scan every private guide against every prepared transcript target."""
    matches: list[TranscriptMatch] = []
    panel_closest_matches: list[TranscriptMatch] = []
    summaries: list[QueryTargetSummary] = []
    comparison_results: list[ComparisonResult] = []
    region_names = ";".join(region.name for region in scan_regions)

    for input_order, query in enumerate(queries, start=1):
        sequence_type = normalize_sequence_type(query.sequence_type)
        for target in targets:
            if target.status != "ready":
                for scan_region in scan_regions:
                    comparison_results.append(
                        comparison_result_for_region(
                            input_order=input_order,
                            query=query,
                            target_accession=target.requested_accession,
                            scan_region=scan_region,
                            target_error=target.error,
                        )
                    )
                summaries.append(
                    QueryTargetSummary(
                        query_name=query.name,
                        sequence_type=sequence_type,
                        requested_accession=target.requested_accession,
                        retrieved_accession=target.retrieved_accession,
                        target_status=target.status,
                        scan_status="target_error",
                        scan_regions=region_names,
                        match_count=0,
                        exact_match_count=0,
                        best_mismatches=None,
                        error=target.error,
                    )
                )
                continue

            pair_matches: list[TranscriptMatch] = []
            pair_all_matches: list[TranscriptMatch] = []
            for scan_region in scan_regions:
                region_matches = scan_antisense_against_transcript(
                    query.sequence_5to3,
                    target.sequence_5to3,
                    transcript_name=target.transcript_name or target.retrieved_accession,
                    antisense_name=query.name,
                    scan_region=scan_region,
                    max_mismatches=None,
                    sequence_type=sequence_type,
                )
                pair_all_matches.extend(region_matches)
                passing_region_matches = [
                    match
                    for match in region_matches
                    if match.mismatches <= max_mismatches
                ]
                pair_matches.extend(passing_region_matches)
                comparison_results.append(
                    comparison_result_for_region(
                        input_order=input_order,
                        query=query,
                        target_accession=target.requested_accession,
                        scan_region=scan_region,
                        passing_matches=passing_region_matches,
                        all_matches=region_matches,
                    )
                )
                if closest is not None:
                    panel_closest_matches.extend(
                        closest_transcript_matches(region_matches, closest)
                    )
            matches.extend(pair_matches)
            best_pool = pair_all_matches if closest is not None else pair_matches
            summaries.append(
                QueryTargetSummary(
                    query_name=query.name,
                    sequence_type=sequence_type,
                    requested_accession=target.requested_accession,
                    retrieved_accession=target.retrieved_accession,
                    target_status=target.status,
                    scan_status="matched" if pair_matches else "no_match",
                    scan_regions=region_names,
                    match_count=len(pair_matches),
                    exact_match_count=sum(match.mismatches == 0 for match in pair_matches),
                    best_mismatches=(
                        min(match.mismatches for match in best_pool)
                        if best_pool
                        else None
                    ),
                )
            )

    return PrivatePanelScanResult(
        targets=tuple(targets),
        matches=tuple(matches),
        summaries=tuple(summaries),
        closest_matches=tuple(panel_closest_matches),
        comparison_results=tuple(comparison_results),
    )


def run_single_sequence_scan(
    config: SingleSequenceScanConfig,
    query: AntisenseQuery,
    scan_regions: list[AntisenseRegion],
    *,
    progress_callback: ProgressCallback | None = None,
    client: NcbiHttpClient | None = None,
) -> PrivatePanelScanResult:
    """Prepare one transcript source and run one private guide locally."""
    source = config.target_source
    if isinstance(source, AccessionTargetSource):
        assert source.cache_dir is not None
        targets = retrieve_transcript_targets(
            [source.accession],
            email=config.email,
            tool=config.tool,
            cache_dir=source.cache_dir,
            offline=source.offline,
            refresh=source.refresh,
            request_seconds=config.request_seconds,
            client=client,
            progress_callback=progress_callback,
        )
    else:
        targets = [local_transcript_target(source)]
        if progress_callback:
            progress_callback(1, 1, targets[0].transcript_name, targets[0].cache_status)

    target = targets[0]
    if target.status != "ready":
        target_label = (
            source.accession
            if isinstance(source, AccessionTargetSource)
            else target.transcript_name
        )
        raise ValueError(target.error or f"Transcript {target_label} could not be prepared.")
    return run_private_panel_scan(
        [query],
        targets,
        scan_regions,
        config.max_mismatches,
        closest=config.closest,
    )


def run_private_panel_workflow(
    config: PrivatePanelWorkflowConfig,
    queries: list[AntisenseQuery],
    scan_regions: list[AntisenseRegion],
    *,
    progress_callback: ProgressCallback | None = None,
    cancel_check: CancelCheck | None = None,
) -> PrivatePanelScanResult:
    """Retrieve public transcript references, then scan private guides locally."""
    targets = retrieve_transcript_targets(
        config.accessions,
        email=config.email,
        tool=config.tool,
        cache_dir=config.cache_dir,
        offline=config.offline,
        refresh=config.refresh,
        request_seconds=config.request_seconds,
        progress_callback=progress_callback,
        cancel_check=cancel_check,
    )
    if config.download_targets_only:
        return PrivatePanelScanResult(targets=tuple(targets), matches=(), summaries=())
    return run_private_panel_scan(
        queries,
        targets,
        scan_regions,
        config.max_mismatches,
        closest=config.closest,
    )
