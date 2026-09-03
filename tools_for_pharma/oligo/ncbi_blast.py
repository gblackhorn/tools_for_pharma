"""Compatibility facade for local transcript scanning and remote NCBI BLAST.

Local transcript scans compare AS/SS sequences on this computer. Fetching a
public transcript accession sends the accession and NCBI contact metadata, but
does not send the oligo sequence. Remote BLAST is a separate, explicit CLI
workflow that submits the query sequence to NCBI.

New code should import focused implementations from ``transcript_scan``.
Existing callers may continue to use the stable names re-exported here, and
``python -m tools_for_pharma.oligo.ncbi_blast`` remains supported.
"""

from __future__ import annotations

import sys

from tools_for_pharma.oligo.core import get_complementary_sequence
from tools_for_pharma.oligo.ncbi_transport import (
    BlastSubmission,
    NcbiBlastClient,
    NcbiHttpClient,
    parse_blast_field,
    require_email,
)
from tools_for_pharma.oligo.transcript_accessions import (
    normalize_versioned_refseq_accession,
)
from tools_for_pharma.oligo.transcript_scan.app_services import (
    application_base_dir,
    application_data_dir,
    gui_log_path,
    gui_settings_path,
    load_gui_settings,
    save_gui_settings,
    shared_gui_transcript_cache_dir,
)
from tools_for_pharma.oligo.transcript_scan.cli import (
    args_antisense_queries,
    build_parser,
    local_scan_config_from_args,
    local_transcript_target_from_args,
    panel_accessions_from_args,
    private_panel_cache_dir,
    private_panel_requested,
    read_antisense_queries,
    read_antisense_table,
    read_target_accession_table,
    run_blast_batches,
    run_local_scan,
    run_local_scan_with_comparison,
    run_private_panel_workflow,
    target_accession_values,
    validate_runtime_args,
    write_text,
    main as run_cli,
)
from tools_for_pharma.oligo.transcript_scan.gui import (
    choose_ncbi_gui_mode,
    choose_ncbi_gui_settings,
    choose_sheet_gui,
    choose_single_sequence_gui_settings,
    default_header,
    excel_headers,
    gui_args,
    prompt_and_save_ncbi_email,
    run_gui,
    run_single_sequence_gui,
    run_single_sequence_scan,
    saved_or_prompted_ncbi_email,
    show_single_sequence_result_gui,
    single_sequence_gui_args,
    single_sequence_gui_draft,
)
from tools_for_pharma.oligo.transcript_scan.models import (
    AntisenseQuery,
    AntisenseRegion,
    ComparisonResult,
    PrivatePanelScanResult,
    QueryTargetSummary,
    TranscriptMatch,
    TranscriptTargetResult,
)
from tools_for_pharma.oligo.transcript_scan.queries import (
    antisense_region_sequence,
    assign_unique_blast_query_ids,
    batch_antisense_queries,
    clean_text_for_id,
    default_query_name,
    duplicate_sequence_groups,
    normalize_sequence_type,
    parse_fasta_records,
    parse_plain_antisense_lines,
    parse_scan_region,
    parse_scan_regions,
    read_antisense_file,
    sanitize_fasta_name,
)
from tools_for_pharma.oligo.transcript_scan.reporting import (
    comparison_result_rows,
    default_gui_result_workbook,
    default_private_panel_workbook,
    default_result_workbook,
    filter_blast_rows,
    format_closest_transcript_matches_for_terminal,
    format_single_sequence_scan_result,
    format_transcript_matches_for_terminal,
    input_query_rows,
    metadata_rows,
    transcript_match_rows,
    transcript_matches_to_csv,
    write_excel_workbook,
    write_result_workbook,
)
from tools_for_pharma.oligo.transcript_scan.remote_blast import (
    BlastBatchResult,
    combine_blast_csv,
    multi_fasta,
    normalize_dna,
)
from tools_for_pharma.oligo.transcript_scan.scanner import (
    closest_transcript_matches,
    comparison_result_for_region,
    mismatch_positions,
    scan_antisense_against_transcript,
    scan_sense_against_transcript,
)
from tools_for_pharma.oligo.transcript_scan.targets import (
    AccessionTargetSource,
    LocalFileTargetSource,
    PastedTargetSource,
    TranscriptTargetSource,
    fetch_transcript_fasta,
    format_cached_transcript_fasta,
    prepare_pasted_transcript_sequence,
    read_transcript_input,
    retrieve_transcript_targets,
    transcript_cache_path,
    transcript_target_from_fasta,
    transcript_target_source,
    validate_single_transcript_record,
)
from tools_for_pharma.oligo.transcript_scan.workflows import (
    SingleSequenceScanConfig,
    run_private_panel_scan,
)


__all__ = (
    "AccessionTargetSource", "AntisenseQuery", "AntisenseRegion",
    "BlastBatchResult", "BlastSubmission", "ComparisonResult",
    "LocalFileTargetSource", "NcbiBlastClient", "NcbiHttpClient",
    "PastedTargetSource", "PrivatePanelScanResult", "QueryTargetSummary",
    "SingleSequenceScanConfig", "TranscriptMatch", "TranscriptTargetResult",
    "TranscriptTargetSource", "antisense_region_sequence", "application_base_dir",
    "application_data_dir", "args_antisense_queries",
    "assign_unique_blast_query_ids", "batch_antisense_queries", "build_parser",
    "choose_ncbi_gui_mode", "choose_ncbi_gui_settings", "choose_sheet_gui",
    "choose_single_sequence_gui_settings", "clean_text_for_id",
    "closest_transcript_matches", "combine_blast_csv",
    "comparison_result_for_region", "comparison_result_rows", "default_header",
    "default_gui_result_workbook", "default_private_panel_workbook",
    "default_query_name", "default_result_workbook", "duplicate_sequence_groups",
    "excel_headers", "fetch_transcript_fasta", "filter_blast_rows",
    "format_cached_transcript_fasta",
    "format_closest_transcript_matches_for_terminal",
    "format_single_sequence_scan_result", "format_transcript_matches_for_terminal",
    "get_complementary_sequence", "gui_args", "gui_log_path",
    "gui_settings_path", "input_query_rows", "load_gui_settings",
    "local_scan_config_from_args", "local_transcript_target_from_args", "main",
    "metadata_rows", "mismatch_positions", "multi_fasta", "normalize_dna",
    "normalize_sequence_type", "normalize_versioned_refseq_accession",
    "panel_accessions_from_args", "parse_blast_field", "parse_fasta_records",
    "parse_plain_antisense_lines", "parse_scan_region", "parse_scan_regions",
    "prepare_pasted_transcript_sequence", "private_panel_cache_dir",
    "private_panel_requested", "prompt_and_save_ncbi_email",
    "read_antisense_file", "read_antisense_queries", "read_antisense_table",
    "read_target_accession_table", "read_transcript_input", "require_email",
    "retrieve_transcript_targets", "run_blast_batches", "run_cli", "run_gui",
    "run_local_scan", "run_local_scan_with_comparison", "run_private_panel_scan",
    "run_private_panel_workflow", "run_single_sequence_gui",
    "run_single_sequence_scan", "sanitize_fasta_name", "save_gui_settings",
    "saved_or_prompted_ncbi_email", "scan_antisense_against_transcript",
    "scan_sense_against_transcript", "shared_gui_transcript_cache_dir",
    "show_single_sequence_result_gui", "single_sequence_gui_args",
    "single_sequence_gui_draft", "target_accession_values",
    "transcript_cache_path", "transcript_match_rows", "transcript_matches_to_csv",
    "transcript_target_from_fasta", "transcript_target_source",
    "validate_runtime_args", "validate_single_transcript_record",
    "write_excel_workbook", "write_result_workbook", "write_text",
)


def main() -> int:
    """Run the extracted CLI while preserving the historical GUI launch flag."""
    return run_cli(gui_runner=run_gui)


if __name__ == "__main__":
    raise SystemExit(main())
