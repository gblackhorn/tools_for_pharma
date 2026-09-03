"""Tkinter interface for private local transcript scanning.

This module has two related workflows:

1. Specific transcript check:
   Fetch an NM/XM/NR/XR accession with NCBI EFetch, then scan the transcript for
   the reverse-complement target of an antisense sequence or the direct target
   of a sense sequence.

2. BLAST database search:
   Submit the oligo sequence to the NCBI BLAST URL API and retrieve a CSV
   report. This is best for broad searches such as refseq_rna or core_nt.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import logging
import math
import os
from pathlib import Path
import re
from typing import Callable, Iterable

from tools_for_pharma.oligo.core import get_complementary_sequence
from tools_for_pharma.oligo.ncbi_transport import (
    BLAST_URL,
    BLASTN_WORD_SIZES,
    DEFAULT_DATABASE,
    DEFAULT_EMAIL,
    DEFAULT_EXPECT,
    DEFAULT_HITLIST_SIZE,
    DEFAULT_MEGABLAST_WORD_SIZE,
    DEFAULT_POLL_SECONDS,
    DEFAULT_PROGRAM,
    DEFAULT_REQUEST_SECONDS,
    DEFAULT_TOOL,
    DEFAULT_WORD_SIZE,
    EFETCH_URL,
    MEGABLAST_WORD_SIZES,
    BlastSubmission,
    NcbiBlastClient,
    NcbiHttpClient,
    efetch_fasta_params,
    parse_blast_field,
    require_email,
    resolve_blast_word_size,
)
from tools_for_pharma.oligo.transcript import fasta_or_plain_text_to_sequence, get_fasta_header
from tools_for_pharma.oligo.transcript_accessions import (
    VERSIONED_REFSEQ_GENOMIC_RE,
    VERSIONED_REFSEQ_TRANSCRIPT_RE,
    extract_refseq_accession_from_header,
    normalize_versioned_refseq_accession,
)
from tools_for_pharma.oligo.transcript_scan.app_services import (
    APP_DATA_DIR_NAME,
    GUI_LOG_FILE_NAME,
    GUI_SETTINGS_FILE_NAME,
    application_base_dir,
    application_data_dir,
    gui_log_path,
    gui_settings_path,
    load_gui_settings,
    save_gui_settings,
    shared_gui_transcript_cache_dir,
)
from tools_for_pharma.oligo.transcript_scan.cli import (
    DEFAULT_CLOSEST_MATCHES,
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
    DEFAULT_BATCH_BASES,
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
    CSV_COLUMNS,
    blast_batch_rows,
    blast_raw_rows,
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
    parse_blast_csv,
    query_length_by_blast_id,
    query_target_summary_rows,
    terminal_table,
    transcript_match_rows,
    transcript_match_terminal_table,
    transcript_matches_to_csv,
    transcript_target_rows,
    write_excel_workbook,
    write_result_workbook,
)
from tools_for_pharma.oligo.transcript_scan.remote_blast import (
    BlastBatchResult,
    append_rid_log,
    combine_blast_csv,
    fasta_record,
    multi_fasta,
    normalize_dna,
)
from tools_for_pharma.oligo.transcript_scan.scanner import (
    DEFAULT_MAX_MISMATCHES,
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
    local_transcript_target,
    prepare_pasted_transcript_sequence,
    read_transcript_input,
    read_transcript_source,
    retrieve_transcript_targets,
    transcript_cache_path,
    transcript_target_from_fasta,
    transcript_target_source,
    validate_single_transcript_record,
)
from tools_for_pharma.oligo.transcript_scan.workflows import (
    SingleSequenceScanConfig,
    run_private_panel_scan,
    run_single_sequence_scan as run_single_sequence_domain_workflow,
)
from tools_for_pharma.sequence.nucleotides import normalize_rna
from tools_for_pharma.shared.excel_utils import list_excel_sheets


DEFAULT_SINGLE_GUI_CLOSEST_MATCHES = 5


def excel_headers(input_file: Path, sheet_name: str | None = None) -> list[str]:
    import pandas as pd

    table = pd.read_excel(input_file, sheet_name=sheet_name or 0, nrows=0)
    return [str(column) for column in table.columns]


def choose_sheet_gui(root, input_file: Path) -> str | None:
    import tkinter as tk
    from tkinter import ttk

    sheets = list_excel_sheets(input_file)
    if len(sheets) <= 1:
        return None

    selected = {"value": sheets[0]}
    window = tk.Toplevel(root)
    window.title("Select transcript scan sheet")
    window.resizable(False, False)
    window.columnconfigure(1, weight=1)

    ttk.Label(window, text="Worksheet").grid(
        row=0, column=0, padx=16, pady=(16, 8), sticky="w"
    )
    sheet_var = tk.StringVar(value=sheets[0])
    sheet_box = ttk.Combobox(
        window,
        textvariable=sheet_var,
        values=sheets,
        state="readonly",
        width=max(30, min(60, max(len(sheet) for sheet in sheets) + 2)),
    )
    sheet_box.grid(row=0, column=1, padx=16, pady=(16, 8), sticky="ew")

    buttons = ttk.Frame(window)
    buttons.grid(row=1, column=0, columnspan=2, padx=16, pady=(8, 16), sticky="e")

    def use_sheet() -> None:
        selected["value"] = sheet_var.get()
        window.destroy()

    def cancel() -> None:
        selected["value"] = None
        window.destroy()

    ttk.Button(buttons, text="Cancel", command=cancel).grid(row=0, column=0, padx=(0, 8))
    ttk.Button(buttons, text="Continue", command=use_sheet).grid(row=0, column=1)
    window.protocol("WM_DELETE_WINDOW", cancel)
    window.bind("<Return>", lambda _event: use_sheet())
    window.bind("<Escape>", lambda _event: cancel())
    window.grab_set()
    sheet_box.focus_set()
    window.wait_window()
    return selected["value"]


def default_header(headers: list[str], candidates: list[str], fallback: str | None = None) -> str:
    by_lower = {header.lower(): header for header in headers}
    for candidate in candidates:
        if candidate.lower() in by_lower:
            return by_lower[candidate.lower()]
    return fallback if fallback is not None else headers[0]


def prompt_and_save_ncbi_email(root, current_email: str = "") -> str | None:
    """Ask for a valid NCBI contact email and save it in portable settings."""
    from tkinter import messagebox, simpledialog

    initial_value = clean_text_for_id(current_email)
    while True:
        entered = simpledialog.askstring(
            "NCBI contact email",
            "Enter your contact email for NCBI transcript requests.\n\n"
            "It is saved locally in TranscriptScanData\\settings.json.",
            initialvalue=initial_value,
            parent=root,
        )
        if entered is None:
            return None
        try:
            email = require_email(entered)
        except ValueError as error:
            messagebox.showerror("Invalid email", str(error), parent=root)
            initial_value = clean_text_for_id(entered)
            continue

        settings = load_gui_settings()
        settings["ncbi_email"] = email
        try:
            save_gui_settings(settings)
        except OSError as error:
            messagebox.showerror(
                "Cannot save settings",
                "The app folder must be writable so the email and transcript "
                f"cache can be saved.\n\n{error}",
                parent=root,
            )
            return None
        return email


def saved_or_prompted_ncbi_email(root) -> str | None:
    """Return the saved email, prompting on first use or invalid settings."""
    saved = clean_text_for_id(load_gui_settings().get("ncbi_email", ""))
    try:
        return require_email(saved)
    except ValueError:
        return prompt_and_save_ncbi_email(root, saved)


def choose_ncbi_gui_mode(root, ncbi_email: str) -> tuple[str | None, str]:
    """Choose the simple single-sequence workflow or the existing Excel workflow."""
    import tkinter as tk
    from tkinter import ttk

    selected = {"mode": None, "email": ncbi_email}
    window = tk.Toplevel(root)
    window.title("Private local transcript scan")
    window.resizable(False, False)

    ttk.Label(
        window,
        text="Choose an input workflow",
        font=("TkDefaultFont", 11, "bold"),
    ).grid(row=0, column=0, columnspan=2, padx=20, pady=(20, 12), sticky="w")
    ttk.Label(
        window,
        text=(
            "Both workflows compare locally. Transcript accessions may be sent "
            "to NCBI, but oligo sequences are not."
        ),
        wraplength=480,
    ).grid(row=1, column=0, columnspan=2, padx=20, pady=(0, 16), sticky="w")

    email_var = tk.StringVar(value=ncbi_email)
    email_frame = ttk.Frame(window)
    email_frame.grid(row=2, column=0, columnspan=2, padx=20, pady=(0, 16), sticky="ew")
    ttk.Label(email_frame, text="NCBI contact email:").grid(row=0, column=0, sticky="w")
    ttk.Label(email_frame, textvariable=email_var).grid(
        row=0, column=1, padx=(6, 12), sticky="w"
    )

    def change_email() -> None:
        changed = prompt_and_save_ncbi_email(window, str(selected["email"]))
        if changed:
            selected["email"] = changed
            email_var.set(changed)

    ttk.Button(email_frame, text="Change", command=change_email).grid(
        row=0, column=2, sticky="e"
    )

    def open_app_data_folder() -> None:
        try:
            data_dir = application_data_dir()
            data_dir.mkdir(parents=True, exist_ok=True)
            os.startfile(str(data_dir))
        except Exception as error:
            from tkinter import messagebox

            messagebox.showerror(
                "Cannot open app data folder",
                str(error),
                parent=window,
            )

    ttk.Button(
        email_frame,
        text="Open app data",
        command=open_app_data_folder,
    ).grid(row=0, column=3, padx=(8, 0), sticky="e")

    def choose(mode: str) -> None:
        selected["mode"] = mode
        window.destroy()

    ttk.Button(
        window,
        text="Single sequence and one transcript",
        command=lambda: choose("single"),
        width=34,
    ).grid(row=3, column=0, padx=(20, 8), pady=(0, 20), sticky="ew")
    ttk.Button(
        window,
        text="Excel sequence table",
        command=lambda: choose("excel"),
        width=26,
    ).grid(row=3, column=1, padx=(8, 20), pady=(0, 20), sticky="ew")

    window.protocol("WM_DELETE_WINDOW", window.destroy)
    window.bind("<Escape>", lambda _event: window.destroy())
    window.grab_set()
    window.wait_window()
    return selected["mode"], str(selected["email"])


def single_sequence_gui_draft(
    previous: dict[str, object] | None = None,
) -> dict[str, object]:
    """Return session-only single-scan values, retaining a previous form draft."""
    draft: dict[str, object] = {
        "sequence_type": "AS",
        "sequence_name": "",
        "sequence": "",
        "target_mode": "accession",
        "target_accession": "",
        "target_name": "",
        "target_sequence": "",
        "target_file": "",
        "scan_regions": ["full"],
        "max_mismatches": DEFAULT_MAX_MISMATCHES,
        "closest": DEFAULT_SINGLE_GUI_CLOSEST_MATCHES,
        "refresh_targets": False,
        "cache_dir": shared_gui_transcript_cache_dir(),
    }
    if previous:
        draft.update(previous)
        draft["scan_regions"] = list(previous.get("scan_regions", ["full"]))
    if draft.get("target_mode") not in {"accession", "paste", "file"}:
        draft["target_mode"] = "accession"
    return draft


def choose_single_sequence_gui_settings(
    root,
    previous: dict[str, object] | None = None,
) -> dict[str, object] | None:
    """Collect one local AS/SS-versus-transcript scan from the user."""
    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk

    selected: dict[str, object] = {}
    draft = single_sequence_gui_draft(previous)
    cache_dir = Path(draft.get("cache_dir") or shared_gui_transcript_cache_dir())
    scan_region_values = {str(value) for value in draft["scan_regions"]}

    window = tk.Toplevel(root)
    window.title("Single sequence transcript scan")
    window.minsize(760, 760)
    window.resizable(True, True)
    window.columnconfigure(1, weight=1)

    sequence_type_var = tk.StringVar(value=str(draft.get("sequence_type") or "AS"))
    sequence_name_var = tk.StringVar(value=str(draft.get("sequence_name") or ""))
    target_mode_var = tk.StringVar(value=str(draft.get("target_mode") or "accession"))
    accession_var = tk.StringVar(value=str(draft.get("target_accession") or ""))
    target_name_var = tk.StringVar(value=str(draft.get("target_name") or ""))
    target_file_var = tk.StringVar(value=str(draft.get("target_file") or ""))
    full_region_var = tk.BooleanVar(value="full" in scan_region_values)
    seed_region_var = tk.BooleanVar(value="seed:2-8" in scan_region_values)
    core_region_var = tk.BooleanVar(value="core:2-18" in scan_region_values)
    max_mismatches_var = tk.StringVar(value=str(draft["max_mismatches"]))
    closest_var = tk.StringVar(value=str(draft["closest"]))
    refresh_var = tk.BooleanVar(value=bool(draft.get("refresh_targets", False)))

    ttk.Label(window, text="Sequence type").grid(
        row=0, column=0, padx=16, pady=(16, 8), sticky="w"
    )
    ttk.Combobox(
        window,
        textvariable=sequence_type_var,
        values=["AS", "SS"],
        state="readonly",
        width=10,
    ).grid(row=0, column=1, padx=16, pady=(16, 8), sticky="w")

    ttk.Label(window, text="Sequence name (optional)").grid(
        row=1, column=0, padx=16, pady=8, sticky="w"
    )
    ttk.Entry(window, textvariable=sequence_name_var, width=52).grid(
        row=1, column=1, padx=16, pady=8, sticky="ew"
    )

    ttk.Label(window, text="Sequence, 5' to 3'").grid(
        row=2, column=0, padx=16, pady=8, sticky="nw"
    )
    sequence_text = tk.Text(window, width=56, height=4, wrap="word")
    sequence_text.grid(row=2, column=1, padx=16, pady=8, sticky="ew")
    sequence_text.insert("1.0", str(draft.get("sequence") or ""))

    target_frame = ttk.LabelFrame(window, text="Transcript target source")
    target_frame.grid(row=3, column=0, columnspan=2, padx=16, pady=8, sticky="nsew")
    target_frame.columnconfigure(1, weight=1)

    ttk.Radiobutton(
        target_frame,
        text="Download/reuse a RefSeq transcript accession",
        variable=target_mode_var,
        value="accession",
    ).grid(row=0, column=0, columnspan=3, padx=12, pady=(10, 4), sticky="w")
    ttk.Label(target_frame, text="Exact NM/XM/NR/XR version").grid(
        row=1, column=0, padx=(28, 8), pady=4, sticky="w"
    )
    accession_entry = ttk.Entry(target_frame, textvariable=accession_var, width=32)
    accession_entry.grid(row=1, column=1, padx=8, pady=4, sticky="ew")
    ttk.Label(target_frame, text="Example: NM_002439.5").grid(
        row=2, column=1, padx=8, pady=(0, 8), sticky="w"
    )

    ttk.Radiobutton(
        target_frame,
        text="Paste one transcript sequence (plain sequence or one FASTA record)",
        variable=target_mode_var,
        value="paste",
    ).grid(row=3, column=0, columnspan=3, padx=12, pady=(8, 4), sticky="w")
    ttk.Label(target_frame, text="Target name (optional)").grid(
        row=4, column=0, padx=(28, 8), pady=4, sticky="w"
    )
    target_name_entry = ttk.Entry(target_frame, textvariable=target_name_var, width=32)
    target_name_entry.grid(row=4, column=1, padx=8, pady=4, sticky="ew")
    target_sequence_text = tk.Text(target_frame, width=56, height=6, wrap="word")
    target_sequence_text.grid(
        row=5,
        column=0,
        columnspan=3,
        padx=(28, 12),
        pady=(4, 8),
        sticky="ew",
    )
    target_sequence_text.insert("1.0", str(draft.get("target_sequence") or ""))

    ttk.Radiobutton(
        target_frame,
        text="Load one transcript FASTA/text file",
        variable=target_mode_var,
        value="file",
    ).grid(row=6, column=0, columnspan=3, padx=12, pady=(8, 4), sticky="w")
    target_file_entry = ttk.Entry(target_frame, textvariable=target_file_var, width=42)
    target_file_entry.grid(row=7, column=0, columnspan=2, padx=(28, 8), pady=(4, 10), sticky="ew")

    def choose_target_file() -> None:
        path = filedialog.askopenfilename(
            parent=window,
            title="Select one transcript FASTA or text file",
            filetypes=[
                ("Sequence text files", "*.txt *.fa *.fasta *.fna *.ffn"),
                ("All files", "*.*"),
            ],
        )
        if path:
            target_file_var.set(path)
            target_mode_var.set("file")
            update_target_controls()

    browse_button = ttk.Button(target_frame, text="Browse", command=choose_target_file)
    browse_button.grid(row=7, column=2, padx=(0, 12), pady=(4, 10), sticky="e")

    ttk.Label(window, text="Scan regions").grid(
        row=4, column=0, padx=16, pady=8, sticky="nw"
    )
    region_frame = ttk.Frame(window)
    region_frame.grid(row=4, column=1, padx=16, pady=8, sticky="w")
    ttk.Checkbutton(region_frame, text="Full sequence", variable=full_region_var).grid(
        row=0, column=0, padx=(0, 16), sticky="w"
    )
    ttk.Checkbutton(
        region_frame,
        text="Seed (positions 2-8)",
        variable=seed_region_var,
    ).grid(row=0, column=1, padx=(0, 16), sticky="w")
    ttk.Checkbutton(
        region_frame,
        text="Core (positions 2-18)",
        variable=core_region_var,
    ).grid(row=0, column=2, sticky="w")

    ttk.Label(window, text="Maximum mismatches").grid(
        row=5, column=0, padx=16, pady=8, sticky="w"
    )
    ttk.Entry(window, textvariable=max_mismatches_var, width=8).grid(
        row=5, column=1, padx=16, pady=8, sticky="w"
    )

    ttk.Label(window, text="Closest windows per region").grid(
        row=6, column=0, padx=16, pady=8, sticky="w"
    )
    ttk.Entry(window, textvariable=closest_var, width=8).grid(
        row=6, column=1, padx=16, pady=8, sticky="w"
    )

    option_frame = ttk.Frame(window)
    option_frame.grid(row=7, column=0, columnspan=2, padx=16, pady=8, sticky="w")
    ttk.Label(
        option_frame,
        text=(
            "Accession mode uses the saved transcript when available and downloads "
            "it automatically when missing. Pasted and file targets stay local."
        ),
        wraplength=680,
    ).grid(row=0, column=0, sticky="w")
    refresh_check = ttk.Checkbutton(
        option_frame,
        text="Refresh this transcript from NCBI once",
        variable=refresh_var,
    )
    refresh_check.grid(row=1, column=0, pady=(4, 0), sticky="w")

    cache_frame = ttk.Frame(window)
    cache_frame.grid(row=8, column=0, columnspan=2, padx=16, pady=8, sticky="ew")
    ttk.Label(
        cache_frame,
        text=f"Shared transcript cache: {cache_dir}",
        wraplength=560,
    ).grid(row=0, column=0, padx=(0, 8), sticky="w")

    def open_cache_folder() -> None:
        try:
            cache_dir.mkdir(parents=True, exist_ok=True)
            os.startfile(str(cache_dir))
        except Exception as error:
            messagebox.showerror("Cannot open cache folder", str(error), parent=window)

    ttk.Button(cache_frame, text="Open cache folder", command=open_cache_folder).grid(
        row=0, column=1, sticky="e"
    )
    ttk.Label(
        window,
        text=(
            "Privacy: the oligo sequence and any pasted/local transcript stay on "
            "this computer. In accession mode, only accession/contact metadata are "
            "sent to NCBI. Form values are discarded when the app closes."
        ),
        wraplength=680,
    ).grid(row=9, column=0, columnspan=2, padx=16, pady=(4, 8), sticky="w")

    buttons = ttk.Frame(window)
    buttons.grid(row=10, column=0, columnspan=2, padx=16, pady=(8, 16), sticky="e")

    def update_target_controls() -> None:
        mode = target_mode_var.get()
        accession_entry.configure(state="normal" if mode == "accession" else "disabled")
        target_name_entry.configure(state="normal" if mode == "paste" else "disabled")
        target_sequence_text.configure(state="normal" if mode == "paste" else "disabled")
        target_file_entry.configure(state="normal" if mode == "file" else "disabled")
        browse_button.configure(state="normal" if mode == "file" else "disabled")
        refresh_check.configure(state="normal" if mode == "accession" else "disabled")

    target_mode_var.trace_add("write", lambda *_args: update_target_controls())
    update_target_controls()

    def use_settings() -> None:
        try:
            sequence = normalize_rna(sequence_text.get("1.0", "end").strip())
            mode = target_mode_var.get()
            if mode == "accession":
                accession = normalize_versioned_refseq_accession(accession_var.get())
            elif mode == "paste":
                prepare_pasted_transcript_sequence(
                    target_sequence_text.get("1.0", "end").strip(),
                    target_name_var.get(),
                )
                accession = accession_var.get().strip()
            elif mode == "file":
                target_file = Path(target_file_var.get().strip())
                if not target_file.exists() or not target_file.is_file():
                    raise ValueError(f"Transcript file does not exist: {target_file}")
                file_text = target_file.read_text(encoding="utf-8-sig")
                validate_single_transcript_record(file_text, str(target_file))
                fasta_or_plain_text_to_sequence(file_text)
                accession = accession_var.get().strip()
            else:
                raise ValueError("Choose a transcript target source.")

            scan_regions = []
            if full_region_var.get():
                scan_regions.append("full")
            if seed_region_var.get():
                scan_regions.append("seed:2-8")
            if core_region_var.get():
                scan_regions.append("core:2-18")
            if not scan_regions:
                raise ValueError("Select at least one scan region.")
            for region in parse_scan_regions(scan_regions):
                antisense_region_sequence(sequence, region)
            max_mismatches = int(max_mismatches_var.get())
            if max_mismatches < 0:
                raise ValueError("Maximum mismatches must be 0 or greater.")
            closest = int(closest_var.get())
            if closest < 1:
                raise ValueError("Closest windows must be 1 or greater.")
        except (OSError, UnicodeError, ValueError) as error:
            messagebox.showerror("Invalid settings", str(error), parent=window)
            return

        selected.update(
            {
                "sequence_type": sequence_type_var.get(),
                "sequence_name": sequence_name_var.get().strip(),
                "sequence": sequence,
                "target_mode": mode,
                "target_accession": accession,
                "target_name": target_name_var.get().strip(),
                "target_sequence": target_sequence_text.get("1.0", "end").strip(),
                "target_file": target_file_var.get().strip(),
                "scan_regions": scan_regions,
                "max_mismatches": max_mismatches,
                "closest": closest,
                "refresh_targets": refresh_var.get() if mode == "accession" else False,
                "cache_dir": cache_dir,
            }
        )
        window.destroy()

    ttk.Button(buttons, text="Cancel", command=window.destroy).grid(
        row=0, column=0, padx=(0, 8)
    )
    ttk.Button(buttons, text="Run local scan", command=use_settings).grid(
        row=0, column=1
    )
    window.protocol("WM_DELETE_WINDOW", window.destroy)
    window.bind("<Escape>", lambda _event: window.destroy())
    window.grab_set()
    window.wait_window()
    return selected or None


def single_sequence_gui_args(settings: dict[str, object]) -> argparse.Namespace:
    """Translate validated single-sequence GUI settings into normal CLI arguments."""
    sequence_type = normalize_sequence_type(str(settings["sequence_type"]))
    sequence_flag = "--as-sequence" if sequence_type == "AS" else "--ss-sequence"
    name_flag = "--as-name" if sequence_type == "AS" else "--ss-name"
    target_mode = str(settings.get("target_mode") or "accession")
    argv = [
        sequence_flag,
        str(settings["sequence"]),
        "--cache-dir",
        str(settings.get("cache_dir") or shared_gui_transcript_cache_dir()),
        "--email",
        str(settings.get("email", "")),
        "--max-mismatches",
        str(settings["max_mismatches"]),
        "--closest",
        str(settings.get("closest", DEFAULT_SINGLE_GUI_CLOSEST_MATCHES)),
    ]
    if target_mode == "accession":
        argv.extend(
            [
                "--private-panel",
                "--target-accession",
                normalize_versioned_refseq_accession(settings.get("target_accession", "")),
            ]
        )
    elif target_mode == "paste":
        argv.extend(
            [
                "--target-sequence",
                prepare_pasted_transcript_sequence(
                    str(settings.get("target_sequence") or ""),
                    str(settings.get("target_name") or ""),
                ),
            ]
        )
    elif target_mode == "file":
        target_file = Path(str(settings.get("target_file") or ""))
        if not str(settings.get("target_file") or "").strip():
            raise ValueError("Choose one transcript FASTA/text file.")
        argv.extend(["--target-file", str(target_file)])
    else:
        raise ValueError(f"Unsupported single-sequence target mode: {target_mode}")
    if settings.get("sequence_name"):
        argv.extend([name_flag, str(settings["sequence_name"])])
    for region in settings.get("scan_regions", ["full"]):
        argv.extend(["--scan-region", str(region)])
    if target_mode == "accession" and settings.get("refresh_targets"):
        argv.append("--refresh-targets")
    return build_parser().parse_args(argv)


def run_single_sequence_scan(
    args: argparse.Namespace,
    *,
    progress_callback: Callable[[int, int, str, str], None] | None = None,
    client: NcbiHttpClient | None = None,
) -> tuple[list[AntisenseQuery], list[AntisenseRegion], PrivatePanelScanResult]:
    """Run one private local query against one accession, pasted, or file target."""
    accessions = target_accession_values(args.target_accession)
    if accessions:
        if len(accessions) != 1:
            raise ValueError("Single-sequence mode requires exactly one transcript accession.")
        source = AccessionTargetSource(
            accession=normalize_versioned_refseq_accession(accessions[0]),
            cache_dir=private_panel_cache_dir(args),
            offline=False,
            refresh=bool(getattr(args, "refresh_targets", False)),
        )
    else:
        source = transcript_target_source(
            transcript_sequence=args.target_sequence,
            transcript_file=args.target_file,
        )

    queries = args_antisense_queries(args)
    if len(queries) != 1:
        raise ValueError("Single-sequence mode requires exactly one AS or SS sequence.")
    scan_regions = parse_scan_regions(args.scan_region)
    result = run_single_sequence_domain_workflow(
        SingleSequenceScanConfig(
            target_source=source,
            email=args.email,
            tool=args.tool,
            request_seconds=args.request_seconds,
            max_mismatches=args.max_mismatches,
            closest=args.closest,
        ),
        queries[0],
        scan_regions,
        progress_callback=progress_callback,
        client=client,
    )
    return queries, scan_regions, result


def show_single_sequence_result_gui(root, result_text: str) -> bool:
    """Show a scrollable result and return True when the user requests another scan."""
    import tkinter as tk
    from tkinter import ttk

    new_scan = {"value": False}
    window = tk.Toplevel(root)
    window.title("Single sequence transcript scan result")
    window.geometry("1100x700")
    window.minsize(760, 480)
    window.rowconfigure(0, weight=1)
    window.columnconfigure(0, weight=1)

    result_frame = ttk.Frame(window)
    result_frame.grid(row=0, column=0, padx=12, pady=12, sticky="nsew")
    result_frame.rowconfigure(0, weight=1)
    result_frame.columnconfigure(0, weight=1)
    output_text = tk.Text(
        result_frame,
        wrap="none",
        font="TkFixedFont",
        padx=8,
        pady=8,
    )
    vertical = ttk.Scrollbar(result_frame, orient="vertical", command=output_text.yview)
    horizontal = ttk.Scrollbar(result_frame, orient="horizontal", command=output_text.xview)
    output_text.configure(yscrollcommand=vertical.set, xscrollcommand=horizontal.set)
    output_text.grid(row=0, column=0, sticky="nsew")
    vertical.grid(row=0, column=1, sticky="ns")
    horizontal.grid(row=1, column=0, sticky="ew")
    output_text.insert("1.0", result_text)
    output_text.configure(state="disabled")

    buttons = ttk.Frame(window)
    buttons.grid(row=1, column=0, padx=12, pady=(0, 12), sticky="e")

    def copy_all() -> None:
        root.clipboard_clear()
        root.clipboard_append(result_text)
        root.update()

    def request_new_scan() -> None:
        new_scan["value"] = True
        window.destroy()

    ttk.Button(buttons, text="Copy all", command=copy_all).grid(
        row=0, column=0, padx=(0, 8)
    )
    ttk.Button(buttons, text="Edit and run again", command=request_new_scan).grid(
        row=0, column=1, padx=(0, 8)
    )
    ttk.Button(buttons, text="Close", command=window.destroy).grid(row=0, column=2)

    window.protocol("WM_DELETE_WINDOW", window.destroy)
    window.bind("<Escape>", lambda _event: window.destroy())
    window.grab_set()
    window.wait_window()
    return new_scan["value"]


def run_single_sequence_gui(root, ncbi_email: str) -> int:
    """Run the one-sequence/one-transcript private local GUI workflow."""
    import tkinter as tk
    from tkinter import messagebox, ttk

    draft = single_sequence_gui_draft()
    while True:
        settings = choose_single_sequence_gui_settings(root, draft)
        if not settings:
            return 0
        settings["email"] = ncbi_email
        draft = {key: value for key, value in settings.items() if key != "email"}
        try:
            args = single_sequence_gui_args(settings)
            validate_runtime_args(args)
        except Exception as error:
            logging.exception("Invalid single-sequence scan settings")
            messagebox.showerror(
                "Single sequence transcript scan failed",
                str(error),
                parent=root,
            )
            continue

        progress_window = tk.Toplevel(root)
        progress_window.title("Single sequence transcript scan")
        progress_window.resizable(False, False)
        status_var = tk.StringVar(value="Preparing transcript reference...")
        ttk.Label(progress_window, textvariable=status_var, width=64).grid(
            row=0, column=0, padx=16, pady=(16, 8), sticky="w"
        )
        progress_bar = ttk.Progressbar(
            progress_window,
            mode="determinate",
            length=440,
            maximum=1,
        )
        progress_bar.grid(row=1, column=0, padx=16, pady=(8, 16), sticky="ew")

        def update_progress(
            completed: int,
            total: int,
            accession: str,
            status: str,
        ) -> None:
            progress_bar.configure(maximum=max(total, 1), value=completed)
            status_var.set(f"{accession}: {status}")
            progress_window.update()

        try:
            try:
                queries, scan_regions, result = run_single_sequence_scan(
                    args,
                    progress_callback=update_progress,
                )
            finally:
                progress_window.destroy()
        except Exception as error:
            logging.exception("Single-sequence transcript scan failed")
            messagebox.showerror(
                "Single sequence transcript scan failed",
                f"{error}\n\nYour form entries have been retained.",
                parent=root,
            )
            continue

        result_text = format_single_sequence_scan_result(
            args,
            queries,
            scan_regions,
            result,
        )
        draft["refresh_targets"] = False
        if not show_single_sequence_result_gui(root, result_text):
            return 0


def choose_ncbi_gui_settings(root, headers: list[str]) -> dict[str, object] | None:
    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk

    selected: dict[str, object] = {}
    target_accession_default = default_header(
        headers,
        ["target_accession", "target accession", "refseq", "refseq accession"],
        "",
    )

    window = tk.Toplevel(root)
    window.title("Private local transcript scan settings")
    window.resizable(False, False)
    window.columnconfigure(1, weight=1)

    sequence_type_var = tk.StringVar(value="AS")
    as_column_var = tk.StringVar(
        value=default_header(headers, ["AS_5to3", "antisense", "as", "sequence"])
    )
    name_column_var = tk.StringVar(
        value=default_header(headers, ["oligo_id", "oligo id", "id", "name"], "")
    )
    target_mode_var = tk.StringVar(
        value="column" if target_accession_default else "refseq"
    )
    target_column_var = tk.StringVar(value=target_accession_default or headers[0])
    refseq_var = tk.StringVar()
    transcript_file_var = tk.StringVar()
    panel_accessions_var = tk.StringVar()
    target_table_var = tk.StringVar()
    scan_regions_var = tk.StringVar(value="full")
    max_mismatches_var = tk.StringVar(value=str(DEFAULT_MAX_MISMATCHES))
    offline_var = tk.BooleanVar(value=False)

    ttk.Label(window, text="Sequence type").grid(
        row=0, column=0, padx=16, pady=(16, 8), sticky="w"
    )
    ttk.Combobox(
        window,
        textvariable=sequence_type_var,
        values=["AS", "SS"],
        state="readonly",
        width=12,
    ).grid(row=0, column=1, padx=16, pady=(16, 8), sticky="w")

    ttk.Label(window, text="Sequence column").grid(
        row=1, column=0, padx=16, pady=8, sticky="w"
    )
    ttk.Combobox(
        window,
        textvariable=as_column_var,
        values=headers,
        state="readonly",
        width=max(30, min(60, max(len(header) for header in headers) + 2)),
    ).grid(row=1, column=1, padx=16, pady=8, sticky="ew")

    ttk.Label(window, text="Sequence name column").grid(
        row=2, column=0, padx=16, pady=8, sticky="w"
    )
    ttk.Combobox(
        window,
        textvariable=name_column_var,
        values=["", *headers],
        state="readonly",
        width=36,
    ).grid(row=2, column=1, padx=16, pady=8, sticky="ew")

    ttk.Label(window, text="Transcript source").grid(
        row=3, column=0, padx=16, pady=8, sticky="nw"
    )
    source_frame = ttk.Frame(window)
    source_frame.grid(row=3, column=1, padx=16, pady=8, sticky="ew")
    source_frame.columnconfigure(1, weight=1)

    ttk.Radiobutton(
        source_frame,
        text="Use target_accession column",
        variable=target_mode_var,
        value="column",
    ).grid(row=0, column=0, columnspan=2, sticky="w")
    ttk.Combobox(
        source_frame,
        textvariable=target_column_var,
        values=headers,
        state="readonly",
        width=36,
    ).grid(row=1, column=0, columnspan=2, pady=(4, 8), sticky="ew")

    ttk.Radiobutton(
        source_frame,
        text="Use this RefSeq accession for all rows",
        variable=target_mode_var,
        value="refseq",
    ).grid(row=2, column=0, columnspan=2, sticky="w")
    ttk.Entry(source_frame, textvariable=refseq_var, width=38).grid(
        row=3, column=0, columnspan=2, pady=(4, 8), sticky="ew"
    )

    ttk.Radiobutton(
        source_frame,
        text="Use transcript FASTA/text file for all rows",
        variable=target_mode_var,
        value="file",
    ).grid(row=4, column=0, columnspan=2, sticky="w")

    def choose_transcript_file() -> None:
        path = filedialog.askopenfilename(
            parent=window,
            title="Select transcript FASTA or text file",
            filetypes=[
                ("Sequence text files", "*.txt *.fa *.fasta *.fna *.ffn"),
                ("FASTA files", "*.fa *.fasta *.fna *.ffn"),
                ("Text files", "*.txt"),
                ("All files", "*.*"),
            ],
        )
        if path:
            transcript_file_var.set(path)
            target_mode_var.set("file")

    ttk.Entry(source_frame, textvariable=transcript_file_var, width=38).grid(
        row=5, column=0, pady=(4, 0), sticky="ew"
    )
    ttk.Button(source_frame, text="Browse", command=choose_transcript_file).grid(
        row=5, column=1, padx=(8, 0), pady=(4, 0), sticky="e"
    )

    ttk.Radiobutton(
        source_frame,
        text="Private panel: versioned RefSeq accessions",
        variable=target_mode_var,
        value="panel",
    ).grid(row=6, column=0, columnspan=2, pady=(8, 0), sticky="w")
    ttk.Entry(source_frame, textvariable=panel_accessions_var, width=38).grid(
        row=7, column=0, columnspan=2, pady=(4, 8), sticky="ew"
    )

    ttk.Radiobutton(
        source_frame,
        text="Private panel: accession table",
        variable=target_mode_var,
        value="panel_table",
    ).grid(row=8, column=0, columnspan=2, sticky="w")

    def choose_target_table() -> None:
        path = filedialog.askopenfilename(
            parent=window,
            title="Select target accession table",
            filetypes=[
                ("Target lists", "*.txt *.list *.csv *.tsv *.xlsx *.xlsm *.xls"),
                ("All files", "*.*"),
            ],
        )
        if path:
            target_table_var.set(path)
            target_mode_var.set("panel_table")

    ttk.Entry(source_frame, textvariable=target_table_var, width=38).grid(
        row=9, column=0, pady=(4, 0), sticky="ew"
    )
    ttk.Button(source_frame, text="Browse", command=choose_target_table).grid(
        row=9, column=1, padx=(8, 0), pady=(4, 0), sticky="e"
    )

    ttk.Label(window, text="Scan regions").grid(
        row=4, column=0, padx=16, pady=8, sticky="w"
    )
    ttk.Entry(window, textvariable=scan_regions_var, width=38).grid(
        row=4, column=1, padx=16, pady=8, sticky="ew"
    )

    ttk.Label(window, text="Max mismatches").grid(
        row=5, column=0, padx=16, pady=8, sticky="w"
    )
    ttk.Entry(window, textvariable=max_mismatches_var, width=12).grid(
        row=5, column=1, padx=16, pady=8, sticky="w"
    )

    ttk.Checkbutton(
        window,
        text="Offline: require every panel transcript in the local cache",
        variable=offline_var,
    ).grid(row=6, column=0, columnspan=2, padx=16, pady=8, sticky="w")
    ttk.Label(
        window,
        text="Private panel mode sends transcript accessions to NCBI; guide sequences remain local.",
    ).grid(row=7, column=0, columnspan=2, padx=16, pady=(0, 8), sticky="w")

    buttons = ttk.Frame(window)
    buttons.grid(row=8, column=0, columnspan=2, padx=16, pady=(8, 16), sticky="e")

    def use_settings() -> None:
        try:
            max_mismatches = int(max_mismatches_var.get())
            if max_mismatches < 0:
                raise ValueError
            scan_regions = [
                value.strip()
                for value in re.split(r"[;,]", scan_regions_var.get())
                if value.strip()
            ]
            parse_scan_regions(scan_regions)
        except ValueError as error:
            messagebox.showerror("Invalid settings", str(error), parent=window)
            return

        mode = target_mode_var.get()
        if mode == "column" and not target_column_var.get():
            messagebox.showerror("Missing target accession column", "Choose a target accession column.", parent=window)
            return
        if mode == "refseq" and not refseq_var.get().strip():
            messagebox.showerror("Missing RefSeq", "Enter a RefSeq accession.", parent=window)
            return
        if mode == "file" and not transcript_file_var.get().strip():
            messagebox.showerror("Missing transcript file", "Choose a transcript FASTA/text file.", parent=window)
            return
        panel_accessions = [
            value
            for value in re.split(r"[,;\s]+", panel_accessions_var.get().strip())
            if value
        ]
        if mode == "panel":
            if not panel_accessions:
                messagebox.showerror(
                    "Missing panel accessions",
                    "Enter one or more versioned RefSeq accessions.",
                    parent=window,
                )
                return
            try:
                panel_accessions = [
                    normalize_versioned_refseq_accession(value)
                    for value in panel_accessions
                ]
            except ValueError as error:
                messagebox.showerror("Invalid panel accession", str(error), parent=window)
                return
        if mode == "panel_table" and not target_table_var.get().strip():
            messagebox.showerror(
                "Missing target table",
                "Choose a transcript accession table.",
                parent=window,
            )
            return
        if offline_var.get() and mode not in {"panel", "panel_table"}:
            messagebox.showerror(
                "Offline panel mode",
                "Offline mode is available only for private transcript panels.",
                parent=window,
            )
            return

        selected.update(
            {
                "sequence_type": sequence_type_var.get(),
                "as_column": as_column_var.get(),
                "as_name_column": name_column_var.get() or None,
                "target_mode": mode,
                "target_accession_column": target_column_var.get() if mode == "column" else None,
                "target_accession": (
                    refseq_var.get().strip()
                    if mode == "refseq"
                    else panel_accessions if mode == "panel" else None
                ),
                "target_file": Path(transcript_file_var.get()) if mode == "file" else None,
                "target_table": (
                    Path(target_table_var.get()) if mode == "panel_table" else None
                ),
                "private_panel": mode in {"panel", "panel_table"},
                "offline": offline_var.get(),
                "scan_regions": scan_regions,
                "max_mismatches": max_mismatches,
            }
        )
        window.destroy()

    def cancel() -> None:
        selected.clear()
        window.destroy()

    ttk.Button(buttons, text="Cancel", command=cancel).grid(row=0, column=0, padx=(0, 8))
    ttk.Button(buttons, text="Run", command=use_settings).grid(row=0, column=1)

    window.protocol("WM_DELETE_WINDOW", cancel)
    window.bind("<Return>", lambda _event: use_settings())
    window.bind("<Escape>", lambda _event: cancel())
    window.grab_set()
    window.wait_window()
    return selected or None


def gui_args(input_file: Path, sheet_name: str | None, settings: dict[str, object]) -> argparse.Namespace:
    sequence_type = normalize_sequence_type(str(settings.get("sequence_type", "AS")))
    private_panel = bool(settings.get("private_panel", False))
    output_path = (
        input_file.with_name(f"{input_file.stem}_private_transcript_panel_results.xlsx")
        if private_panel
        else default_gui_result_workbook(input_file)
    )
    return argparse.Namespace(
        as_sequence=None,
        as_name=None,
        as_file=None,
        as_table=input_file if sequence_type == "AS" else None,
        as_column=settings["as_column"] if sequence_type == "AS" else None,
        as_name_column=settings["as_name_column"] if sequence_type == "AS" else None,
        as_sheet=sheet_name if sequence_type == "AS" else None,
        ss_sequence=None,
        ss_name=None,
        ss_file=None,
        ss_table=input_file if sequence_type == "SS" else None,
        ss_column=settings["as_column"] if sequence_type == "SS" else None,
        ss_name_column=settings["as_name_column"] if sequence_type == "SS" else None,
        ss_sheet=sheet_name if sequence_type == "SS" else None,
        target_accession=settings["target_accession"],
        target_accession_column=settings["target_accession_column"],
        target_file=settings["target_file"],
        target_sequence=None,
        target_table=settings.get("target_table"),
        target_column=None,
        target_sheet=None,
        private_panel=private_panel,
        offline=bool(settings.get("offline", False)),
        refresh_targets=False,
        download_targets_only=False,
        scan_region=settings["scan_regions"],
        max_mismatches=settings["max_mismatches"],
        email=str(settings.get("email", "")),
        tool=DEFAULT_TOOL,
        blast=False,
        blast_only=False,
        database=DEFAULT_DATABASE,
        expect=DEFAULT_EXPECT,
        word_size=DEFAULT_WORD_SIZE,
        hitlist_size=DEFAULT_HITLIST_SIZE,
        megablast=False,
        timeout_seconds=1800,
        max_batch_bases=DEFAULT_BATCH_BASES,
        request_seconds=DEFAULT_REQUEST_SECONDS,
        poll_seconds=DEFAULT_POLL_SECONDS,
        filter_max_mismatches=DEFAULT_MAX_MISMATCHES,
        filter_max_gap_opens=0,
        filter_min_alignment_fraction=0.8,
        cache_dir=shared_gui_transcript_cache_dir(),
        rid_log=None,
        output=None,
        blast_output=None,
        terminal=False,
        stdout_csv=False,
        closest=None,
        result_workbook=output_path,
        gui=False,
    )


def run_gui() -> int:
    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk

    root = tk.Tk()
    root.withdraw()

    try:
        ncbi_email = saved_or_prompted_ncbi_email(root)
        if ncbi_email is None:
            return 0
        gui_mode, ncbi_email = choose_ncbi_gui_mode(root, ncbi_email)
        if gui_mode is None:
            return 0
        if gui_mode == "single":
            return run_single_sequence_gui(root, ncbi_email)

        input_path = filedialog.askopenfilename(
            title="Select local transcript scan Excel input file",
            filetypes=[
                ("Excel files", "*.xlsx *.xlsm *.xls"),
                ("All files", "*.*"),
            ],
        )
        if not input_path:
            return 0

        input_file = Path(input_path)
        sheet_name = choose_sheet_gui(root, input_file)
        if sheet_name is None and len(list_excel_sheets(input_file)) > 1:
            return 0

        headers = excel_headers(input_file, sheet_name)
        if not headers:
            messagebox.showerror("No columns", "The selected sheet has no headers.")
            return 1

        settings = choose_ncbi_gui_settings(root, headers)
        if not settings:
            return 0
        settings["email"] = ncbi_email

        args = gui_args(input_file, sheet_name, settings)
        validate_runtime_args(args)
        started_at = datetime.now(timezone.utc).isoformat()
        if private_panel_requested(args):
            progress_window = tk.Toplevel(root)
            progress_window.title("Private transcript panel")
            progress_window.resizable(False, False)
            status_var = tk.StringVar(value="Preparing transcript references...")
            ttk.Label(progress_window, textvariable=status_var, width=64).grid(
                row=0,
                column=0,
                padx=16,
                pady=(16, 8),
                sticky="w",
            )
            progress_bar = ttk.Progressbar(
                progress_window,
                mode="determinate",
                length=440,
            )
            progress_bar.grid(row=1, column=0, padx=16, pady=8, sticky="ew")
            cancelled = {"value": False}

            def request_cancel() -> None:
                cancelled["value"] = True
                status_var.set("Cancelling after the current retrieval request...")

            ttk.Button(
                progress_window,
                text="Cancel",
                command=request_cancel,
            ).grid(row=2, column=0, padx=16, pady=(8, 16), sticky="e")
            progress_window.protocol("WM_DELETE_WINDOW", request_cancel)

            def update_progress(
                completed: int,
                total: int,
                accession: str,
                status: str,
            ) -> None:
                progress_bar.configure(maximum=max(total, 1), value=completed)
                status_var.set(
                    f"{completed}/{total} targets - {accession}: {status}"
                )
                progress_window.update()

            try:
                exit_code = run_private_panel_workflow(
                    args,
                    started_at,
                    progress_callback=update_progress,
                    cancel_check=lambda: cancelled["value"],
                    include_comparison_results=True,
                )
            finally:
                progress_window.destroy()
            if cancelled["value"]:
                messagebox.showwarning(
                    "Private transcript panel cancelled",
                    "The partial status workbook was saved. Guide sequences "
                    f"remained local.\n\n{args.result_workbook}",
                )
            elif exit_code == 0:
                messagebox.showinfo(
                    "Private transcript panel complete",
                    "Guide sequences remained local.\n\n"
                    f"Wrote panel workbook to:\n{args.result_workbook}",
                )
            elif exit_code == 2:
                messagebox.showwarning(
                    "Private transcript panel completed with target errors",
                    "Some targets could not be scanned. Review comparison_results "
                    f"and transcript_targets in:\n{args.result_workbook}",
                )
            else:
                messagebox.showwarning(
                    "Private transcript panel incomplete",
                    "No target transcript was ready. Review the transcript_targets "
                    f"sheet in:\n{args.result_workbook}",
                )
            return exit_code
        queries = args_antisense_queries(args)
        scan_regions = parse_scan_regions(args.scan_region)
        local_matches, comparison_results = run_local_scan_with_comparison(
            args,
            queries,
            scan_regions,
        )
        completed_at = datetime.now(timezone.utc).isoformat()
        write_result_workbook(
            args.result_workbook,
            args,
            queries,
            scan_regions,
            local_matches,
            [],
            started_at,
            completed_at,
            include_blast_sheets=False,
            comparison_results=comparison_results,
        )
        messagebox.showinfo(
            "Local transcript scan complete",
            f"Wrote local transcript scan workbook to:\n{args.result_workbook}",
        )
        return 0
    except Exception as error:
        logging.exception("Local transcript scan failed")
        messagebox.showerror(
            "Local transcript scan failed",
            f"{error}\n\nDiagnostic log:\n{gui_log_path()}",
        )
        return 1
    finally:
        root.destroy()
