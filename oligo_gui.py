"""Small tkinter GUI for oligo processing."""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Tuple

from oligo_core import DEFAULT_END, DEFAULT_START
from oligo_table import process_table, read_table
from oligo_transcript import format_transcript_oligo_result, transcript_file_to_oligo
from oligo_transcript_table import (
    default_transcript_range_output_path,
    process_transcript_range_table,
)


def _choose_settings_gui(
    root,
    headers: List[str],
) -> Optional[Tuple[str, int, int]]:
    """Ask the user for column, start, and end positions."""
    import tkinter as tk
    from tkinter import messagebox, ttk

    selected = {"value": None}
    default_column = next(
        (header for header in headers if header.strip().lower() == "antisense"),
        next((header for header in headers if "antisense" in header.lower()), headers[0]),
    )

    window = tk.Toplevel(root)
    window.title("Oligo settings")
    window.resizable(False, False)
    window.columnconfigure(1, weight=1)

    ttk.Label(window, text="Antisense column").grid(
        row=0, column=0, padx=16, pady=(16, 8), sticky="w"
    )
    column_var = tk.StringVar(value=default_column)
    column_box = ttk.Combobox(
        window,
        textvariable=column_var,
        values=headers,
        state="readonly",
        width=max(30, min(60, max(len(header) for header in headers) + 2)),
    )
    column_box.grid(row=0, column=1, padx=16, pady=(16, 8), sticky="ew")

    ttk.Label(window, text="Start position").grid(
        row=1, column=0, padx=16, pady=8, sticky="w"
    )
    start_var = tk.StringVar(value=str(DEFAULT_START))
    ttk.Entry(window, textvariable=start_var, width=10).grid(
        row=1, column=1, padx=16, pady=8, sticky="w"
    )

    ttk.Label(window, text="End position").grid(
        row=2, column=0, padx=16, pady=8, sticky="w"
    )
    end_var = tk.StringVar(value=str(DEFAULT_END))
    ttk.Entry(window, textvariable=end_var, width=10).grid(
        row=2, column=1, padx=16, pady=8, sticky="w"
    )

    buttons = ttk.Frame(window)
    buttons.grid(row=3, column=0, columnspan=2, padx=16, pady=(8, 16), sticky="e")

    def use_settings() -> None:
        try:
            start = int(start_var.get())
            end = int(end_var.get())
            if start < 1 or end < start:
                raise ValueError
        except ValueError:
            messagebox.showerror(
                "Invalid positions",
                "Start and end must be whole numbers, with end >= start.",
                parent=window,
            )
            return

        selected["value"] = (column_var.get(), start, end)
        window.destroy()

    def cancel() -> None:
        window.destroy()

    ttk.Button(buttons, text="Cancel", command=cancel).grid(row=0, column=0, padx=(0, 8))
    ttk.Button(buttons, text="Continue", command=use_settings).grid(row=0, column=1)

    window.protocol("WM_DELETE_WINDOW", cancel)
    window.bind("<Return>", lambda _event: use_settings())
    window.bind("<Escape>", lambda _event: cancel())
    window.grab_set()
    column_box.focus_set()
    window.wait_window()
    return selected["value"]


def _choose_workflow_gui(root) -> Optional[str]:
    """Ask which oligo workflow to run."""
    import tkinter as tk
    from tkinter import ttk

    selected = {"value": None}

    window = tk.Toplevel(root)
    window.title("Oligo workflow")
    window.resizable(False, False)

    ttk.Label(window, text="Choose oligo workflow").grid(
        row=0, column=0, padx=16, pady=(16, 12), sticky="w"
    )

    def choose(value: str) -> None:
        selected["value"] = value
        window.destroy()

    ttk.Button(
        window,
        text="Transcript range",
        command=lambda: choose("transcript"),
        width=28,
    ).grid(row=1, column=0, padx=16, pady=(0, 8), sticky="ew")
    ttk.Button(
        window,
        text="Transcript range table",
        command=lambda: choose("transcript_table"),
        width=28,
    ).grid(row=2, column=0, padx=16, pady=(0, 8), sticky="ew")
    ttk.Button(
        window,
        text="Batch antisense table",
        command=lambda: choose("table"),
        width=28,
    ).grid(row=3, column=0, padx=16, pady=(0, 16), sticky="ew")

    window.protocol("WM_DELETE_WINDOW", window.destroy)
    window.bind("<Escape>", lambda _event: window.destroy())
    window.grab_set()
    window.wait_window()
    return selected["value"]


def _choose_transcript_workflow_gui(root) -> Optional[str]:
    """Ask which transcript sequence workflow to run."""
    import tkinter as tk
    from tkinter import ttk

    selected = {"value": None}

    window = tk.Toplevel(root)
    window.title("Transcript sequence workflow")
    window.resizable(False, False)

    ttk.Label(window, text="Choose transcript workflow").grid(
        row=0, column=0, padx=16, pady=(16, 12), sticky="w"
    )

    def choose(value: str) -> None:
        selected["value"] = value
        window.destroy()

    ttk.Button(
        window,
        text="Single transcript range",
        command=lambda: choose("transcript"),
        width=30,
    ).grid(row=1, column=0, padx=16, pady=(0, 8), sticky="ew")
    ttk.Button(
        window,
        text="Range table from Excel/CSV",
        command=lambda: choose("transcript_table"),
        width=30,
    ).grid(row=2, column=0, padx=16, pady=(0, 16), sticky="ew")

    window.protocol("WM_DELETE_WINDOW", window.destroy)
    window.bind("<Escape>", lambda _event: window.destroy())
    window.grab_set()
    window.wait_window()
    return selected["value"]


def _choose_transcript_range_gui(root) -> Optional[Tuple[int, int]]:
    """Ask for transcript start/end coordinates."""
    import tkinter as tk
    from tkinter import messagebox, ttk

    selected = {"value": None}

    window = tk.Toplevel(root)
    window.title("Transcript range")
    window.resizable(False, False)

    ttk.Label(window, text="Start position").grid(
        row=0, column=0, padx=16, pady=(16, 8), sticky="w"
    )
    start_var = tk.StringVar(value="1")
    ttk.Entry(window, textvariable=start_var, width=12).grid(
        row=0, column=1, padx=16, pady=(16, 8), sticky="w"
    )

    ttk.Label(window, text="End position").grid(
        row=1, column=0, padx=16, pady=8, sticky="w"
    )
    end_var = tk.StringVar(value="21")
    ttk.Entry(window, textvariable=end_var, width=12).grid(
        row=1, column=1, padx=16, pady=8, sticky="w"
    )

    buttons = ttk.Frame(window)
    buttons.grid(row=2, column=0, columnspan=2, padx=16, pady=(8, 16), sticky="e")

    def use_settings() -> None:
        try:
            start = int(start_var.get())
            end = int(end_var.get())
            if start < 1 or end < start:
                raise ValueError
        except ValueError:
            messagebox.showerror(
                "Invalid positions",
                "Start and end must be whole numbers, with end >= start.",
                parent=window,
            )
            return

        selected["value"] = (start, end)
        window.destroy()

    def cancel() -> None:
        window.destroy()

    ttk.Button(buttons, text="Cancel", command=cancel).grid(row=0, column=0, padx=(0, 8))
    ttk.Button(buttons, text="Continue", command=use_settings).grid(row=0, column=1)

    window.protocol("WM_DELETE_WINDOW", cancel)
    window.bind("<Return>", lambda _event: use_settings())
    window.bind("<Escape>", lambda _event: cancel())
    window.grab_set()
    window.wait_window()
    return selected["value"]


def _default_transcript_output_path(input_path: str, start: int, end: int) -> Path:
    """Return the default text output path for transcript SS/AS results."""
    source = Path(input_path)
    return source.with_name(f"{source.stem}_oligo_{start}_{end}.txt")


def _run_table_gui(root) -> int:
    """Run the existing batch antisense-table workflow."""
    from tkinter import filedialog, messagebox

    input_path = filedialog.askopenfilename(
        title="Select oligo CSV or Excel file",
        filetypes=[
            ("Excel or CSV files", "*.xlsx *.xls *.csv"),
            ("Excel files", "*.xlsx *.xls"),
            ("CSV files", "*.csv"),
            ("All files", "*.*"),
        ],
    )
    if not input_path:
        return 0

    data = read_table(input_path)
    if data.empty:
        messagebox.showerror(
            "No data rows",
            "The selected file has headers but no data rows to process.",
        )
        return 1

    settings = _choose_settings_gui(root, [str(column) for column in data.columns])
    if not settings:
        return 0

    column, start, end = settings
    result_path = process_table(
        input_path,
        column=column,
        start=start,
        end=end,
    )
    messagebox.showinfo(
        "Done",
        f"Appended oligo off-target columns to:\n{result_path}",
    )
    return 0


def _run_transcript_gui(root) -> int:
    """Run a transcript FASTA/plain-text file-picker workflow."""
    from tkinter import filedialog, messagebox

    input_path = filedialog.askopenfilename(
        title="Select transcript FASTA or text file",
        filetypes=[
            ("Sequence text files", "*.txt *.fa *.fasta *.fna *.ffn"),
            ("FASTA files", "*.fa *.fasta *.fna *.ffn"),
            ("Text files", "*.txt"),
            ("All files", "*.*"),
        ],
    )
    if not input_path:
        return 0

    settings = _choose_transcript_range_gui(root)
    if not settings:
        return 0

    start, end = settings
    result = transcript_file_to_oligo(input_path, start=start, end=end)
    output_text = format_transcript_oligo_result(result)
    default_output = _default_transcript_output_path(input_path, start, end)
    output_path = filedialog.asksaveasfilename(
        title="Save transcript SS/AS oligo results",
        initialdir=str(default_output.parent),
        initialfile=default_output.name,
        defaultextension=".txt",
        filetypes=[("Text files", "*.txt"), ("All files", "*.*")],
    )
    if not output_path:
        return 0

    Path(output_path).write_text(f"{output_text}\n", encoding="utf-8")
    messagebox.showinfo("Done", f"{output_text}\n\nWrote results to:\n{output_path}")
    return 0


def _run_transcript_table_gui(root) -> int:
    """Run transcript matching from an Excel/CSV start-end table."""
    from tkinter import filedialog, messagebox

    transcript_path = filedialog.askopenfilename(
        title="Select transcript FASTA or text file",
        filetypes=[
            ("Sequence text files", "*.txt *.fa *.fasta *.fna *.ffn"),
            ("FASTA files", "*.fa *.fasta *.fna *.ffn"),
            ("Text files", "*.txt"),
            ("All files", "*.*"),
        ],
    )
    if not transcript_path:
        return 0

    table_path = filedialog.askopenfilename(
        title="Select Excel table with start/end columns",
        filetypes=[
            ("Excel or CSV files", "*.xlsx *.xls *.csv"),
            ("Excel files", "*.xlsx *.xls"),
            ("CSV files", "*.csv"),
            ("All files", "*.*"),
        ],
    )
    if not table_path:
        return 0

    default_output = default_transcript_range_output_path(table_path)
    output_path = filedialog.asksaveasfilename(
        title="Save transcript range table results",
        initialdir=str(default_output.parent),
        initialfile=default_output.name,
        defaultextension=".xlsx",
        filetypes=[("Excel files", "*.xlsx"), ("CSV files", "*.csv")],
    )
    if not output_path:
        return 0

    result_path = process_transcript_range_table(
        transcript_path=transcript_path,
        table_path=table_path,
        output_path=output_path,
    )
    messagebox.showinfo("Done", f"Wrote results to:\n{result_path}")
    return 0


def run_gui() -> int:
    """Run a file-picker workflow and save processed results."""
    import tkinter as tk
    from tkinter import messagebox

    root = tk.Tk()
    root.withdraw()

    try:
        workflow = _choose_workflow_gui(root)
        if workflow == "transcript":
            return _run_transcript_gui(root)
        if workflow == "transcript_table":
            return _run_transcript_table_gui(root)
        if workflow == "table":
            return _run_table_gui(root)
        if not workflow:
            return 0
        raise ValueError(f"Unknown workflow: {workflow}")
    except Exception as error:
        messagebox.showerror("Oligo processing failed", str(error))
        return 1
    finally:
        root.destroy()


def run_offtarget_gui() -> int:
    """Run only the off-target antisense-table workflow."""
    import tkinter as tk
    from tkinter import messagebox

    root = tk.Tk()
    root.withdraw()

    try:
        return _run_table_gui(root)
    except Exception as error:
        messagebox.showerror("Off-target table processing failed", str(error))
        return 1
    finally:
        root.destroy()


def run_transcript_sequence_gui() -> int:
    """Run only transcript sequence extraction workflows."""
    import tkinter as tk
    from tkinter import messagebox

    root = tk.Tk()
    root.withdraw()

    try:
        workflow = _choose_transcript_workflow_gui(root)
        if workflow == "transcript":
            return _run_transcript_gui(root)
        if workflow == "transcript_table":
            return _run_transcript_table_gui(root)
        if not workflow:
            return 0
        raise ValueError(f"Unknown workflow: {workflow}")
    except Exception as error:
        messagebox.showerror("Transcript sequence processing failed", str(error))
        return 1
    finally:
        root.destroy()
