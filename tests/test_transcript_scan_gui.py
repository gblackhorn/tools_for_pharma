"""Contracts for the extracted Transcript Scan GUI and portable services."""

from __future__ import annotations

import ast
from pathlib import Path

import tkinter.messagebox
import tkinter.simpledialog

import transcript_scan_app
from tools_for_pharma.oligo import ncbi_blast
from tools_for_pharma.oligo.transcript_scan import app_services, gui


GUI_COMPATIBILITY_EXPORTS = (
    "choose_ncbi_gui_mode",
    "choose_ncbi_gui_settings",
    "choose_sheet_gui",
    "choose_single_sequence_gui_settings",
    "default_header",
    "excel_headers",
    "gui_args",
    "prompt_and_save_ncbi_email",
    "run_gui",
    "run_single_sequence_gui",
    "run_single_sequence_scan",
    "saved_or_prompted_ncbi_email",
    "show_single_sequence_result_gui",
    "single_sequence_gui_args",
    "single_sequence_gui_draft",
)

APP_SERVICE_COMPATIBILITY_EXPORTS = (
    "application_base_dir",
    "application_data_dir",
    "gui_log_path",
    "gui_settings_path",
    "load_gui_settings",
    "save_gui_settings",
    "shared_gui_transcript_cache_dir",
)


def imported_modules(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    modules = {
        node.module
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom) and node.module
    }
    modules.update(
        alias.name
        for node in ast.walk(tree)
        if isinstance(node, ast.Import)
        for alias in node.names
    )
    return modules


def test_gui_and_portable_entry_point_do_not_depend_on_compatibility_facade() -> None:
    assert "tools_for_pharma.oligo.ncbi_blast" not in imported_modules(
        Path(gui.__file__)
    )
    assert "tools_for_pharma.oligo.ncbi_blast" not in imported_modules(
        Path(transcript_scan_app.__file__)
    )


def test_portable_services_have_no_tkinter_dependency() -> None:
    imports = imported_modules(Path(app_services.__file__))

    assert "tkinter" not in imports
    assert all(not name.startswith("tkinter.") for name in imports)


def test_facade_reexports_gui_and_portable_service_contracts() -> None:
    assert all(
        getattr(ncbi_blast, name) is getattr(gui, name)
        for name in GUI_COMPATIBILITY_EXPORTS
    )
    assert all(
        getattr(ncbi_blast, name) is getattr(app_services, name)
        for name in APP_SERVICE_COMPATIBILITY_EXPORTS
    )


def test_portable_entry_point_imports_extracted_gui_and_workflows_for_self_test() -> None:
    source = Path(transcript_scan_app.__file__).read_text(encoding="utf-8")

    assert "from tools_for_pharma.oligo.transcript_scan.gui import run_gui" in source
    assert "from tools_for_pharma.oligo.transcript_scan import gui, workflows" in source


def test_portable_settings_paths_and_invalid_settings_fallback(
    tmp_path: Path,
    monkeypatch,
) -> None:
    data_dir = tmp_path / "TranscriptScanData"
    monkeypatch.setattr(app_services, "application_data_dir", lambda: data_dir)

    assert app_services.gui_settings_path() == data_dir / "settings.json"
    assert app_services.gui_log_path() == data_dir / "logs" / "transcript_scan.log"
    assert app_services.shared_gui_transcript_cache_dir() == data_dir / "transcript_cache"
    assert app_services.load_gui_settings() == {}

    app_services.gui_settings_path().parent.mkdir(parents=True)
    app_services.gui_settings_path().write_text("not JSON", encoding="utf-8")
    assert app_services.load_gui_settings() == {}


def test_saved_email_is_reused_without_prompting(monkeypatch) -> None:
    monkeypatch.setattr(
        gui,
        "load_gui_settings",
        lambda: {"ncbi_email": " saved@example.com "},
    )

    def unexpected_prompt(*_args, **_kwargs):
        raise AssertionError("A valid saved email should not prompt again.")

    monkeypatch.setattr(gui, "prompt_and_save_ncbi_email", unexpected_prompt)

    assert gui.saved_or_prompted_ncbi_email(object()) == "saved@example.com"


def test_first_use_email_prompt_persists_only_application_settings(monkeypatch) -> None:
    saved: dict[str, object] = {}
    monkeypatch.setattr(gui, "load_gui_settings", lambda: {})
    monkeypatch.setattr(
        tkinter.simpledialog,
        "askstring",
        lambda *_args, **_kwargs: "first.user@example.com",
    )
    monkeypatch.setattr(
        tkinter.messagebox,
        "showerror",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("A valid email should not show an error.")
        ),
    )

    def record_settings(settings: dict[str, object]) -> Path:
        saved.update(settings)
        return Path("settings.json")

    monkeypatch.setattr(gui, "save_gui_settings", record_settings)

    assert gui.prompt_and_save_ncbi_email(object()) == "first.user@example.com"
    assert saved == {"ncbi_email": "first.user@example.com"}
