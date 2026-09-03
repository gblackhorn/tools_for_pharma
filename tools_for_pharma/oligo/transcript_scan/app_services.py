"""Portable filesystem and settings services for the Transcript Scan app."""

from __future__ import annotations

import json
from pathlib import Path
import sys


APP_DATA_DIR_NAME = "TranscriptScanData"
GUI_SETTINGS_FILE_NAME = "settings.json"
GUI_LOG_FILE_NAME = "transcript_scan.log"


def application_base_dir() -> Path:
    """Return the executable folder, or the repository root during development."""
    if getattr(sys, "frozen", False):
        return Path(sys.executable).resolve().parent
    return Path(__file__).resolve().parents[3]


def application_data_dir() -> Path:
    """Return the writable data folder kept beside the portable application."""
    return application_base_dir() / APP_DATA_DIR_NAME


def gui_settings_path() -> Path:
    return application_data_dir() / GUI_SETTINGS_FILE_NAME


def gui_log_path() -> Path:
    return application_data_dir() / "logs" / GUI_LOG_FILE_NAME


def load_gui_settings() -> dict[str, object]:
    """Load portable per-user settings, returning an empty mapping if unavailable."""
    path = gui_settings_path()
    if not path.exists():
        return {}
    try:
        loaded = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return loaded if isinstance(loaded, dict) else {}


def save_gui_settings(settings: dict[str, object]) -> Path:
    """Atomically save portable settings beside the packaged application."""
    path = gui_settings_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = path.with_suffix(path.suffix + ".tmp")
    temporary_path.write_text(
        json.dumps(settings, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary_path.replace(path)
    return path


def shared_gui_transcript_cache_dir() -> Path:
    """Return the persistent transcript cache inside the portable data folder."""
    return application_data_dir() / "transcript_cache"
