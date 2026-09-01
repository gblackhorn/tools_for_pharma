"""Windows GUI entry point for the portable Transcript Scan application."""

from __future__ import annotations

import logging
from logging.handlers import RotatingFileHandler
import sys

from tools_for_pharma.oligo.ncbi_blast import (
    application_data_dir,
    gui_log_path,
    run_gui,
    shared_gui_transcript_cache_dir,
)


APP_NAME = "Transcript Scan"
APP_VERSION = "1.1.0"


class _LoggerStream:
    """Send print output to the application log in a windowed executable."""

    def __init__(self, logger: logging.Logger, level: int) -> None:
        self.logger = logger
        self.level = level

    def write(self, message: str) -> int:
        cleaned = message.rstrip()
        if cleaned:
            self.logger.log(self.level, cleaned)
        return len(message)

    def flush(self) -> None:
        return None


def configure_logging() -> None:
    log_path = gui_log_path()
    log_path.parent.mkdir(parents=True, exist_ok=True)
    handler = RotatingFileHandler(
        log_path,
        maxBytes=2_000_000,
        backupCount=3,
        encoding="utf-8",
    )
    handler.setFormatter(
        logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s")
    )
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)
    root_logger.handlers.clear()
    root_logger.addHandler(handler)

    if sys.stdout is None:
        sys.stdout = _LoggerStream(logging.getLogger("stdout"), logging.INFO)
    if sys.stderr is None:
        sys.stderr = _LoggerStream(logging.getLogger("stderr"), logging.ERROR)


def run_self_test() -> int:
    """Verify packaged imports and writable portable data paths without opening GUI."""
    import openpyxl  # noqa: F401
    import pandas
    import tkinter

    data_dir = application_data_dir()
    data_dir.mkdir(parents=True, exist_ok=True)
    shared_gui_transcript_cache_dir().mkdir(parents=True, exist_ok=True)
    workbook_path = data_dir / ".packaged_self_test.xlsx"
    try:
        pandas.DataFrame([{"status": "ok"}]).to_excel(
            workbook_path,
            index=False,
        )
        round_trip = pandas.read_excel(workbook_path)
        if round_trip.to_dict(orient="records") != [{"status": "ok"}]:
            raise RuntimeError("Packaged Excel round-trip returned unexpected data.")
    finally:
        workbook_path.unlink(missing_ok=True)

    root = tkinter.Tk()
    root.withdraw()
    root.update_idletasks()
    root.destroy()
    logging.info("Packaged self-test passed for %s %s", APP_NAME, APP_VERSION)
    return 0


def main() -> int:
    configure_logging()
    logging.info("Starting %s %s", APP_NAME, APP_VERSION)
    if "--self-test" in sys.argv[1:]:
        try:
            return run_self_test()
        except Exception:
            logging.exception("Packaged self-test failed")
            return 1

    try:
        return run_gui()
    except Exception:
        logging.exception("Unhandled GUI error")
        try:
            import tkinter as tk
            from tkinter import messagebox

            root = tk.Tk()
            root.withdraw()
            messagebox.showerror(
                "Transcript Scan failed",
                "An unexpected error occurred. Diagnostic details were saved to:\n\n"
                f"{gui_log_path()}",
            )
            root.destroy()
        except Exception:
            logging.exception("Could not display the fatal-error dialog")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
