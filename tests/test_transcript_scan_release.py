"""Release packaging contracts for the standalone Transcript Scan app."""

from __future__ import annotations

from pathlib import Path
import re

import transcript_scan_app


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEPLOYMENT_DIR = REPOSITORY_ROOT / "deployment"


def test_release_version_is_consistent() -> None:
    expected = f"{transcript_scan_app.APP_NAME} {transcript_scan_app.APP_VERSION}"
    version_text = (DEPLOYMENT_DIR / "VERSION.txt").read_text(encoding="utf-8").strip()
    readme_heading = (
        (DEPLOYMENT_DIR / "README_TRANSCRIPT_SCAN.txt")
        .read_text(encoding="utf-8")
        .splitlines()[0]
    )

    assert re.fullmatch(r"\d+\.\d+\.\d+", transcript_scan_app.APP_VERSION)
    assert version_text == expected
    assert readme_heading == expected


def test_packaged_self_test_and_release_script_share_required_files() -> None:
    package_script = (
        DEPLOYMENT_DIR / "package_transcript_scan_release.ps1"
    ).read_text(encoding="utf-8")
    build_script = (DEPLOYMENT_DIR / "build_transcript_scan.ps1").read_text(
        encoding="utf-8"
    )

    for name in transcript_scan_app.REQUIRED_DISTRIBUTION_FILES:
        assert name in package_script
        if name != "TranscriptScan.exe":
            assert name in build_script

    assert "TranscriptScanData" in package_script
    assert "settings.json" in package_script
    assert 'Extension -eq ".log"' in package_script


def test_release_documentation_covers_portability_privacy_and_workbooks() -> None:
    readme = (DEPLOYMENT_DIR / "README_TRANSCRIPT_SCAN.txt").read_text(
        encoding="utf-8"
    )
    notices = (DEPLOYMENT_DIR / "THIRD_PARTY_NOTICES.txt").read_text(
        encoding="utf-8"
    )

    assert "No Python installation is required" in readme
    assert "comparison_results sheet is second" in readme
    assert "Pasted target sequences" in readme
    assert "do not send oligo sequences to NCBI" in readme
    for dependency in ("Python", "PyInstaller", "pandas", "NumPy", "openpyxl"):
        assert dependency in notices
