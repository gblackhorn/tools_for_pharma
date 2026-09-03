"""Contracts for the extracted transcript-scan command-line interface."""

from __future__ import annotations

import ast
from pathlib import Path

import pandas as pd
import pytest

from tools_for_pharma.oligo import ncbi_blast
from tools_for_pharma.oligo.ncbi_transport import BlastSubmission
from tools_for_pharma.oligo.transcript_scan import cli
from tools_for_pharma.oligo.transcript_scan.models import AntisenseQuery
from tools_for_pharma.oligo.transcript_scan.remote_blast import BlastBatchResult
from tools_for_pharma.oligo.transcript_scan.workflows import (
    LocalScanConfig,
    PrivatePanelWorkflowConfig,
)


CLI_COMPATIBILITY_EXPORTS = (
    "args_antisense_queries",
    "build_parser",
    "local_scan_config_from_args",
    "local_transcript_target_from_args",
    "panel_accessions_from_args",
    "private_panel_cache_dir",
    "private_panel_requested",
    "read_antisense_queries",
    "read_antisense_table",
    "read_target_accession_table",
    "run_blast_batches",
    "run_local_scan",
    "run_local_scan_with_comparison",
    "run_private_panel_workflow",
    "target_accession_values",
    "validate_runtime_args",
    "write_text",
)


def test_facade_reexports_extracted_cli_contracts() -> None:
    assert {
        name: getattr(ncbi_blast, name) is getattr(cli, name)
        for name in CLI_COMPATIBILITY_EXPORTS
    } == {name: True for name in CLI_COMPATIBILITY_EXPORTS}


def test_cli_module_has_no_gui_or_compatibility_facade_dependency() -> None:
    tree = ast.parse(Path(cli.__file__).read_text(encoding="utf-8"))
    imported_modules = {
        node.module
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom) and node.module
    }
    imported_names = {
        alias.name
        for node in ast.walk(tree)
        if isinstance(node, ast.Import)
        for alias in node.names
    }

    assert "tools_for_pharma.oligo.ncbi_blast" not in imported_modules
    assert "tkinter" not in imported_names
    assert all(not name.startswith("tkinter") for name in imported_modules)


def test_facade_main_delegates_to_cli_with_gui_runner(monkeypatch) -> None:
    captured = {}

    def fake_run_cli(*, gui_runner):
        captured["gui_runner"] = gui_runner
        return 7

    monkeypatch.setattr(ncbi_blast, "run_cli", fake_run_cli)

    assert ncbi_blast.main() == 7
    assert captured["gui_runner"] is ncbi_blast.run_gui


def test_cli_gui_dispatch_uses_injected_runner() -> None:
    assert cli.main(["--gui"], gui_runner=lambda: 9) == 9


def test_cli_adapts_arguments_to_explicit_workflow_configs(tmp_path) -> None:
    local_args = cli.build_parser().parse_args(
        ["--as-sequence", "AUGC", "--target-sequence", "GCAU"]
    )
    cli.validate_runtime_args(local_args)
    local_config = cli.local_scan_config_from_args(local_args)

    panel_args = cli.build_parser().parse_args(
        [
            "--as-sequence",
            "AUGC",
            "--private-panel",
            "--target-accession",
            "NM_000001.1",
            "--cache-dir",
            str(tmp_path),
            "--offline",
        ]
    )
    cli.validate_runtime_args(panel_args)
    panel_config = cli.private_panel_config_from_args(panel_args)

    assert isinstance(local_config, LocalScanConfig)
    assert local_config.target_sequence == "GCAU"
    assert isinstance(panel_config, PrivatePanelWorkflowConfig)
    assert panel_config.accessions == ("NM_000001.1",)
    assert panel_config.offline is True


@pytest.mark.parametrize(
    "arguments",
    [
        ["--as-sequence", "AUGC", "--target-sequence", "GCAU"],
        ["--ss-sequence", "GCAU", "--target-sequence", "GCAU"],
    ],
)
def test_cli_runs_direct_as_and_ss_local_scans(arguments, capsys) -> None:
    exit_code = cli.main([*arguments, "--max-mismatches", "0", "--terminal"])

    assert exit_code == 0
    assert "Matches: 1" in capsys.readouterr().out


def test_cli_runs_local_file_target_and_writes_csv(tmp_path) -> None:
    target = tmp_path / "target.fasta"
    target.write_text(">local target\nGCAU\n", encoding="utf-8")
    output = tmp_path / "matches.csv"

    exit_code = cli.main(
        [
            "--as-sequence",
            "AUGC",
            "--target-file",
            str(target),
            "--max-mismatches",
            "0",
            "--output",
            str(output),
        ]
    )

    assert exit_code == 0
    assert "local target" in output.read_text(encoding="utf-8")


def test_cli_reuses_cached_accession_without_network(tmp_path, capsys) -> None:
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    (cache_dir / "NM_000001.1.fasta").write_text(
        ">NM_000001.1 cached target\nGCAU\n",
        encoding="utf-8",
    )

    exit_code = cli.main(
        [
            "--as-sequence",
            "AUGC",
            "--target-accession",
            "NM_000001.1",
            "--cache-dir",
            str(cache_dir),
            "--max-mismatches",
            "0",
            "--terminal",
        ]
    )

    assert exit_code == 0
    assert "Matches: 1" in capsys.readouterr().out


def test_cli_runs_table_input_and_preserves_source_columns(tmp_path) -> None:
    source = tmp_path / "queries.csv"
    source.write_text(
        "oligo_id,sequence,Pos20\nAS_one,AUGC,20\n",
        encoding="utf-8",
    )
    workbook = tmp_path / "results.xlsx"

    exit_code = cli.main(
        [
            "--as-table",
            str(source),
            "--as-column",
            "sequence",
            "--as-name-column",
            "oligo_id",
            "--target-sequence",
            "GCAU",
            "--max-mismatches",
            "0",
            "--result-workbook",
            str(workbook),
        ]
    )

    assert exit_code == 0
    input_queries = pd.read_excel(workbook, sheet_name="input_queries")
    assert input_queries.loc[0, "Pos20"] == 20


def test_cli_remote_dispatch_can_be_tested_without_network(monkeypatch, capsys) -> None:
    captured = {}

    def fake_run_blast_batches(args, queries, **_kwargs):
        captured["args"] = args
        captured["queries"] = queries
        return [
            BlastBatchResult(
                batch_index=1,
                submission=BlastSubmission("RID_FAKE", None),
                queries=(AntisenseQuery("antisense_query", "AUGC"),),
                csv_text="query,subject,100,4,0,0,1,4,1,4,0,8\n",
            )
        ]

    monkeypatch.setattr(cli, "run_blast_batches", fake_run_blast_batches)

    exit_code = cli.main(
        [
            "--as-sequence",
            "AUGC",
            "--blast-only",
            "--email",
            "test@example.com",
        ]
    )

    output = capsys.readouterr().out
    assert exit_code == 0
    assert captured["args"].blast is True
    assert captured["queries"][0].sequence_5to3 == "AUGC"
    assert "RID_FAKE" in output


def test_gui_batch_files_keep_historical_module_launch_contract() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    expected = "python -m tools_for_pharma.oligo.ncbi_blast --gui"

    for filename in ["run_ncbi_blast_gui.bat", "run_ncbi_transcript_scan_gui.bat"]:
        text = (repo_root / filename).read_text(encoding="utf-8")
        assert expected in text
