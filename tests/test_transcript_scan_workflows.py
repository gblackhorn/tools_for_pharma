"""Contracts for explicit local workflows and remote-BLAST separation."""

from __future__ import annotations

import ast
from datetime import datetime, timezone
from pathlib import Path

from tools_for_pharma.oligo import ncbi_blast
from tools_for_pharma.oligo.ncbi_transport import BlastSubmission
from tools_for_pharma.oligo.transcript_scan import remote_blast, workflows
from tools_for_pharma.oligo.transcript_scan.models import (
    AntisenseQuery,
    AntisenseRegion,
)
from tools_for_pharma.oligo.transcript_scan.targets import PastedTargetSource


def test_facade_reexports_moved_workflow_and_remote_types() -> None:
    assert ncbi_blast.run_private_panel_scan is workflows.run_private_panel_scan
    assert ncbi_blast.BlastBatchResult is remote_blast.BlastBatchResult
    assert ncbi_blast.multi_fasta is remote_blast.multi_fasta
    assert ncbi_blast.combine_blast_csv is remote_blast.combine_blast_csv
    assert ncbi_blast.filter_blast_rows is remote_blast.filter_blast_rows


def test_explicit_local_scan_config_does_not_require_argparse() -> None:
    matches = workflows.run_local_scan(
        workflows.LocalScanConfig(
            target_sequence=">local target\nGCAU\n",
            max_mismatches=0,
        ),
        [AntisenseQuery("AS_one", "AUGC")],
        [AntisenseRegion("full")],
    )

    assert len(matches) == 1
    assert matches[0].transcript_name == "local target"
    assert matches[0].mismatches == 0


def test_explicit_single_sequence_workflow_keeps_pasted_target_local() -> None:
    result = workflows.run_single_sequence_scan(
        workflows.SingleSequenceScanConfig(
            target_source=PastedTargetSource(">pasted\nGCAU\n"),
            max_mismatches=0,
            closest=2,
        ),
        AntisenseQuery("AS_one", "AUGC"),
        [AntisenseRegion("full")],
    )

    assert result.targets[0].cache_status == "pasted sequence"
    assert result.matches[0].mismatches == 0
    assert len(result.closest_matches) == 1


def test_explicit_private_panel_workflow_reuses_offline_cache(tmp_path) -> None:
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    (cache_dir / "NM_000001.1.fasta").write_text(
        ">NM_000001.1 cached target\nGCAU\n",
        encoding="utf-8",
    )

    result = workflows.run_private_panel_workflow(
        workflows.PrivatePanelWorkflowConfig(
            accessions=("NM_000001.1",),
            cache_dir=cache_dir,
            offline=True,
            max_mismatches=0,
        ),
        [AntisenseQuery("AS_one", "AUGC")],
        [AntisenseRegion("full")],
    )

    assert result.targets[0].cache_status == "cache"
    assert result.matches[0].mismatches == 0


def test_remote_blast_config_preserves_batches_rids_and_privacy_warning(tmp_path) -> None:
    clients = []

    class FakeBlastClient:
        def __init__(self, **kwargs) -> None:
            self.kwargs = kwargs
            self.submissions = []
            self.waits = []
            self.fetches = []
            clients.append(self)

        def submit_blastn(self, **kwargs) -> BlastSubmission:
            self.submissions.append(kwargs)
            return BlastSubmission(f"RID_{len(self.submissions)}", None)

        def wait_for_result(self, rid, **kwargs) -> None:
            self.waits.append((rid, kwargs))

        def fetch_csv(self, rid, **kwargs) -> str:
            self.fetches.append((rid, kwargs))
            return f"query,{rid},100,4,0,0,1,4,1,4,0,8\n"

    messages = []
    rid_log = tmp_path / "rid_log.csv"
    results = remote_blast.run_blast_batches(
        remote_blast.RemoteBlastConfig(
            email="test@example.com",
            max_batch_bases=4,
            request_seconds=0,
            poll_seconds=0,
            timeout_seconds=9,
            hitlist_size=7,
            rid_log=rid_log,
        ),
        [AntisenseQuery("A B", "AUGC"), AntisenseQuery("A_B", "CCGA")],
        client_factory=FakeBlastClient,
        sleeper=lambda _seconds: None,
        message_callback=messages.append,
    )

    assert len(results) == 2
    assert [result.submission.rid for result in results] == ["RID_1", "RID_2"]
    assert [result.queries[0].blast_query_id for result in results] == ["A_B", "A_B_2"]
    assert "Privacy warning" in messages[0]
    assert "transmit 2 oligo sequence(s)" in messages[0]
    assert clients[0].kwargs["email"] == "test@example.com"
    assert clients[0].waits == [
        ("RID_1", {"poll_seconds": 75, "timeout_seconds": 9}),
        ("RID_2", {"poll_seconds": 75, "timeout_seconds": 9}),
    ]
    assert rid_log.read_text(encoding="utf-8").count("RID_") == 2


def test_rid_log_clock_is_injectable(tmp_path) -> None:
    path = tmp_path / "rid.csv"
    result = remote_blast.BlastBatchResult(
        batch_index=1,
        submission=BlastSubmission("RID_CLOCK", 12),
        queries=(AntisenseQuery("AS_one", "AUGC"),),
        csv_text="",
    )

    remote_blast.append_rid_log(
        path,
        result,
        now=lambda: datetime(2026, 8, 12, tzinfo=timezone.utc),
    )

    text = path.read_text(encoding="utf-8")
    assert "RID_CLOCK" in text
    assert "2026-08-12T00:00:00+00:00" in text


def test_remote_blast_filtering_is_kept_with_remote_service() -> None:
    rows = [
        {"mismatches": "1", "gap_opens": "0", "alignment_fraction": 1.0},
        {"mismatches": "4", "gap_opens": "0", "alignment_fraction": 1.0},
        {"mismatches": "0", "gap_opens": "1", "alignment_fraction": 1.0},
    ]

    assert remote_blast.filter_blast_rows(rows, 3, 0, 0.8) == rows[:1]


def test_local_workflows_do_not_import_remote_blast_client() -> None:
    source = ast.parse(Path(workflows.__file__).read_text(encoding="utf-8"))
    imported_modules = {
        alias.name
        for node in ast.walk(source)
        if isinstance(node, ast.Import)
        for alias in node.names
    }
    imported_from = {
        node.module
        for node in ast.walk(source)
        if isinstance(node, ast.ImportFrom) and node.module
    }
    imported_names = {
        alias.name
        for node in ast.walk(source)
        if isinstance(node, ast.ImportFrom)
        for alias in node.names
    }

    assert "tools_for_pharma.oligo.transcript_scan.remote_blast" not in (
        imported_modules | imported_from
    )
    assert "NcbiBlastClient" not in imported_names
    assert "argparse" not in imported_modules
