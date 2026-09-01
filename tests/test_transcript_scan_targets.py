from __future__ import annotations

import ast
from pathlib import Path

import pytest

from tools_for_pharma.oligo import ncbi_blast
from tools_for_pharma.oligo.transcript_scan import targets


def test_legacy_facade_reexports_target_operations() -> None:
    reexports = {
        "transcript_cache_path": targets.transcript_cache_path,
        "validate_single_transcript_record": (
            targets.validate_single_transcript_record
        ),
        "format_cached_transcript_fasta": (
            targets.format_cached_transcript_fasta
        ),
        "transcript_target_from_fasta": targets.transcript_target_from_fasta,
        "retrieve_transcript_targets": targets.retrieve_transcript_targets,
        "fetch_transcript_fasta": targets.fetch_transcript_fasta,
        "read_transcript_input": targets.read_transcript_input,
        "prepare_pasted_transcript_sequence": (
            targets.prepare_pasted_transcript_sequence
        ),
    }

    assert {
        name: getattr(ncbi_blast, name) is implementation
        for name, implementation in reexports.items()
    } == {name: True for name in reexports}


def test_target_sources_are_explicit_and_pasted_content_is_not_represented(
    tmp_path,
) -> None:
    accession = targets.transcript_target_source(
        accession="NM_000001.1",
        cache_dir=tmp_path / "cache",
        refresh=True,
    )
    pasted = targets.transcript_target_source(
        transcript_sequence="PRIVATE_TRANSCRIPT_SEQUENCE",
        target_name="manual target",
    )
    local_file = targets.transcript_target_source(
        transcript_file=tmp_path / "target.fasta"
    )

    assert accession == targets.AccessionTargetSource(
        "NM_000001.1",
        tmp_path / "cache",
        refresh=True,
    )
    assert pasted == targets.PastedTargetSource(
        "PRIVATE_TRANSCRIPT_SEQUENCE",
        "manual target",
    )
    assert "PRIVATE_TRANSCRIPT_SEQUENCE" not in repr(pasted)
    assert local_file == targets.LocalFileTargetSource(tmp_path / "target.fasta")


def test_target_source_requires_exactly_one_source(tmp_path) -> None:
    with pytest.raises(ValueError, match="Provide exactly one"):
        targets.transcript_target_source()
    with pytest.raises(ValueError, match="Provide exactly one"):
        targets.transcript_target_source(
            transcript_sequence="AUGC",
            transcript_file=tmp_path / "target.fasta",
        )


def test_pasted_and_local_file_targets_do_not_create_cache(tmp_path) -> None:
    cache_dir = tmp_path / "cache"
    target_file = tmp_path / "manual.fasta"
    target_file.write_text(">local target\nGCAU\n", encoding="utf-8")

    pasted_result = targets.local_transcript_target(
        targets.PastedTargetSource("CCCAUGCUUU", "manual target")
    )
    file_result = targets.local_transcript_target(
        targets.LocalFileTargetSource(target_file)
    )

    assert pasted_result.transcript_name == "manual target"
    assert pasted_result.cache_status == "pasted sequence"
    assert pasted_result.cache_path == ""
    assert file_result.transcript_name == "local target"
    assert file_result.cache_status == "local file"
    assert file_result.cache_path == str(target_file)
    assert not cache_dir.exists()


def test_explicit_sources_preserve_single_record_validation(tmp_path) -> None:
    target_file = tmp_path / "multiple.fasta"
    target_file.write_text(">one\nAUGC\n>two\nGCAU\n", encoding="utf-8")

    with pytest.raises(ValueError, match="contains 2 FASTA records"):
        targets.read_transcript_source(
            targets.LocalFileTargetSource(target_file)
        )
    with pytest.raises(ValueError, match="contains 2 FASTA records"):
        targets.read_transcript_source(
            targets.PastedTargetSource(">one\nAUGC\n>two\nGCAU\n")
        )


def test_accession_cache_miss_downloads_then_reuses_without_network(tmp_path) -> None:
    class RecordingClient:
        email = "test@example.com"

        def __init__(self) -> None:
            self.requests: list[dict[str, object]] = []

        def get_text(self, _url: str, params: dict[str, object]) -> str:
            self.requests.append(params)
            return ">NM_000001.1 downloaded target\nGCAU\n"

    class UnexpectedClient:
        email = "test@example.com"

        def get_text(self, *_args, **_kwargs) -> str:
            raise AssertionError("Cached accession must not be downloaded again.")

    cache_dir = tmp_path / "cache"
    source = targets.AccessionTargetSource("NM_000001.1", cache_dir)
    client = RecordingClient()

    first = targets.read_transcript_source(
        source,
        email="test@example.com",
        client=client,
    )
    second = targets.read_transcript_source(
        source,
        email="test@example.com",
        client=UnexpectedClient(),
    )

    assert first == second == ("NM_000001.1 downloaded target", "GCAU")
    assert [request["id"] for request in client.requests] == ["NM_000001.1"]
    assert all("QUERY" not in request for request in client.requests)
    assert targets.transcript_cache_path(
        cache_dir,
        "NM_000001.1",
    ).exists()


def test_explicit_accession_offline_and_refresh_policies_are_per_run(tmp_path) -> None:
    class RefreshClient:
        email = "test@example.com"

        def get_text(self, *_args, **_kwargs) -> str:
            return ">NM_000001.1 refreshed target\nAAAA\n"

    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    cache_path = targets.transcript_cache_path(cache_dir, "NM_000001.1")
    cache_path.write_text(">NM_000001.1 cached target\nGCAU\n", encoding="utf-8")

    refreshed = targets.read_transcript_source(
        targets.AccessionTargetSource(
            "NM_000001.1",
            cache_dir,
            refresh=True,
        ),
        email="test@example.com",
        client=RefreshClient(),
    )
    reused = targets.read_transcript_source(
        targets.AccessionTargetSource("NM_000001.1", cache_dir),
        email="test@example.com",
    )

    assert refreshed == ("NM_000001.1 refreshed target", "AAAA")
    assert reused == refreshed
    with pytest.raises(ValueError, match="Offline mode requires cached transcript"):
        targets.read_transcript_source(
            targets.AccessionTargetSource(
                "NM_000002.2",
                cache_dir,
                offline=True,
            ),
            email="test@example.com",
        )


def test_targets_module_avoids_interface_and_reporting_dependencies() -> None:
    tree = ast.parse(Path(targets.__file__).read_text(encoding="utf-8"))
    imported_names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported_names.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported_names.add(node.module)

    forbidden = {
        "argparse",
        "openpyxl",
        "pandas",
        "tkinter",
        "tools_for_pharma.shared.excel_utils",
    }
    assert imported_names.isdisjoint(forbidden)
