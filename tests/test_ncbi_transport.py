from __future__ import annotations

from urllib.error import URLError
from urllib.parse import parse_qs, urlsplit

import pytest

from tools_for_pharma.oligo.ncbi_transport import (
    BLAST_URL,
    EFETCH_URL,
    NcbiBlastClient,
    NcbiHttpClient,
    efetch_fasta_params,
)


class _Response:
    def __init__(self, text: str) -> None:
        self._payload = text.encode("utf-8")

    def __enter__(self) -> "_Response":
        return self

    def __exit__(self, *_args: object) -> None:
        return None

    def read(self) -> bytes:
        return self._payload


def test_efetch_params_identify_public_record_without_query_sequence() -> None:
    params = efetch_fasta_params(
        "NM_000041.4",
        email="user@example.com",
        tool="test_tool",
    )

    assert params == {
        "db": "nuccore",
        "id": "NM_000041.4",
        "rettype": "fasta",
        "retmode": "text",
        "tool": "test_tool",
        "email": "user@example.com",
    }
    assert "QUERY" not in params


def test_http_client_uses_injected_opener_clock_and_spacing() -> None:
    now = [100.0]
    sleeps: list[float] = []
    requests: list[tuple[str, int]] = []

    def sleeper(seconds: float) -> None:
        sleeps.append(seconds)
        now[0] += seconds

    def opener(request, *, timeout: int):
        requests.append((request.full_url, timeout))
        return _Response("ok")

    client = NcbiHttpClient(
        email="user@example.com",
        request_seconds=15,
        opener=opener,
        monotonic=lambda: now[0],
        sleeper=sleeper,
        timeout_seconds=30,
    )

    assert client.get_text(EFETCH_URL, {"id": "NM_000041.4"}) == "ok"
    now[0] += 5
    assert client.get_text(EFETCH_URL, {"id": "NM_000041.4"}) == "ok"

    assert sleeps == [10.0]
    assert [timeout for _, timeout in requests] == [30, 30]
    assert parse_qs(urlsplit(requests[0][0]).query) == {"id": ["NM_000041.4"]}


def test_http_client_translates_injected_network_failure() -> None:
    def fail(*_args, **_kwargs):
        raise URLError("offline")

    client = NcbiHttpClient(
        email="user@example.com",
        request_seconds=0,
        opener=fail,
    )

    with pytest.raises(ValueError, match="NCBI request failed: offline"):
        client.get_text(EFETCH_URL, {"id": "NM_000041.4"})


def test_remote_blast_submission_is_the_flow_that_contains_query() -> None:
    captured: list[tuple[str, dict[str, object]]] = []

    class StubBlastClient(NcbiBlastClient):
        def get_text(self, url: str, params: dict[str, object]) -> str:
            captured.append((url, params))
            return "RID = RID_TEST\nRTOE = 1\n"

    submission = StubBlastClient(
        email="user@example.com",
        request_seconds=0,
    ).submit_blastn(query_sequence="AUGC")

    assert submission.rid == "RID_TEST"
    assert captured[0][0] == BLAST_URL
    assert captured[0][1]["QUERY"] == ">oligo_query\nATGC"
