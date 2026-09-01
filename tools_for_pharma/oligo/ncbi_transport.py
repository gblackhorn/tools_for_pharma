"""NCBI request transport with explicit EFetch and remote-BLAST flows.

EFetch requests identify a public database record and never contain an oligo
query. ``NcbiBlastClient`` is deliberately separate because its ``CMD=Put``
request transmits the query sequence to NCBI.
"""

from __future__ import annotations

from dataclasses import dataclass
import re
import time
from typing import Callable
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from tools_for_pharma.sequence.fasta import FastaRecord, format_fasta
from tools_for_pharma.sequence.nucleotides import normalize_dna


BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
DEFAULT_TOOL = "tools_for_pharma_oligo"
DEFAULT_EMAIL = ""
DEFAULT_DATABASE = "refseq_rna"
DEFAULT_PROGRAM = "blastn"
DEFAULT_EXPECT = "1000"
DEFAULT_WORD_SIZE = 7
DEFAULT_MEGABLAST_WORD_SIZE = 28
DEFAULT_HITLIST_SIZE = 50
DEFAULT_POLL_SECONDS = 75
DEFAULT_REQUEST_SECONDS = 15
BLASTN_WORD_SIZES = {7, 11, 15}
MEGABLAST_WORD_SIZES = {16, 20, 24, 28, 32, 48, 64}
CONTACT_EMAIL_RE = re.compile(r"^[^@\s]+@[^@\s]+\.[^@\s]+$")


@dataclass(frozen=True)
class BlastSubmission:
    """NCBI BLAST request metadata returned after ``CMD=Put``."""

    rid: str
    rtoe_seconds: int | None


def _clean_text(value: object) -> str:
    return re.sub(r"\s+", " ", str(value)).strip()


def require_email(email: str | None) -> str:
    """Return a usable NCBI contact email or raise a clear error."""
    normalized = _clean_text(email or DEFAULT_EMAIL)
    if not CONTACT_EMAIL_RE.fullmatch(normalized):
        raise ValueError(
            "A valid contact email is required before downloading transcripts "
            "from NCBI."
        )
    return normalized


def efetch_fasta_params(
    accession: str,
    *,
    email: str,
    tool: str = DEFAULT_TOOL,
    rettype: str = "fasta",
) -> dict[str, object]:
    """Build EFetch parameters for one public NCBI nucleotide record."""
    return {
        "db": "nuccore",
        "id": accession,
        "rettype": rettype,
        "retmode": "text",
        "tool": tool,
        "email": require_email(email),
    }


def resolve_blast_word_size(word_size: int, megablast: bool) -> int:
    """Return a BLAST-mode-compatible word size or raise a clear error."""
    if megablast and word_size == DEFAULT_WORD_SIZE:
        return DEFAULT_MEGABLAST_WORD_SIZE
    allowed = MEGABLAST_WORD_SIZES if megablast else BLASTN_WORD_SIZES
    if word_size not in allowed:
        mode = "megablast" if megablast else "blastn"
        allowed_text = ", ".join(str(value) for value in sorted(allowed))
        raise ValueError(
            f"Invalid --word-size {word_size} for {mode}. "
            f"Allowed values: {allowed_text}."
        )
    return word_size


def parse_blast_field(text: str, field_name: str) -> str | None:
    """Parse BLAST API fields such as RID, RTOE, or Status."""
    pattern = re.compile(
        rf"^\s*{re.escape(field_name)}\s*=\s*(\S+)\s*$",
        re.MULTILINE,
    )
    match = pattern.search(text)
    return match.group(1) if match else None


def _single_query_fasta(sequence: str) -> str:
    return format_fasta(
        FastaRecord("oligo_query", "", normalize_dna(sequence)),
        width=80,
        trailing_newline=False,
    )


class NcbiHttpClient:
    """Small HTTP client with NCBI-friendly request spacing."""

    def __init__(
        self,
        email: str,
        tool: str = DEFAULT_TOOL,
        request_seconds: int = DEFAULT_REQUEST_SECONDS,
        *,
        opener: Callable[..., object] = urlopen,
        monotonic: Callable[[], float] = time.monotonic,
        sleeper: Callable[[float], None] = time.sleep,
        timeout_seconds: int = 120,
    ) -> None:
        self.email = require_email(email)
        self.tool = tool
        self.request_seconds = request_seconds
        self.timeout_seconds = timeout_seconds
        self._opener = opener
        self._monotonic = monotonic
        self._sleeper = sleeper
        self._last_request_time = 0.0

    def _wait_if_needed(self) -> None:
        elapsed = self._monotonic() - self._last_request_time
        wait_seconds = self.request_seconds - elapsed
        if wait_seconds > 0:
            self._sleeper(wait_seconds)

    def get_text(self, url: str, params: dict[str, object]) -> str:
        self._wait_if_needed()
        clean_params = {
            key: value
            for key, value in params.items()
            if value is not None and value != ""
        }
        query = urlencode(clean_params)
        request = Request(
            f"{url}?{query}",
            headers={"User-Agent": f"{self.tool}/1.0 ({self.email})"},
        )
        try:
            with self._opener(request, timeout=self.timeout_seconds) as response:
                text = response.read().decode("utf-8", errors="replace")
        except HTTPError as error:
            raise ValueError(f"NCBI HTTP error {error.code}: {error.reason}") from error
        except URLError as error:
            raise ValueError(f"NCBI request failed: {error.reason}") from error
        finally:
            self._last_request_time = self._monotonic()
        return text


class NcbiBlastClient(NcbiHttpClient):
    """Client whose submission method sends query sequence to remote BLAST."""

    def submit_blastn(
        self,
        query_sequence: str | None = None,
        query_fasta: str | None = None,
        database: str = DEFAULT_DATABASE,
        expect: str = DEFAULT_EXPECT,
        word_size: int = DEFAULT_WORD_SIZE,
        hitlist_size: int = DEFAULT_HITLIST_SIZE,
        megablast: bool = False,
        short_query_adjust: bool = True,
    ) -> BlastSubmission:
        word_size = resolve_blast_word_size(word_size, megablast)
        if query_fasta is None:
            if query_sequence is None:
                raise ValueError(
                    "Provide query_sequence or query_fasta for BLAST submission."
                )
            query_fasta = _single_query_fasta(query_sequence)
        params = {
            "CMD": "Put",
            "PROGRAM": DEFAULT_PROGRAM,
            "DATABASE": database,
            "QUERY": query_fasta,
            "EXPECT": expect,
            "WORD_SIZE": word_size,
            "HITLIST_SIZE": hitlist_size,
            "SHORT_QUERY_ADJUST": str(short_query_adjust).lower(),
            "FILTER": "F",
            "MEGABLAST": "on" if megablast else None,
            "tool": self.tool,
            "email": self.email,
        }
        text = self.get_text(BLAST_URL, params)
        rid = parse_blast_field(text, "RID")
        if not rid:
            raise ValueError(
                f"NCBI BLAST submission did not return an RID:\n{text[:500]}"
            )
        rtoe = parse_blast_field(text, "RTOE")
        return BlastSubmission(
            rid=rid,
            rtoe_seconds=int(rtoe) if rtoe and rtoe.isdigit() else None,
        )

    def blast_status(self, rid: str) -> str:
        text = self.get_text(
            BLAST_URL,
            {
                "CMD": "Get",
                "RID": rid,
                "FORMAT_OBJECT": "SearchInfo",
                "tool": self.tool,
                "email": self.email,
            },
        )
        status = parse_blast_field(text, "Status")
        return status or "UNKNOWN"

    def wait_for_result(
        self,
        rid: str,
        poll_seconds: int = DEFAULT_POLL_SECONDS,
        timeout_seconds: int = 1800,
    ) -> None:
        deadline = self._monotonic() + timeout_seconds
        while self._monotonic() < deadline:
            status = self.blast_status(rid)
            if status == "READY":
                return
            if status in {"FAILED", "UNKNOWN", "EXPIRED"}:
                raise ValueError(f"NCBI BLAST RID {rid} returned status {status}.")
            self._sleeper(max(poll_seconds, DEFAULT_POLL_SECONDS))
        raise TimeoutError(f"Timed out waiting for NCBI BLAST RID {rid}.")

    def fetch_csv(self, rid: str, alignments: int = DEFAULT_HITLIST_SIZE) -> str:
        return self.get_text(
            BLAST_URL,
            {
                "CMD": "Get",
                "RID": rid,
                "FORMAT_TYPE": "CSV",
                "ALIGNMENT_VIEW": "Tabular",
                "ALIGNMENTS": alignments,
                "DESCRIPTIONS": alignments,
                "tool": self.tool,
                "email": self.email,
            },
        )

    def run_blastn(
        self,
        query_sequence: str | None = None,
        query_fasta: str | None = None,
        database: str = DEFAULT_DATABASE,
        expect: str = DEFAULT_EXPECT,
        word_size: int = DEFAULT_WORD_SIZE,
        hitlist_size: int = DEFAULT_HITLIST_SIZE,
        megablast: bool = False,
        timeout_seconds: int = 1800,
    ) -> tuple[BlastSubmission, str]:
        submission = self.submit_blastn(
            query_sequence=query_sequence,
            query_fasta=query_fasta,
            database=database,
            expect=expect,
            word_size=word_size,
            hitlist_size=hitlist_size,
            megablast=megablast,
        )
        if submission.rtoe_seconds:
            self._sleeper(max(submission.rtoe_seconds, DEFAULT_REQUEST_SECONDS))
        self.wait_for_result(submission.rid, timeout_seconds=timeout_seconds)
        return submission, self.fetch_csv(submission.rid, alignments=hitlist_size)
