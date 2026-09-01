"""Transcript target sources, validation, retrieval, and cache handling."""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
import hashlib
from pathlib import Path
from typing import Callable

from tools_for_pharma.oligo.ncbi_transport import (
    DEFAULT_REQUEST_SECONDS,
    DEFAULT_TOOL,
    EFETCH_URL,
    NcbiHttpClient,
    efetch_fasta_params,
    require_email,
)
from tools_for_pharma.oligo.transcript import (
    fasta_or_plain_text_to_sequence,
    get_fasta_header,
)
from tools_for_pharma.oligo.transcript_accessions import (
    extract_refseq_accession_from_header,
)
from tools_for_pharma.oligo.transcript_scan.models import TranscriptTargetResult
from tools_for_pharma.oligo.transcript_scan.queries import (
    clean_text_for_id,
    sanitize_fasta_name,
)
from tools_for_pharma.sequence.fasta import FastaRecord, format_fasta
from tools_for_pharma.sequence.nucleotides import normalize_dna


@dataclass(frozen=True)
class AccessionTargetSource:
    """One public NCBI accession and its per-run cache policy."""

    accession: str
    cache_dir: Path | None = None
    offline: bool = False
    refresh: bool = False


@dataclass(frozen=True)
class PastedTargetSource:
    """One locally pasted transcript that must not be persisted automatically."""

    text: str = field(repr=False)
    name: str = ""


@dataclass(frozen=True)
class LocalFileTargetSource:
    """One local FASTA or plain-text transcript file."""

    path: Path


TranscriptTargetSource = (
    AccessionTargetSource | PastedTargetSource | LocalFileTargetSource
)


def transcript_target_source(
    *,
    transcript_sequence: str | None = None,
    transcript_file: Path | None = None,
    accession: str | None = None,
    target_name: str | None = None,
    cache_dir: Path | None = None,
    offline: bool = False,
    refresh: bool = False,
) -> TranscriptTargetSource:
    """Model exactly one accession, pasted, or local-file target source."""
    provided = [
        transcript_sequence is not None,
        transcript_file is not None,
        accession is not None,
    ]
    if sum(provided) != 1:
        raise ValueError(
            "Provide exactly one of --target-sequence, --target-file, "
            "or --target-accession."
        )
    if accession is not None:
        return AccessionTargetSource(
            accession=accession,
            cache_dir=cache_dir,
            offline=offline,
            refresh=refresh,
        )
    if transcript_file is not None:
        return LocalFileTargetSource(Path(transcript_file))
    assert transcript_sequence is not None
    return PastedTargetSource(
        text=transcript_sequence,
        name=clean_text_for_id(target_name or ""),
    )


def transcript_cache_path(cache_dir: Path, accession: str) -> Path:
    return cache_dir / f"{sanitize_fasta_name(accession)}.fasta"


def validate_single_transcript_record(text: str, source_label: str) -> None:
    """Reject multi-record FASTA instead of joining transcript records."""
    record_count = sum(
        1
        for raw_line in str(text).splitlines()
        if raw_line.lstrip().startswith(">")
    )
    if record_count > 1:
        raise ValueError(
            f"{source_label} contains {record_count} FASTA records. "
            "The current local transcript scanner accepts exactly one transcript "
            "record per target; use separate one-record FASTA files."
        )


def format_cached_transcript_fasta(
    header: str,
    sequence: str,
    width: int = 80,
) -> str:
    """Format a verified transcript while preserving its descriptive header."""
    dna = normalize_dna(sequence)
    return format_fasta(
        FastaRecord.from_header(clean_text_for_id(header), dna),
        width=width,
    )


def transcript_target_from_fasta(
    requested_accession: str,
    fasta_text: str,
    cache_path: Path,
    cache_status: str,
    retrieved_at_utc: str,
) -> TranscriptTargetResult:
    """Validate one fetched/cached exact-version transcript target."""
    validate_single_transcript_record(
        fasta_text,
        f"Transcript {requested_accession}",
    )
    header = get_fasta_header(fasta_text)
    if not header:
        raise ValueError(f"Transcript {requested_accession} FASTA header is missing.")
    retrieved_accession = extract_refseq_accession_from_header(header)
    if retrieved_accession != requested_accession:
        raise ValueError(
            f"Requested exact RefSeq version {requested_accession}, but retrieved "
            f"{retrieved_accession}."
        )
    sequence = fasta_or_plain_text_to_sequence(fasta_text)
    sequence_dna = normalize_dna(sequence)
    return TranscriptTargetResult(
        requested_accession=requested_accession,
        retrieved_accession=retrieved_accession,
        transcript_name=header,
        sequence_5to3=sequence,
        sequence_length_nt=len(sequence),
        cache_path=str(cache_path),
        cache_status=cache_status,
        exact_version_match=True,
        sequence_sha256=hashlib.sha256(sequence_dna.encode("ascii")).hexdigest(),
        retrieved_at_utc=retrieved_at_utc,
        status="ready",
    )


def retrieve_transcript_targets(
    accessions: list[str],
    *,
    email: str,
    tool: str = DEFAULT_TOOL,
    cache_dir: Path,
    offline: bool = False,
    refresh: bool = False,
    request_seconds: int = DEFAULT_REQUEST_SECONDS,
    client: NcbiHttpClient | None = None,
    progress_callback: Callable[[int, int, str, str], None] | None = None,
    cancel_check: Callable[[], bool] | None = None,
) -> list[TranscriptTargetResult]:
    """Retrieve public references without transmitting guide sequences."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    http_client = client

    results = []
    total = len(accessions)
    for index, accession in enumerate(accessions, start=1):
        if progress_callback:
            progress_callback(index - 1, total, accession, "starting")
        if cancel_check and cancel_check():
            for cancelled_index, cancelled_accession in enumerate(
                accessions[index - 1 :],
                start=index,
            ):
                cancelled_path = transcript_cache_path(
                    cache_dir,
                    cancelled_accession,
                )
                results.append(
                    TranscriptTargetResult(
                        requested_accession=cancelled_accession,
                        cache_path=str(cancelled_path),
                        cache_status=(
                            "cache" if cancelled_path.exists() else "missing"
                        ),
                        status="error",
                        error="Transcript retrieval cancelled by user.",
                    )
                )
                if progress_callback:
                    progress_callback(
                        cancelled_index,
                        total,
                        cancelled_accession,
                        "cancelled",
                    )
            break
        cache_path = transcript_cache_path(cache_dir, accession)
        cache_existed = cache_path.exists()
        try:
            if cache_existed and not refresh:
                fasta_text = cache_path.read_text(encoding="utf-8-sig")
                retrieved_at = datetime.fromtimestamp(
                    cache_path.stat().st_mtime,
                    tz=timezone.utc,
                ).isoformat()
                target = transcript_target_from_fasta(
                    accession,
                    fasta_text,
                    cache_path,
                    "cache",
                    retrieved_at,
                )
            elif offline:
                raise ValueError(
                    f"Offline mode requires cached transcript {accession} "
                    f"at {cache_path}."
                )
            else:
                if http_client is None:
                    http_client = NcbiHttpClient(
                        email=require_email(email),
                        tool=tool,
                        request_seconds=max(
                            request_seconds,
                            DEFAULT_REQUEST_SECONDS,
                        ),
                    )
                request_email = require_email(
                    getattr(http_client, "email", None) or email
                )
                fasta_text = http_client.get_text(
                    EFETCH_URL,
                    efetch_fasta_params(
                        accession,
                        email=request_email,
                        tool=tool,
                    ),
                )
                retrieved_at = datetime.now(timezone.utc).isoformat()
                target = transcript_target_from_fasta(
                    accession,
                    fasta_text,
                    cache_path,
                    "refreshed" if cache_existed else "downloaded",
                    retrieved_at,
                )
                cache_path.write_text(
                    format_cached_transcript_fasta(
                        target.transcript_name,
                        target.sequence_5to3,
                    ),
                    encoding="utf-8",
                )
            results.append(target)
        except Exception as error:
            results.append(
                TranscriptTargetResult(
                    requested_accession=accession,
                    cache_path=str(cache_path),
                    cache_status=(
                        "refresh_failed"
                        if refresh and cache_path.exists()
                        else "missing"
                        if not cache_path.exists()
                        else "invalid"
                    ),
                    status="error",
                    error=str(error),
                )
            )
        if progress_callback:
            progress_status = (
                results[-1].cache_status
                if results[-1].status == "ready"
                else results[-1].status
            )
            progress_callback(index, total, accession, progress_status)
    return results


def fetch_transcript_fasta(
    accession: str,
    email: str | None,
    tool: str = DEFAULT_TOOL,
    cache_dir: Path | None = None,
    *,
    client: NcbiHttpClient | None = None,
    offline: bool = False,
    refresh: bool = False,
) -> str:
    """Fetch one NCBI FASTA, optionally reusing or refreshing its cache."""
    cache_path = (
        transcript_cache_path(cache_dir, accession)
        if cache_dir is not None
        else None
    )
    if cache_path is not None and cache_path.exists() and not refresh:
        return cache_path.read_text(encoding="utf-8-sig")
    if offline:
        required_path = cache_path or Path(f"{accession}.fasta")
        raise ValueError(
            f"Offline mode requires cached transcript {accession} at {required_path}."
        )

    request_email = require_email(getattr(client, "email", None) or email)
    http_client = client or NcbiHttpClient(email=request_email, tool=tool)
    text = http_client.get_text(
        EFETCH_URL,
        efetch_fasta_params(
            accession,
            email=getattr(http_client, "email", None) or request_email,
            tool=tool,
        ),
    )
    if not text.lstrip().startswith(">"):
        raise ValueError(
            f"NCBI EFetch did not return FASTA for {accession}:\n{text[:500]}"
        )
    if cache_path is not None:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(text, encoding="utf-8")
    return text


def read_transcript_source(
    source: TranscriptTargetSource,
    *,
    email: str | None = None,
    tool: str = DEFAULT_TOOL,
    client: NcbiHttpClient | None = None,
) -> tuple[str, str]:
    """Return the name and normalized sequence from one explicit source."""
    if isinstance(source, AccessionTargetSource):
        fasta_text = fetch_transcript_fasta(
            source.accession,
            email=email,
            tool=tool,
            cache_dir=source.cache_dir,
            client=client,
            offline=source.offline,
            refresh=source.refresh,
        )
        validate_single_transcript_record(
            fasta_text,
            f"NCBI record {source.accession}",
        )
        return (
            get_fasta_header(fasta_text) or source.accession,
            fasta_or_plain_text_to_sequence(fasta_text),
        )
    if isinstance(source, LocalFileTargetSource):
        text = source.path.read_text(encoding="utf-8-sig")
        validate_single_transcript_record(text, str(source.path))
        return (
            get_fasta_header(text) or source.path.name,
            fasta_or_plain_text_to_sequence(text),
        )
    validate_single_transcript_record(source.text, "--target-sequence")
    return (
        source.name or get_fasta_header(source.text) or "target_transcript",
        fasta_or_plain_text_to_sequence(source.text),
    )


def read_transcript_input(
    transcript_sequence: str | None = None,
    transcript_file: Path | None = None,
    accession: str | None = None,
    email: str | None = None,
    tool: str = DEFAULT_TOOL,
    cache_dir: Path | None = None,
) -> tuple[str, str]:
    """Compatibility adapter for one accession, pasted, or file target."""
    source = transcript_target_source(
        transcript_sequence=transcript_sequence,
        transcript_file=transcript_file,
        accession=accession,
        cache_dir=cache_dir,
    )
    return read_transcript_source(source, email=email, tool=tool)


def prepare_pasted_transcript_sequence(
    text: str,
    target_name: str | None = None,
) -> str:
    """Return one canonical FASTA record for a pasted local target."""
    source = PastedTargetSource(text=text, name=target_name or "")
    transcript_name, sequence = read_transcript_source(source)
    header = clean_text_for_id(transcript_name) or "pasted_transcript"
    return format_cached_transcript_fasta(header, sequence)


def local_transcript_target(
    source: PastedTargetSource | LocalFileTargetSource,
) -> TranscriptTargetResult:
    """Build a ready result without persisting pasted or local-file content."""
    transcript_name, sequence = read_transcript_source(source)
    sequence_dna = normalize_dna(sequence)
    if isinstance(source, LocalFileTargetSource):
        source_label = "local file"
        source_path = str(source.path)
    else:
        source_label = "pasted sequence"
        source_path = ""
    return TranscriptTargetResult(
        requested_accession=transcript_name,
        transcript_name=transcript_name,
        sequence_5to3=sequence,
        sequence_length_nt=len(sequence),
        cache_path=source_path,
        cache_status=source_label,
        exact_version_match=False,
        sequence_sha256=hashlib.sha256(sequence_dna.encode("ascii")).hexdigest(),
        status="ready",
    )
