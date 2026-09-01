"""Download and compare a versioned transcript panel against one reference.

This module is intentionally independent of the single-target scanner in
``ncbi_blast.py``. Multi-record FASTA inputs are parsed as separate records so
sequence boundaries can never be joined accidentally.
"""

from __future__ import annotations

import argparse
import csv
import difflib
import json
import os
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass
from datetime import date
from pathlib import Path
from typing import Iterable, Sequence

from tools_for_pharma.sequence.nucleotides import (
    SequenceNormalizationError,
    normalize_dna as normalize_dna_sequence,
)


ENSEMBL_REST_URL = "https://rest.ensembl.org"
NCBI_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
DEFAULT_TOOL_NAME = "tools_for_pharma_transcript_panel"


@dataclass(frozen=True)
class FastaRecord:
    """One FASTA record."""

    identifier: str
    description: str
    sequence: str


@dataclass(frozen=True)
class TranscriptRequest:
    """A requested transcript and the source used to retrieve it."""

    accession: str
    source: str
    expected_name: str = ""


@dataclass(frozen=True)
class TranscriptManifestRow:
    """Retrieval and validation metadata for one mature transcript."""

    requested_accession: str
    retrieved_accession: str
    source: str
    transcript_name: str
    gene_id: str
    biotype: str
    assembly: str
    source_release: str
    sequence_length_nt: int
    cds_length_nt: int
    cds_start_1based: int | None
    cds_end_1based: int | None
    exact_version_match: bool
    retrieved_date: str
    sequence_source_url: str
    metadata_source_url: str


@dataclass(frozen=True)
class DifferenceBlock:
    """One non-identical block from a reference-anchored comparison."""

    transcript_accession: str
    block_number: int
    difference_type: str
    reference_region: str
    reference_start_1based: int | None
    reference_end_1based: int | None
    reference_boundary_after_1based: int | None
    transcript_start_1based: int | None
    transcript_end_1based: int | None
    transcript_boundary_after_1based: int | None
    reference_length_nt: int
    transcript_length_nt: int
    reference_sequence: str
    transcript_sequence: str


@dataclass(frozen=True)
class ConservedBlock:
    """One identical block shared by the reference and comparison transcript."""

    transcript_accession: str
    block_number: int
    reference_start_1based: int
    reference_end_1based: int
    transcript_start_1based: int
    transcript_end_1based: int
    length_nt: int
    reference_region: str


@dataclass(frozen=True)
class ComparisonSummary:
    """High-level comparison metrics for one transcript."""

    reference_accession: str
    transcript_accession: str
    transcript_name: str
    reference_length_nt: int
    transcript_length_nt: int
    length_delta_nt: int
    exact_full_length_match: bool
    conserved_bases_nt: int
    reference_exactly_conserved_pct: float
    transcript_exactly_conserved_pct: float
    difference_block_count: int
    insertion_block_count: int
    deletion_block_count: int
    substitution_block_count: int
    complex_replacement_block_count: int
    reference_difference_bases_nt: int
    transcript_difference_bases_nt: int


@dataclass(frozen=True)
class PanelResult:
    """All files and comparison data produced for a transcript panel."""

    output_dir: Path
    combined_fasta: Path
    manifest_json: Path
    manifest_csv: Path
    summary_csv: Path
    differences_csv: Path
    conserved_blocks_csv: Path
    workbook_data_json: Path
    manifest: tuple[TranscriptManifestRow, ...]
    summaries: tuple[ComparisonSummary, ...]
    differences: tuple[DifferenceBlock, ...]
    conserved_blocks: tuple[ConservedBlock, ...]


def normalize_dna(sequence: str) -> str:
    """Return an uppercase DNA sequence containing only supported bases."""
    try:
        return normalize_dna_sequence(
            sequence,
            allowed_ambiguity_codes="N",
            cleanup="whitespace",
        )
    except SequenceNormalizationError as error:
        if error.reason == "empty":
            raise ValueError("Sequence is empty.") from None
        raise ValueError(
            f"Unsupported sequence characters: {''.join(error.invalid_bases)}"
        ) from None


def parse_fasta_records(text: str) -> list[FastaRecord]:
    """Parse FASTA text without joining adjacent records."""

    records: list[FastaRecord] = []
    header: str | None = None
    sequence_lines: list[str] = []

    def append_record() -> None:
        if header is None:
            return
        identifier, _, description = header.partition(" ")
        records.append(
            FastaRecord(
                identifier=identifier,
                description=description.strip(),
                sequence=normalize_dna("".join(sequence_lines)),
            )
        )

    for raw_line in str(text).splitlines():
        line = raw_line.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith(">"):
            append_record()
            header = line[1:].strip()
            if not header:
                raise ValueError("FASTA header cannot be blank.")
            sequence_lines = []
        else:
            if header is None:
                raise ValueError("Sequence content appeared before the first FASTA header.")
            sequence_lines.append(line)

    append_record()
    if not records:
        raise ValueError("No FASTA records found.")
    return records


def format_fasta(record: FastaRecord, width: int = 70) -> str:
    """Format one record as FASTA."""

    header = record.identifier
    if record.description:
        header += f" {record.description}"
    sequence_lines = [
        record.sequence[index : index + width]
        for index in range(0, len(record.sequence), width)
    ]
    return f">{header}\n" + "\n".join(sequence_lines) + "\n"


def _request_text(url: str, *, accept: str, timeout_seconds: int = 60) -> str:
    request = urllib.request.Request(
        url,
        headers={
            "Accept": accept,
            "User-Agent": f"{DEFAULT_TOOL_NAME}/1.0",
        },
    )
    with urllib.request.urlopen(request, timeout=timeout_seconds) as response:
        return response.read().decode("utf-8")


def _ensembl_sequence_url(stable_id: str, sequence_type: str) -> str:
    return (
        f"{ENSEMBL_REST_URL}/sequence/id/{urllib.parse.quote(stable_id)}"
        f"?type={urllib.parse.quote(sequence_type)}"
    )


def _ensembl_lookup_url(stable_id: str) -> str:
    return f"{ENSEMBL_REST_URL}/lookup/id/{urllib.parse.quote(stable_id)}?expand=1"


def _ensembl_archive_url(stable_id: str) -> str:
    return f"{ENSEMBL_REST_URL}/archive/id/{urllib.parse.quote(stable_id)}"


def _ensembl_current_release() -> str:
    url = f"{ENSEMBL_REST_URL}/info/data"
    payload = json.loads(_request_text(url, accept="application/json"))
    releases = payload.get("releases") or []
    return str(releases[0]) if releases else ""


def _ncbi_efetch_url(
    accession: str,
    *,
    rettype: str,
    email: str | None,
) -> str:
    params = {
        "db": "nuccore",
        "id": accession,
        "rettype": rettype,
        "retmode": "text",
        "tool": DEFAULT_TOOL_NAME,
    }
    if email:
        params["email"] = email
    return f"{NCBI_EFETCH_URL}?{urllib.parse.urlencode(params)}"


def _single_fasta_record(text: str, expected_accession: str) -> FastaRecord:
    records = parse_fasta_records(text)
    if len(records) != 1:
        raise ValueError(
            f"Expected one FASTA record for {expected_accession}, received {len(records)}."
        )
    return records[0]


def _find_cds_bounds(cdna_sequence: str, cds_sequence: str) -> tuple[int, int]:
    start_index = cdna_sequence.find(cds_sequence)
    if start_index < 0:
        raise ValueError("CDS sequence was not found within the mature cDNA.")
    if cdna_sequence.find(cds_sequence, start_index + 1) >= 0:
        raise ValueError("CDS sequence occurs more than once within the mature cDNA.")
    return start_index + 1, start_index + len(cds_sequence)


def download_ensembl_transcript(
    request: TranscriptRequest,
    output_dir: Path,
    *,
    source_release: str = "",
) -> tuple[FastaRecord, FastaRecord, TranscriptManifestRow]:
    """Download one exact-version Ensembl cDNA and its CDS."""

    requested_base, separator, requested_version = request.accession.partition(".")
    if not separator or not requested_version.isdigit():
        raise ValueError(
            f"Ensembl accession must include a numeric version: {request.accession}"
        )

    cdna_url = _ensembl_sequence_url(requested_base, "cdna")
    cds_url = _ensembl_sequence_url(requested_base, "cds")
    lookup_url = _ensembl_lookup_url(requested_base)

    fasta_dir = output_dir / "fasta"
    fasta_dir.mkdir(parents=True, exist_ok=True)
    cdna_path = fasta_dir / f"{request.accession}.cdna.fasta"
    cds_path = fasta_dir / f"{request.accession}.cds.fasta"
    cdna_text = (
        cdna_path.read_text(encoding="utf-8")
        if cdna_path.exists()
        else _request_text(cdna_url, accept="text/x-fasta")
    )
    cds_text = (
        cds_path.read_text(encoding="utf-8")
        if cds_path.exists()
        else _request_text(cds_url, accept="text/x-fasta")
    )
    cdna = _single_fasta_record(cdna_text, request.accession)
    cds = _single_fasta_record(cds_text, request.accession)
    lookup = json.loads(_request_text(lookup_url, accept="application/json"))

    retrieved_accession = cdna.identifier
    if retrieved_accession != request.accession:
        raise ValueError(
            f"Requested exact Ensembl version {request.accession}, but the live "
            f"server returned {retrieved_accession}. Use the Ensembl archive that "
            "contains the requested version."
        )

    cds_start, cds_end = _find_cds_bounds(cdna.sequence, cds.sequence)
    transcript_name = str(lookup.get("display_name") or request.expected_name)
    if request.expected_name and transcript_name != request.expected_name:
        raise ValueError(
            f"Expected transcript name {request.expected_name} for {request.accession}, "
            f"but Ensembl returned {transcript_name}."
        )

    cdna_path.write_text(
        format_fasta(cdna),
        encoding="utf-8",
    )
    cds_path.write_text(
        format_fasta(cds),
        encoding="utf-8",
    )

    manifest = TranscriptManifestRow(
        requested_accession=request.accession,
        retrieved_accession=retrieved_accession,
        source="Ensembl",
        transcript_name=transcript_name,
        gene_id=str(lookup.get("Parent") or ""),
        biotype=str(lookup.get("biotype") or ""),
        assembly=str(lookup.get("assembly_name") or ""),
        source_release=source_release,
        sequence_length_nt=len(cdna.sequence),
        cds_length_nt=len(cds.sequence),
        cds_start_1based=cds_start,
        cds_end_1based=cds_end,
        exact_version_match=True,
        retrieved_date=date.today().isoformat(),
        sequence_source_url=cdna_url,
        metadata_source_url=lookup_url,
    )
    return cdna, cds, manifest


def download_ncbi_transcript(
    request: TranscriptRequest,
    output_dir: Path,
    *,
    email: str | None = None,
) -> tuple[FastaRecord, FastaRecord, TranscriptManifestRow]:
    """Download one exact-version NCBI RefSeq mature transcript and its CDS."""

    cdna_url = _ncbi_efetch_url(request.accession, rettype="fasta", email=email)
    cds_url = _ncbi_efetch_url(
        request.accession,
        rettype="fasta_cds_na",
        email=email,
    )
    fasta_dir = output_dir / "fasta"
    fasta_dir.mkdir(parents=True, exist_ok=True)
    cdna_path = fasta_dir / f"{request.accession}.cdna.fasta"
    cds_path = fasta_dir / f"{request.accession}.cds.fasta"
    cdna_text = (
        cdna_path.read_text(encoding="utf-8")
        if cdna_path.exists()
        else _request_text(cdna_url, accept="text/plain")
    )
    cds_text = (
        cds_path.read_text(encoding="utf-8")
        if cds_path.exists()
        else _request_text(cds_url, accept="text/plain")
    )
    cdna = _single_fasta_record(cdna_text, request.accession)
    cds = _single_fasta_record(cds_text, request.accession)
    retrieved_accession = cdna.identifier.split("|")[0]
    if retrieved_accession != request.accession:
        raise ValueError(
            f"Requested exact NCBI version {request.accession}, but NCBI returned "
            f"{retrieved_accession}."
        )

    cds_start, cds_end = _find_cds_bounds(cdna.sequence, cds.sequence)
    cdna_path.write_text(
        format_fasta(cdna),
        encoding="utf-8",
    )
    cds_path.write_text(
        format_fasta(cds),
        encoding="utf-8",
    )

    manifest = TranscriptManifestRow(
        requested_accession=request.accession,
        retrieved_accession=retrieved_accession,
        source="NCBI RefSeq",
        transcript_name=request.expected_name,
        gene_id="",
        biotype="protein_coding",
        assembly="GRCh38",
        source_release="RefSeq accession.version",
        sequence_length_nt=len(cdna.sequence),
        cds_length_nt=len(cds.sequence),
        cds_start_1based=cds_start,
        cds_end_1based=cds_end,
        exact_version_match=True,
        retrieved_date=date.today().isoformat(),
        sequence_source_url=cdna_url,
        metadata_source_url=(
            f"https://www.ncbi.nlm.nih.gov/nuccore/{urllib.parse.quote(request.accession)}"
        ),
    )
    return cdna, cds, manifest


def classify_reference_region(
    start_1based: int,
    end_1based: int,
    cds_start_1based: int | None,
    cds_end_1based: int | None,
) -> str:
    """Classify a reference interval relative to its CDS."""

    if cds_start_1based is None or cds_end_1based is None:
        return "unknown"
    if end_1based < cds_start_1based:
        return "5' UTR"
    if start_1based > cds_end_1based:
        return "3' UTR"
    if start_1based >= cds_start_1based and end_1based <= cds_end_1based:
        return "CDS"
    return "CDS/UTR boundary"


def _difference_type(
    tag: str,
    reference_length: int,
    transcript_length: int,
) -> str:
    if tag == "insert":
        return "insertion in transcript"
    if tag == "delete":
        return "deletion from transcript"
    if reference_length == transcript_length == 1:
        return "substitution"
    if reference_length == transcript_length:
        return "replacement"
    return "complex replacement"


def _merge_opcodes(
    opcodes: Iterable[tuple[str, int, int, int, int]],
) -> list[tuple[str, int, int, int, int]]:
    merged: list[tuple[str, int, int, int, int]] = []
    for opcode in opcodes:
        if (
            merged
            and merged[-1][0] == opcode[0]
            and merged[-1][2] == opcode[1]
            and merged[-1][4] == opcode[3]
        ):
            tag, start_a, _, start_b, _ = merged[-1]
            merged[-1] = (tag, start_a, opcode[2], start_b, opcode[4])
        else:
            merged.append(opcode)
    return merged


def _reference_opcodes(
    reference: str,
    transcript: str,
    *,
    kmer_size: int = 12,
    character_refine_limit: int = 500,
) -> list[tuple[str, int, int, int, int]]:
    """Return fast reference-anchored opcodes for similar transcript sequences.

    Character-level ``SequenceMatcher`` is accurate for short intervals but is
    slow for multi-kilobase DNA because the alphabet contains only four common
    symbols. For long inputs, overlapping k-mers provide specific anchors; the
    shorter gaps between anchors are then refined at single-base resolution.
    """

    if max(len(reference), len(transcript)) <= character_refine_limit:
        return difflib.SequenceMatcher(
            None,
            reference,
            transcript,
            autojunk=False,
        ).get_opcodes()

    if len(reference) < kmer_size or len(transcript) < kmer_size:
        tag = "replace" if reference and transcript else ("delete" if reference else "insert")
        return [(tag, 0, len(reference), 0, len(transcript))]

    reference_kmers = [
        reference[index : index + kmer_size]
        for index in range(0, len(reference) - kmer_size + 1)
    ]
    transcript_kmers = [
        transcript[index : index + kmer_size]
        for index in range(0, len(transcript) - kmer_size + 1)
    ]
    anchors = difflib.SequenceMatcher(
        None,
        reference_kmers,
        transcript_kmers,
        autojunk=True,
    ).get_matching_blocks()

    opcodes: list[tuple[str, int, int, int, int]] = []
    previous_reference_end = 0
    previous_transcript_end = 0
    for reference_start, transcript_start, token_length in anchors:
        if token_length == 0:
            continue
        base_length = token_length + kmer_size - 1
        if (
            reference_start < previous_reference_end
            or transcript_start < previous_transcript_end
        ):
            overlap = max(
                previous_reference_end - reference_start,
                previous_transcript_end - transcript_start,
            )
            reference_start += overlap
            transcript_start += overlap
            base_length -= overlap
        if base_length <= 0:
            continue

        gap_reference = reference[previous_reference_end:reference_start]
        gap_transcript = transcript[previous_transcript_end:transcript_start]
        for tag, start_a, end_a, start_b, end_b in _reference_opcodes(
            gap_reference,
            gap_transcript,
            kmer_size=kmer_size,
            character_refine_limit=character_refine_limit,
        ):
            opcodes.append(
                (
                    tag,
                    start_a + previous_reference_end,
                    end_a + previous_reference_end,
                    start_b + previous_transcript_end,
                    end_b + previous_transcript_end,
                )
            )

        opcodes.append(
            (
                "equal",
                reference_start,
                reference_start + base_length,
                transcript_start,
                transcript_start + base_length,
            )
        )
        previous_reference_end = reference_start + base_length
        previous_transcript_end = transcript_start + base_length

    tail_reference = reference[previous_reference_end:]
    tail_transcript = transcript[previous_transcript_end:]
    if tail_reference or tail_transcript:
        if (
            max(len(tail_reference), len(tail_transcript))
            > character_refine_limit
        ):
            tag = (
                "replace"
                if tail_reference and tail_transcript
                else ("delete" if tail_reference else "insert")
            )
            opcodes.append(
                (
                    tag,
                    previous_reference_end,
                    len(reference),
                    previous_transcript_end,
                    len(transcript),
                )
            )
        else:
            for tag, start_a, end_a, start_b, end_b in _reference_opcodes(
                tail_reference,
                tail_transcript,
                kmer_size=kmer_size,
                character_refine_limit=character_refine_limit,
            ):
                opcodes.append(
                    (
                        tag,
                        start_a + previous_reference_end,
                        end_a + previous_reference_end,
                        start_b + previous_transcript_end,
                        end_b + previous_transcript_end,
                    )
                )
    return _merge_opcodes(opcodes)


def compare_transcript_to_reference(
    reference: FastaRecord,
    transcript: FastaRecord,
    *,
    transcript_name: str = "",
    reference_cds_start_1based: int | None = None,
    reference_cds_end_1based: int | None = None,
) -> tuple[ComparisonSummary, list[DifferenceBlock], list[ConservedBlock]]:
    """Compare one mature transcript to a reference using exact shared blocks.

    ``SequenceMatcher`` is used with ``autojunk=False`` because DNA uses a small
    alphabet. The output is a deterministic reference-anchored difference map,
    not a local BLAST hit list.
    """

    differences: list[DifferenceBlock] = []
    conserved_blocks: list[ConservedBlock] = []
    conserved_bases = 0

    for tag, ref_start, ref_end, tx_start, tx_end in _reference_opcodes(
        reference.sequence,
        transcript.sequence,
    ):
        reference_length = ref_end - ref_start
        transcript_length = tx_end - tx_start
        if tag == "equal":
            if reference_length:
                conserved_bases += reference_length
                conserved_blocks.append(
                    ConservedBlock(
                        transcript_accession=transcript.identifier,
                        block_number=len(conserved_blocks) + 1,
                        reference_start_1based=ref_start + 1,
                        reference_end_1based=ref_end,
                        transcript_start_1based=tx_start + 1,
                        transcript_end_1based=tx_end,
                        length_nt=reference_length,
                        reference_region=classify_reference_region(
                            ref_start + 1,
                            ref_end,
                            reference_cds_start_1based,
                            reference_cds_end_1based,
                        ),
                    )
                )
            continue

        if reference_length:
            region_start = ref_start + 1
            region_end = ref_end
        else:
            # An insertion lies between reference bases ref_start and ref_start+1.
            region_start = max(1, ref_start)
            region_end = max(1, ref_start)

        differences.append(
            DifferenceBlock(
                transcript_accession=transcript.identifier,
                block_number=len(differences) + 1,
                difference_type=_difference_type(
                    tag,
                    reference_length,
                    transcript_length,
                ),
                reference_region=classify_reference_region(
                    region_start,
                    region_end,
                    reference_cds_start_1based,
                    reference_cds_end_1based,
                ),
                reference_start_1based=(ref_start + 1 if reference_length else None),
                reference_end_1based=(ref_end if reference_length else None),
                reference_boundary_after_1based=(
                    ref_start if not reference_length else None
                ),
                transcript_start_1based=(tx_start + 1 if transcript_length else None),
                transcript_end_1based=(tx_end if transcript_length else None),
                transcript_boundary_after_1based=(
                    tx_start if not transcript_length else None
                ),
                reference_length_nt=reference_length,
                transcript_length_nt=transcript_length,
                reference_sequence=reference.sequence[ref_start:ref_end],
                transcript_sequence=transcript.sequence[tx_start:tx_end],
            )
        )

    difference_types = [block.difference_type for block in differences]
    summary = ComparisonSummary(
        reference_accession=reference.identifier,
        transcript_accession=transcript.identifier,
        transcript_name=transcript_name,
        reference_length_nt=len(reference.sequence),
        transcript_length_nt=len(transcript.sequence),
        length_delta_nt=len(transcript.sequence) - len(reference.sequence),
        exact_full_length_match=reference.sequence == transcript.sequence,
        conserved_bases_nt=conserved_bases,
        reference_exactly_conserved_pct=conserved_bases / len(reference.sequence),
        transcript_exactly_conserved_pct=conserved_bases / len(transcript.sequence),
        difference_block_count=len(differences),
        insertion_block_count=difference_types.count("insertion in transcript"),
        deletion_block_count=difference_types.count("deletion from transcript"),
        substitution_block_count=difference_types.count("substitution"),
        complex_replacement_block_count=sum(
            item in {"replacement", "complex replacement"}
            for item in difference_types
        ),
        reference_difference_bases_nt=sum(
            block.reference_length_nt for block in differences
        ),
        transcript_difference_bases_nt=sum(
            block.transcript_length_nt for block in differences
        ),
    )
    return summary, differences, conserved_blocks


def _write_csv(path: Path, rows: Iterable[object], fieldnames: Sequence[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def build_transcript_panel(
    reference_request: TranscriptRequest,
    transcript_requests: Sequence[TranscriptRequest],
    output_dir: Path,
    *,
    ncbi_email: str | None = None,
) -> PanelResult:
    """Download, validate, and compare a transcript panel."""

    output_dir.mkdir(parents=True, exist_ok=True)
    reference, _, reference_manifest = download_ncbi_transcript(
        reference_request,
        output_dir,
        email=ncbi_email,
    )

    records = [reference]
    manifest = [reference_manifest]
    summaries: list[ComparisonSummary] = []
    differences: list[DifferenceBlock] = []
    conserved_blocks: list[ConservedBlock] = []

    ensembl_release = _ensembl_current_release()
    with ThreadPoolExecutor(max_workers=4) as executor:
        downloads = [
            executor.submit(
                download_ensembl_transcript,
                request,
                output_dir,
                source_release=ensembl_release,
            )
            for request in transcript_requests
        ]

    for download in downloads:
        transcript, _, transcript_manifest = download.result()
        records.append(transcript)
        manifest.append(transcript_manifest)
        summary, transcript_differences, transcript_conserved = (
            compare_transcript_to_reference(
                reference,
                transcript,
                transcript_name=transcript_manifest.transcript_name,
                reference_cds_start_1based=reference_manifest.cds_start_1based,
                reference_cds_end_1based=reference_manifest.cds_end_1based,
            )
        )
        summaries.append(summary)
        differences.extend(transcript_differences)
        conserved_blocks.extend(transcript_conserved)

    combined_fasta = output_dir / "msh3_transcript_panel.cdna.fasta"
    combined_fasta.write_text(
        "".join(format_fasta(record) for record in records),
        encoding="utf-8",
    )

    manifest_json = output_dir / "manifest.json"
    manifest_json.write_text(
        json.dumps([asdict(row) for row in manifest], indent=2),
        encoding="utf-8",
    )
    manifest_csv = output_dir / "manifest.csv"
    summary_csv = output_dir / "comparison_summary.csv"
    differences_csv = output_dir / "difference_blocks.csv"
    conserved_blocks_csv = output_dir / "conserved_blocks.csv"

    _write_csv(
        manifest_csv,
        manifest,
        tuple(TranscriptManifestRow.__dataclass_fields__),
    )
    _write_csv(
        summary_csv,
        summaries,
        tuple(ComparisonSummary.__dataclass_fields__),
    )
    _write_csv(
        differences_csv,
        differences,
        tuple(DifferenceBlock.__dataclass_fields__),
    )
    _write_csv(
        conserved_blocks_csv,
        conserved_blocks,
        tuple(ConservedBlock.__dataclass_fields__),
    )

    workbook_data_json = output_dir / "workbook_data.json"
    workbook_data_json.write_text(
        json.dumps(
            {
                "reference": asdict(reference_manifest),
                "manifest": [asdict(row) for row in manifest],
                "summaries": [asdict(row) for row in summaries],
                "differences": [asdict(row) for row in differences],
                "conserved_blocks": [asdict(row) for row in conserved_blocks],
                "sequences": [
                    {
                        "accession": record.identifier,
                        "description": record.description,
                        "length_nt": len(record.sequence),
                        "sequence": record.sequence,
                    }
                    for record in records
                ],
                "method": {
                    "reference": reference.identifier,
                    "comparison": (
                        "Reference-anchored exact-block comparison using overlapping "
                        "12-nt anchors, with short gaps refined at single-base "
                        "resolution."
                    ),
                    "interpretation": (
                        "Conserved percentages count bases in exact shared blocks. "
                        "Difference blocks are descriptive sequence differences, "
                        "not BLAST local alignments."
                    ),
                },
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    return PanelResult(
        output_dir=output_dir,
        combined_fasta=combined_fasta,
        manifest_json=manifest_json,
        manifest_csv=manifest_csv,
        summary_csv=summary_csv,
        differences_csv=differences_csv,
        conserved_blocks_csv=conserved_blocks_csv,
        workbook_data_json=workbook_data_json,
        manifest=tuple(manifest),
        summaries=tuple(summaries),
        differences=tuple(differences),
        conserved_blocks=tuple(conserved_blocks),
    )


def default_msh3_requests() -> tuple[TranscriptRequest, list[TranscriptRequest]]:
    """Return the requested MSH3 reference and Ensembl transcript panel."""

    reference = TranscriptRequest(
        accession="NM_002439.5",
        source="NCBI RefSeq",
        expected_name="MSH3 MANE Select",
    )
    transcripts = [
        TranscriptRequest("ENST00000265081.7", "Ensembl", "MSH3-201"),
        TranscriptRequest("ENST00000658259.2", "Ensembl", "MSH3-204"),
        TranscriptRequest("ENST00000667069.2", "Ensembl", "MSH3-206"),
        TranscriptRequest("ENST00000875834.2", "Ensembl", "MSH3-208"),
        TranscriptRequest("ENST00000875835.2", "Ensembl", "MSH3-209"),
        TranscriptRequest("ENST00000875836.2", "Ensembl", "MSH3-210"),
        TranscriptRequest("ENST00000933917.2", "Ensembl", "MSH3-211"),
        TranscriptRequest("ENST00000969133.2", "Ensembl", "MSH3-212"),
        TranscriptRequest("ENST00000969134.2", "Ensembl", "MSH3-213"),
        TranscriptRequest("ENST00001133246.1", "Ensembl", "MSH3-225"),
    ]
    return reference, transcripts


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Download the requested versioned MSH3 mature transcripts and compare "
            "each Ensembl cDNA with NM_002439.5."
        )
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs") / "msh3_transcript_comparison",
        help="Directory for FASTA, manifest, and comparison files.",
    )
    parser.add_argument(
        "--email",
        default=os.environ.get("NCBI_EMAIL"),
        help="Optional contact email sent to NCBI EFetch.",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    reference, transcripts = default_msh3_requests()
    result = build_transcript_panel(
        reference,
        transcripts,
        args.output_dir,
        ncbi_email=args.email,
    )
    print(f"Reference: {reference.accession}")
    print(f"Transcripts compared: {len(result.summaries)}")
    print(
        "Exact full-length matches: "
        f"{sum(row.exact_full_length_match for row in result.summaries)}"
    )
    print(f"Output directory: {result.output_dir.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
