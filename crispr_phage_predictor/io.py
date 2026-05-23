from __future__ import annotations

import gzip
import hashlib
import re
from dataclasses import dataclass
from io import StringIO
from typing import Iterable

import pandas as pd
from Bio import SeqIO


@dataclass(frozen=True)
class FastaRecord:
    source_file: str
    record_id: str
    description: str
    sequence: str

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def gc_fraction(self) -> float:
        if not self.sequence:
            return 0.0
        seq = self.sequence.upper()
        gc_count = seq.count("G") + seq.count("C")
        return gc_count / len(seq)

    @property
    def sequence_hash(self) -> str:
        return sequence_hash(self.sequence)

    @property
    def accession(self) -> str:
        return extract_accession(self.record_id, self.description)


@dataclass(frozen=True)
class UploadedFileSummary:
    file_name: str
    size_bytes: int
    parsed_records: int
    total_bp: int
    status: str


def parse_uploaded_fastas(uploaded_files: Iterable) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    for uploaded_file in uploaded_files or []:
        text = _decode_uploaded_file(uploaded_file)
        handle = StringIO(text)
        for seq_record in SeqIO.parse(handle, "fasta"):
            records.append(
                FastaRecord(
                    source_file=uploaded_file.name,
                    record_id=seq_record.id,
                    description=seq_record.description,
                    sequence=str(seq_record.seq).upper(),
                )
            )
    return records


def summarize_uploaded_files(uploaded_files: Iterable, records: list[FastaRecord]) -> pd.DataFrame:
    rows = []
    for uploaded_file in uploaded_files or []:
        file_records = [record for record in records if record.source_file == uploaded_file.name]
        rows.append(
            {
                "file_name": uploaded_file.name,
                "size_mb": round(len(uploaded_file.getvalue()) / 1_000_000, 3),
                "parsed_records": len(file_records),
                "total_bp": sum(record.length for record in file_records),
                "status": "parsed" if file_records else "no FASTA records found",
            }
        )
    return pd.DataFrame(rows)


def summarize_records(records: list[FastaRecord]) -> pd.DataFrame:
    rows = [
        {
            "source_file": record.source_file,
            "record_id": record.record_id,
            "length_bp": record.length,
            "gc_percent": round(record.gc_fraction * 100, 2),
            "accession": record.accession,
            "sequence_hash": record.sequence_hash,
            "description": record.description,
        }
        for record in records
    ]
    return pd.DataFrame(rows)


def deduplicate_records(records: list[FastaRecord]) -> list[FastaRecord]:
    unique_records: list[FastaRecord] = []
    seen_hashes: set[str] = set()
    for record in records:
        digest = record.sequence_hash
        if digest in seen_hashes:
            continue
        seen_hashes.add(digest)
        unique_records.append(record)
    return unique_records


def summarize_duplicate_records(records: list[FastaRecord]) -> pd.DataFrame:
    groups: dict[str, list[FastaRecord]] = {}
    for record in records:
        groups.setdefault(record.sequence_hash, []).append(record)

    rows = []
    for digest, group in groups.items():
        if len(group) < 2:
            continue
        kept = group[0]
        duplicates = group[1:]
        rows.append(
            {
                "sequence_hash": digest,
                "accession": kept.accession,
                "length_bp": kept.length,
                "kept_source_file": kept.source_file,
                "kept_record_id": kept.record_id,
                "duplicate_count": len(duplicates),
                "duplicate_records": "; ".join(
                    f"{record.source_file}:{record.record_id}" for record in duplicates
                ),
            }
        )
    return pd.DataFrame(rows)


def summarize_accession_conflicts(records: list[FastaRecord]) -> pd.DataFrame:
    groups: dict[str, list[FastaRecord]] = {}
    for record in records:
        if record.accession:
            groups.setdefault(record.accession, []).append(record)

    rows = []
    for accession, group in groups.items():
        hashes = {record.sequence_hash for record in group}
        if len(hashes) <= 1:
            continue
        rows.append(
            {
                "accession": accession,
                "record_count": len(group),
                "distinct_sequence_hashes": len(hashes),
                "records": "; ".join(
                    f"{record.source_file}:{record.record_id}:{record.sequence_hash[:12]}"
                    for record in group
                ),
            }
        )
    return pd.DataFrame(rows)


def sequence_hash(sequence: str) -> str:
    normalized = "".join(str(sequence).upper().split())
    return hashlib.sha256(normalized.encode("ascii", errors="ignore")).hexdigest()


def extract_accession(record_id: str, description: str = "") -> str:
    candidates = [str(record_id or ""), str(description or "")]
    for text in candidates:
        accession = _extract_accession_from_text(text)
        if accession:
            return accession
    return ""


def _extract_accession_from_text(text: str) -> str:
    normalized = text.strip()
    if not normalized:
        return ""
    pipe_parts = normalized.split("|")
    for part in pipe_parts:
        accession = _match_accession_token(part)
        if accession:
            return accession
    for token in re.split(r"\s+", normalized):
        accession = _match_accession_token(token)
        if accession:
            return accession
    return ""


def _match_accession_token(token: str) -> str:
    cleaned = token.strip().strip(",;()[]")
    patterns = [
        r"^[A-Z]{1,2}_\d{5,9}(?:\.\d+)?$",
        r"^[A-Z]{1,4}\d{5,9}(?:\.\d+)?$",
    ]
    for pattern in patterns:
        if re.match(pattern, cleaned):
            return cleaned
    return ""


def _decode_uploaded_file(uploaded_file) -> str:
    raw_bytes = uploaded_file.getvalue()
    if uploaded_file.name.endswith(".gz") or raw_bytes.startswith(b"\x1f\x8b"):
        raw_bytes = gzip.decompress(raw_bytes)
    return raw_bytes.decode("utf-8", errors="replace")
