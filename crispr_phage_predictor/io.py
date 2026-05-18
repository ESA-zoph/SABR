from __future__ import annotations

import gzip
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
            "description": record.description,
        }
        for record in records
    ]
    return pd.DataFrame(rows)


def _decode_uploaded_file(uploaded_file) -> str:
    raw_bytes = uploaded_file.getvalue()
    if uploaded_file.name.endswith(".gz") or raw_bytes.startswith(b"\x1f\x8b"):
        raw_bytes = gzip.decompress(raw_bytes)
    return raw_bytes.decode("utf-8", errors="replace")
