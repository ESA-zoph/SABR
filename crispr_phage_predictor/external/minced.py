from __future__ import annotations

import subprocess
import tempfile
from importlib.util import find_spec
from pathlib import Path
from typing import Callable

from crispr_phage_predictor.crispr import CrisprArray
from crispr_phage_predictor.io import FastaRecord
from crispr_phage_predictor.external.tools import tool_available


MINCED_COMMAND = "minced"
DICED_PACKAGE = "diced"


def minced_available() -> bool:
    return diced_available() or tool_available(MINCED_COMMAND)


def diced_available() -> bool:
    return find_spec(DICED_PACKAGE) is not None


def active_minced_backend() -> str:
    if diced_available():
        return "diced"
    if tool_available(MINCED_COMMAND):
        return "minced"
    return "unavailable"


ProgressCallback = Callable[[int, int, FastaRecord], None]


def detect_arrays_with_minced(
    records: list[FastaRecord],
    progress_callback: ProgressCallback | None = None,
) -> list[CrisprArray]:
    """Detect CRISPR arrays with a MinCED-compatible backend."""
    if diced_available():
        return detect_arrays_with_diced(records, progress_callback=progress_callback)
    return detect_arrays_with_minced_command(records, progress_callback=progress_callback)


def detect_arrays_with_diced(
    records: list[FastaRecord],
    progress_callback: ProgressCallback | None = None,
) -> list[CrisprArray]:
    """Detect CRISPR arrays with Diced, a MinCED-compatible Python package."""
    import diced

    arrays = []
    total = len(records)
    for record_index, record in enumerate(records, start=1):
        for crispr in diced.scan(record.sequence.upper()):
            repeats = [str(crispr.repeats[index]).upper() for index in range(len(crispr.repeats))]
            spacers = [str(crispr.spacers[index]).upper() for index in range(len(crispr.spacers))]
            if not repeats or not spacers:
                continue
            arrays.append(
                CrisprArray(
                    array_id="pending",
                    genome_id=record.source_file,
                    contig_id=record.record_id,
                    start=crispr.start + 1,
                    end=crispr.end,
                    repeat_consensus=_consensus_repeat(repeats),
                    spacers=spacers,
                )
            )
        if progress_callback:
            progress_callback(record_index, total, record)

    return [
        CrisprArray(
            array_id=f"{array.genome_id}|{array.contig_id}|array_{index}",
            genome_id=array.genome_id,
            contig_id=array.contig_id,
            start=array.start,
            end=array.end,
            repeat_consensus=array.repeat_consensus,
            spacers=array.spacers,
        )
        for index, array in enumerate(arrays, start=1)
    ]


def detect_arrays_with_minced_command(
    records: list[FastaRecord],
    progress_callback: ProgressCallback | None = None,
) -> list[CrisprArray]:
    """Detect CRISPR arrays with the MinCED command and convert them to local arrays."""
    if not records:
        return []

    arrays: list[CrisprArray] = []
    with tempfile.TemporaryDirectory(prefix="crispr_minced_") as temp_dir:
        work_dir = Path(temp_dir)
        input_fasta = work_dir / "bacteria.fasta"
        output_path = work_dir / "minced.txt"
        record_map = _write_records(input_fasta, records)

        subprocess.run(
            [MINCED_COMMAND, str(input_fasta), str(output_path)],
            check=True,
            capture_output=True,
            text=True,
        )
        arrays = parse_minced_output(output_path.read_text(encoding="utf-8"), record_map)

    if progress_callback:
        for record_index, record in enumerate(records, start=1):
            progress_callback(record_index, len(records), record)

    return [
        CrisprArray(
            array_id=f"{array.genome_id}|{array.contig_id}|array_{index}",
            genome_id=array.genome_id,
            contig_id=array.contig_id,
            start=array.start,
            end=array.end,
            repeat_consensus=array.repeat_consensus,
            spacers=array.spacers,
        )
        for index, array in enumerate(arrays, start=1)
    ]


def parse_minced_output(text: str, record_map: dict[str, FastaRecord]) -> list[CrisprArray]:
    arrays: list[CrisprArray] = []
    current_record: FastaRecord | None = None
    current_start: int | None = None
    current_end: int | None = None
    repeats: list[str] = []
    spacers: list[str] = []

    def flush_array() -> None:
        nonlocal current_start, current_end, repeats, spacers
        if current_record and current_start and current_end and repeats and spacers:
            arrays.append(
                CrisprArray(
                    array_id="pending",
                    genome_id=current_record.source_file,
                    contig_id=current_record.record_id,
                    start=current_start,
                    end=current_end,
                    repeat_consensus=_consensus_repeat(repeats),
                    spacers=spacers,
                )
            )
        current_start = None
        current_end = None
        repeats = []
        spacers = []

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith("Sequence '") or line.startswith("Sequence \""):
            flush_array()
            seq_id = line.split("'", 2)[1] if "'" in line else line.split('"', 2)[1]
            current_record = record_map.get(seq_id)
            continue
        if line.startswith("CRISPR"):
            flush_array()
            if "Range:" in line:
                range_text = line.split("Range:", 1)[1].strip()
                left, right = range_text.split("-", 1)
                current_start = int(left.strip())
                current_end = int(right.strip())
            continue
        if current_record and current_start and _looks_like_minced_repeat_row(line):
            parts = line.split()
            repeat = parts[1].upper()
            repeats.append(repeat)
            if len(parts) >= 3 and set(parts[2].upper()).issubset({"A", "C", "G", "T", "N"}):
                spacers.append(parts[2].upper())

    flush_array()
    return arrays


def _write_records(path: Path, records: list[FastaRecord]) -> dict[str, FastaRecord]:
    record_map = {}
    with path.open("w", encoding="utf-8") as handle:
        for index, record in enumerate(records, start=1):
            record_id = f"record_{index}"
            record_map[record_id] = record
            handle.write(f">{record_id}\n{record.sequence.upper()}\n")
    return record_map


def _looks_like_minced_repeat_row(line: str) -> bool:
    parts = line.split()
    if len(parts) < 2 or not parts[0].isdigit():
        return False
    repeat = parts[1].upper()
    return set(repeat).issubset({"A", "C", "G", "T", "N"})


def _consensus_repeat(repeats: list[str]) -> str:
    counts: dict[str, int] = {}
    for repeat in repeats:
        counts[repeat] = counts.get(repeat, 0) + 1
    return max(counts, key=counts.get)
