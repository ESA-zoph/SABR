from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from crispr_phage_predictor.crispr import CrisprArray, detect_crispr_arrays
from crispr_phage_predictor.external.minced import (
    detect_arrays_with_minced,
    minced_available,
)
from crispr_phage_predictor.io import FastaRecord
from crispr_phage_predictor.ml.dataset import (
    REPEAT_CAS_DATASET_COLUMNS,
    validate_repeat_cas_training_table,
)


REQUIRED_MANIFEST_COLUMNS = [
    "fasta_path",
    "genome_id",
    "organism",
    "cas_type",
    "cas_subtype",
    "label_source",
    "label_confidence",
    "label_scope",
]


def collect_curated_training_table(
    manifest_path: str | Path,
    detector: str = "auto",
) -> pd.DataFrame:
    manifest = _read_manifest(manifest_path)
    rows = []
    for _, entry in manifest.iterrows():
        records = _load_fasta_records(Path(entry["fasta_path"]), str(entry["genome_id"]))
        arrays = _detect_arrays(records, detector=detector)
        for array in _arrays_in_label_scope(arrays, entry):
            rows.append(_training_row(entry, array))

    table = pd.DataFrame(rows, columns=REPEAT_CAS_DATASET_COLUMNS)
    validate_repeat_cas_training_table(table)
    return table


def _read_manifest(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    sep = "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","
    table = pd.read_csv(path, sep=sep).fillna("")
    missing = [column for column in REQUIRED_MANIFEST_COLUMNS if column not in table.columns]
    if missing:
        raise ValueError(f"Curated manifest is missing required columns: {', '.join(missing)}")
    return table


def _load_fasta_records(path: Path, genome_id: str) -> list[FastaRecord]:
    records = []
    for seq_record in SeqIO.parse(path, "fasta"):
        records.append(
            FastaRecord(
                source_file=genome_id,
                record_id=seq_record.id,
                description=seq_record.description,
                sequence=str(seq_record.seq).upper(),
            )
        )
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return records


def _detect_arrays(records: list[FastaRecord], detector: str) -> list[CrisprArray]:
    detector = detector.lower().strip()
    if detector not in {"auto", "minced", "internal"}:
        raise ValueError("detector must be one of: auto, minced, internal")

    if detector == "minced" or (detector == "auto" and minced_available()):
        return detect_arrays_with_minced(records)

    arrays = []
    for record in records:
        arrays.extend(
            detect_crispr_arrays(
                sequence=record.sequence,
                genome_id=record.source_file,
                contig_id=record.record_id,
            )
        )
    return arrays


def _arrays_in_label_scope(arrays: list[CrisprArray], entry: pd.Series) -> list[CrisprArray]:
    scope = str(entry["label_scope"]).strip().lower()
    if scope == "genome":
        return arrays

    contig_id = str(entry.get("contig_id", "")).strip()
    if not contig_id:
        raise ValueError("contig_id is required when label_scope is contig or array_coordinates")

    scoped = [array for array in arrays if array.contig_id == contig_id]
    if scope == "contig":
        return scoped

    if scope == "array_coordinates":
        locus_start = _required_int(entry, "locus_start")
        locus_end = _required_int(entry, "locus_end")
        return [array for array in scoped if _overlaps(array.start, array.end, locus_start, locus_end)]

    raise ValueError("label_scope must be one of: genome, contig, array_coordinates")


def _training_row(entry: pd.Series, array: CrisprArray) -> dict[str, object]:
    return {
        "source": str(entry.get("source", "curated_minced")).strip() or "curated_minced",
        "genome_id": str(entry["genome_id"]),
        "organism": str(entry.get("organism", "")),
        "taxonomy": str(entry.get("taxonomy", "")),
        "assembly_level": str(entry.get("assembly_level", "")),
        "contig_id": array.contig_id,
        "array_start": array.start,
        "array_end": array.end,
        "repeat_sequence": array.repeat_consensus,
        "repeat_length": array.repeat_length,
        "spacer_count": array.spacer_count,
        "mean_spacer_length": array.mean_spacer_length,
        "cas_type": str(entry["cas_type"]),
        "cas_subtype": str(entry["cas_subtype"]),
        "label_source": str(entry["label_source"]),
        "label_confidence": str(entry["label_confidence"]),
        "pam_rule": str(entry.get("pam_rule", "")),
    }


def _required_int(entry: pd.Series, column: str) -> int:
    value = entry.get(column, "")
    if str(value).strip() == "":
        raise ValueError(f"{column} is required when label_scope is array_coordinates")
    return int(float(value))


def _overlaps(left_start: int, left_end: int, right_start: int, right_end: int) -> bool:
    return max(left_start, right_start) <= min(left_end, right_end)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build repeat/Cas training rows from curated Cas labels and FASTA files."
    )
    parser.add_argument("manifest", help="Curated manifest CSV/TSV.")
    parser.add_argument(
        "--output",
        default="data/training/repeats_cas_types.csv",
        help="Output training CSV path.",
    )
    parser.add_argument(
        "--detector",
        choices=["auto", "minced", "internal"],
        default="auto",
        help="CRISPR array detector. Auto uses MinCED/Diced when available, otherwise internal.",
    )
    args = parser.parse_args()

    table = collect_curated_training_table(args.manifest, detector=args.detector)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_path, index=False)
    print(f"Wrote {len(table)} training rows to {output_path}")


if __name__ == "__main__":
    main()
