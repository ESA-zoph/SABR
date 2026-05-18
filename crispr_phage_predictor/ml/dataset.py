from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd


REPEAT_CAS_DATASET_COLUMNS = [
    "source",
    "genome_id",
    "organism",
    "taxonomy",
    "assembly_level",
    "contig_id",
    "array_start",
    "array_end",
    "repeat_sequence",
    "repeat_length",
    "spacer_count",
    "mean_spacer_length",
    "cas_type",
    "cas_subtype",
    "label_source",
    "label_confidence",
    "pam_rule",
]

REQUIRED_REPEAT_CAS_COLUMNS = [
    "source",
    "genome_id",
    "contig_id",
    "repeat_sequence",
    "repeat_length",
    "spacer_count",
    "cas_type",
    "cas_subtype",
    "label_source",
    "label_confidence",
]

HIGH_CONFIDENCE_LABELS = {"high", "curated", "cas_operon_supported"}


@dataclass(frozen=True)
class CasTrainingExample:
    genome_id: str
    organism: str
    repeat_sequence: str
    cas_type: str
    cas_subtype: str | None = None
    pam: str | None = None
    source: str | None = None

    @property
    def repeat_length(self) -> int:
        return len(self.repeat_sequence)


def empty_repeat_cas_training_table() -> pd.DataFrame:
    return pd.DataFrame(columns=REPEAT_CAS_DATASET_COLUMNS)


def load_repeat_cas_training_table(path: str | Path) -> pd.DataFrame:
    table = pd.read_csv(path)
    validate_repeat_cas_training_table(table)
    return table


def validate_repeat_cas_training_table(table: pd.DataFrame) -> None:
    missing = [column for column in REQUIRED_REPEAT_CAS_COLUMNS if column not in table.columns]
    if missing:
        raise ValueError(f"Training table is missing required columns: {', '.join(missing)}")

    invalid_rows: list[str] = []
    for row_number, row in table.reset_index(drop=True).iterrows():
        repeat = str(row["repeat_sequence"]).upper()
        if not repeat or not set(repeat).issubset({"A", "C", "G", "T", "N"}):
            invalid_rows.append(f"row {row_number + 2}: invalid repeat_sequence")
        if int(row["repeat_length"]) != len(repeat):
            invalid_rows.append(f"row {row_number + 2}: repeat_length does not match sequence")
        if int(row["spacer_count"]) < 0:
            invalid_rows.append(f"row {row_number + 2}: spacer_count is negative")
        if not str(row["cas_type"]).strip():
            invalid_rows.append(f"row {row_number + 2}: cas_type is empty")
        if not str(row["cas_subtype"]).strip():
            invalid_rows.append(f"row {row_number + 2}: cas_subtype is empty")

    if invalid_rows:
        raise ValueError("Invalid repeat/Cas training rows: " + "; ".join(invalid_rows))


def filter_high_confidence_labels(table: pd.DataFrame) -> pd.DataFrame:
    validate_repeat_cas_training_table(table)
    confidence = table["label_confidence"].astype(str).str.lower().str.strip()
    return table[confidence.isin(HIGH_CONFIDENCE_LABELS)].copy()
