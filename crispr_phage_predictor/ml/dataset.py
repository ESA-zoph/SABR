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
AMBIGUOUS_SUBTYPE_MARKERS = {"", "na", "nan", "none", "unknown", "ambiguous", "hybrid"}


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


def build_training_table_from_cctyper(
    crisprs_near_cas: pd.DataFrame,
    genome_id: str,
    organism: str = "",
    taxonomy: str = "",
    assembly_level: str = "",
    pam_rules: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Normalize CCTyper CRISPR-near-Cas output into the local training schema.

    CCTyper's `crisprs_near_cas.tab` is the preferred bootstrap table because it
    contains arrays associated with a Cas operon. We keep only rows with a usable
    subtype label and a valid consensus repeat.
    """
    _require_cctyper_columns(crisprs_near_cas)
    rows = []
    for _, row in crisprs_near_cas.iterrows():
        repeat = str(row["Consensus_repeat"]).upper().strip()
        subtype = _first_label(row.get("Subtype", row.get("Prediction", "")))
        if not _is_usable_repeat(repeat) or _is_ambiguous_subtype(subtype):
            continue

        repeat_length = _coerce_int(row.get("Repeat_len"), len(repeat))
        if repeat_length != len(repeat):
            repeat_length = len(repeat)

        rows.append(
            {
                "source": "cctyper",
                "genome_id": genome_id,
                "organism": organism,
                "taxonomy": taxonomy,
                "assembly_level": assembly_level,
                "contig_id": str(row["Contig"]),
                "array_start": _coerce_int(row.get("Start"), None),
                "array_end": _coerce_int(row.get("End"), None),
                "repeat_sequence": repeat,
                "repeat_length": repeat_length,
                "spacer_count": _repeat_count_to_spacer_count(row.get("N_repeats")),
                "mean_spacer_length": _coerce_float(row.get("Spacer_len_avg"), 0.0),
                "cas_type": cas_type_from_subtype(subtype),
                "cas_subtype": subtype,
                "label_source": "nearby_cas_operon",
                "label_confidence": _cctyper_label_confidence(row),
                "pam_rule": (pam_rules or {}).get(subtype, ""),
            }
        )

    table = pd.DataFrame(rows, columns=REPEAT_CAS_DATASET_COLUMNS)
    validate_repeat_cas_training_table(table)
    return table


def load_cctyper_crisprs_near_cas(
    path: str | Path,
    genome_id: str,
    organism: str = "",
    taxonomy: str = "",
    assembly_level: str = "",
    pam_rules: dict[str, str] | None = None,
) -> pd.DataFrame:
    cctyper_table = pd.read_csv(path, sep="\t")
    return build_training_table_from_cctyper(
        crisprs_near_cas=cctyper_table,
        genome_id=genome_id,
        organism=organism,
        taxonomy=taxonomy,
        assembly_level=assembly_level,
        pam_rules=pam_rules,
    )


def cas_type_from_subtype(subtype: str) -> str:
    normalized = str(subtype).strip().upper()
    if normalized.startswith("I-"):
        return "Type I"
    if normalized.startswith("II-"):
        return "Type II"
    if normalized.startswith("III-"):
        return "Type III"
    if normalized.startswith("IV-"):
        return "Type IV"
    if normalized.startswith("V-"):
        return "Type V"
    if normalized.startswith("VI-"):
        return "Type VI"
    return "Unknown"


def _require_cctyper_columns(table: pd.DataFrame) -> None:
    required = ["Contig", "Start", "End", "Consensus_repeat", "N_repeats"]
    missing = [column for column in required if column not in table.columns]
    if "Subtype" not in table.columns and "Prediction" not in table.columns:
        missing.append("Subtype or Prediction")
    if missing:
        raise ValueError(f"CCTyper table is missing required columns: {', '.join(missing)}")


def _is_usable_repeat(repeat: str) -> bool:
    return bool(repeat) and set(repeat).issubset({"A", "C", "G", "T", "N"})


def _is_ambiguous_subtype(subtype: str) -> bool:
    normalized = str(subtype).strip().lower()
    return normalized in AMBIGUOUS_SUBTYPE_MARKERS or "/" in normalized or "," in normalized


def _first_label(value: object) -> str:
    text = str(value).strip()
    if not text:
        return ""
    for separator in [";", "|"]:
        if separator in text:
            return text.split(separator)[0].strip()
    return text


def _repeat_count_to_spacer_count(value: object) -> int:
    repeat_count = _coerce_int(value, 0) or 0
    return max(0, repeat_count - 1)


def _cctyper_label_confidence(row: pd.Series) -> str:
    trusted = str(row.get("Trusted", "")).strip().lower()
    probability = _coerce_float(row.get("Subtype_probability"), None)
    if trusted == "true" and (probability is None or probability >= 0.8):
        return "high"
    if probability is not None and probability >= 0.9:
        return "high"
    return "medium"


def _coerce_int(value: object, default: int | None) -> int | None:
    if pd.isna(value):
        return default
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def _coerce_float(value: object, default: float | None) -> float | None:
    if pd.isna(value):
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default
