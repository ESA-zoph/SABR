from __future__ import annotations

from pathlib import Path
import re

import pandas as pd


INTERACTION_COLUMNS = [
    "interaction_id",
    "source_key",
    "source_type",
    "pmid",
    "doi",
    "assay_type",
    "bacterium",
    "strain",
    "bacterial_accession",
    "phage",
    "phage_accession",
    "reference_host",
    "raw_eop",
    "eop_relation",
    "eop_value",
    "eop_class",
    "susceptibility_label",
    "plaque_result",
    "anti_crispr_status",
    "anti_crispr_genes",
    "crispr_interference_evidence",
    "other_defense_evidence",
    "experimental_conditions",
    "curation_status",
    "curation_confidence",
    "notes",
]

REQUIRED_INTERACTION_COLUMNS = [
    "interaction_id",
    "source_key",
    "assay_type",
    "bacterium",
    "strain",
    "phage",
    "raw_eop",
    "susceptibility_label",
    "curation_status",
    "curation_confidence",
]

ALLOWED_INTERACTION_VALUES = {
    "source_type": {"paper", "database", "supplement", "thesis", "preprint", "unknown", ""},
    "assay_type": {
        "eop",
        "plaque_assay",
        "spot_test",
        "host_range_panel",
        "adsorption",
        "growth_curve",
        "mixed",
        "unknown",
    },
    "eop_relation": {"=", "<", "<=", ">", ">=", "range", "not_reported", ""},
    "eop_class": {"high", "medium", "low", "trace", "none", "mixed", "not_reported", ""},
    "susceptibility_label": {
        "susceptible",
        "reduced_susceptibility",
        "resistant",
        "nonhost",
        "mixed",
        "unknown",
    },
    "plaque_result": {
        "clear_plaques",
        "turbid_plaques",
        "pinpoint_plaques",
        "lysis_from_without",
        "no_plaques",
        "not_reported",
        "mixed",
        "",
    },
    "anti_crispr_status": {
        "present",
        "candidate_present",
        "absent",
        "not_evaluated",
        "not_applicable",
        "unknown",
        "",
    },
    "crispr_interference_evidence": {
        "experimental",
        "spacer_pam_seed",
        "spacer_pam",
        "spacer_only",
        "none",
        "not_evaluated",
        "unknown",
        "",
    },
    "curation_status": {
        "curated",
        "candidate",
        "needs_review",
        "exclude_until_verified",
    },
    "curation_confidence": {"high", "medium", "low", "exclude"},
}


def empty_interaction_table() -> pd.DataFrame:
    return pd.DataFrame(columns=INTERACTION_COLUMNS)


def load_interaction_table(path: str | Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    validate_interaction_table(table)
    return table


def validate_interaction_table(table: pd.DataFrame) -> None:
    missing_columns = [
        column for column in REQUIRED_INTERACTION_COLUMNS if column not in table.columns
    ]
    if missing_columns:
        raise ValueError(
            "Missing required interaction columns: " + ", ".join(missing_columns)
        )

    invalid_rows: list[str] = []
    seen_ids: set[str] = set()
    for index, row in table.iterrows():
        row_number = index + 2
        interaction_id = str(row.get("interaction_id", "")).strip()
        if not interaction_id:
            invalid_rows.append(f"row {row_number}: interaction_id is empty")
        elif interaction_id in seen_ids:
            invalid_rows.append(f"row {row_number}: duplicate interaction_id {interaction_id}")
        seen_ids.add(interaction_id)

        for column in REQUIRED_INTERACTION_COLUMNS:
            if not str(row.get(column, "")).strip():
                invalid_rows.append(f"row {row_number}: {column} is empty")

        for column, allowed in ALLOWED_INTERACTION_VALUES.items():
            if column not in table.columns:
                continue
            value = str(row.get(column, "")).strip()
            if value not in allowed:
                invalid_rows.append(f"row {row_number}: invalid {column} '{value}'")

        eop_value = str(row.get("eop_value", "")).strip()
        if eop_value and _float_or_none(eop_value) is None:
            invalid_rows.append(f"row {row_number}: eop_value is not numeric")

    if invalid_rows:
        raise ValueError("Invalid interaction table:\n" + "\n".join(invalid_rows))


def normalize_interaction_table(table: pd.DataFrame) -> pd.DataFrame:
    """Fill EOP relation, numeric EOP, EOP class, and label when possible."""
    normalized = table.copy()
    for column in INTERACTION_COLUMNS:
        if column not in normalized.columns:
            normalized[column] = ""
    for index, row in normalized.iterrows():
        relation, value = parse_eop(row.get("raw_eop", ""))
        if not str(row.get("eop_relation", "")).strip():
            normalized.at[index, "eop_relation"] = relation
        if not str(row.get("eop_value", "")).strip() and value is not None:
            normalized.at[index, "eop_value"] = _format_float(value)
        if not str(row.get("eop_class", "")).strip():
            normalized.at[index, "eop_class"] = eop_class_from_value(value, relation)
        if str(row.get("susceptibility_label", "")).strip() == "unknown":
            normalized.at[index, "susceptibility_label"] = susceptibility_from_eop_class(
                normalized.at[index, "eop_class"]
            )
    return normalized[INTERACTION_COLUMNS]


def parse_eop(raw_eop: object) -> tuple[str, float | None]:
    text = _normalize_eop_text(raw_eop)
    if not text or text.lower() in {"na", "n/a", "nr", "not reported", "unknown"}:
        return "not_reported", None
    match = re.match(r"^(<=|>=|<|>|=)?\s*([0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)$", text)
    if not match:
        if re.search(r"\d\s*-\s*\d", text):
            return "range", None
        return "not_reported", None
    relation = match.group(1) or "="
    return relation, float(match.group(2))


def eop_class_from_value(value: float | None, relation: str = "=") -> str:
    if value is None:
        return "not_reported"
    if value == 0 or (relation in {"<", "<="} and value <= 1e-3):
        return "none"
    if value < 1e-3 or (relation in {"<", "<="} and value <= 1e-3):
        return "trace"
    if value < 0.1:
        return "low"
    if value < 0.5:
        return "medium"
    return "high"


def susceptibility_from_eop_class(eop_class: str) -> str:
    normalized = str(eop_class).strip()
    if normalized in {"high", "medium"}:
        return "susceptible"
    if normalized == "low":
        return "reduced_susceptibility"
    if normalized in {"trace", "none"}:
        return "resistant"
    if normalized == "mixed":
        return "mixed"
    return "unknown"


def _float_or_none(value: object) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _format_float(value: float) -> str:
    return f"{value:.8g}"


def _normalize_eop_text(raw_eop: object) -> str:
    text = str(raw_eop).replace("\xa0", " ").strip()
    text = text.replace("≤", "<=").replace("≥", ">=")
    text = re.sub(r"\s+", " ", text)
    if "±" in text:
        text = text.split("±", 1)[0].strip()
    return text
