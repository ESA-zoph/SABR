from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from crispr_phage_predictor.ml.dataset import (
    REPEAT_CAS_DATASET_COLUMNS,
    cas_type_from_subtype,
    validate_repeat_cas_training_table,
)


VALID_BASES = {"A", "C", "G", "T", "N"}


def import_vink2021_candidate_repeats(
    input_csv: str | Path,
    max_per_subtype: int = 200,
) -> pd.DataFrame:
    """Import a conservative candidate repeat/type table from Vink et al. 2021.

    The source table is spacer-level. We collapse to accession + repeat +
    subtype rows and keep only rows where the repeat-inferred subtype agrees
    with the CRISPRCasdb proximity subtype metadata. These rows are not treated
    as manually curated gold labels; they are marked as computational candidates.
    """
    table = pd.read_csv(input_csv, low_memory=False)
    table = table.fillna("")
    required = ["spacers", "repeats", "accessionnrs", "subtype", "subtypesinproximity"]
    missing = [column for column in required if column not in table.columns]
    if missing:
        raise ValueError(f"Vink 2021 table is missing required columns: {', '.join(missing)}")

    table["normalized_subtype"] = table["subtype"].map(_normalize_subtype)
    table["repeat_sequence"] = table["repeats"].astype(str).str.upper().str.strip()
    table["spacer_length"] = table["spacers"].astype(str).str.len()
    filtered = table[
        (table["normalized_subtype"] != "")
        & table["repeat_sequence"].map(_is_usable_repeat)
        & table.apply(_has_matching_proximity_subtype, axis=1)
    ].copy()

    grouped = (
        filtered.groupby(
            [
                "accessionnrs",
                "repeat_sequence",
                "normalized_subtype",
                "genus",
                "family",
                "order",
                "class",
                "phylum",
                "superkingdom",
                "PAM",
            ],
            dropna=False,
        )
        .agg(spacer_count=("spacers", "nunique"), mean_spacer_length=("spacer_length", "mean"))
        .reset_index()
    )
    grouped = grouped[grouped["spacer_count"] >= 2].copy()
    grouped["repeat_length"] = grouped["repeat_sequence"].str.len()
    grouped = _balanced_sample(grouped, max_per_subtype=max_per_subtype)

    rows = []
    for _, row in grouped.iterrows():
        subtype = row["normalized_subtype"]
        taxonomy = ";".join(
            str(row[column])
            for column in ["superkingdom", "phylum", "class", "order", "family", "genus"]
            if str(row[column]).strip()
        )
        rows.append(
            {
                "source": "vink2021_crisprcasdb_candidate",
                "genome_id": str(row["accessionnrs"]),
                "organism": str(row["genus"]),
                "taxonomy": taxonomy,
                "assembly_level": "",
                "contig_id": str(row["accessionnrs"]),
                "array_start": "",
                "array_end": "",
                "repeat_sequence": row["repeat_sequence"],
                "repeat_length": int(row["repeat_length"]),
                "spacer_count": int(row["spacer_count"]),
                "mean_spacer_length": round(float(row["mean_spacer_length"]), 6),
                "cas_type": cas_type_from_subtype(subtype),
                "cas_subtype": subtype,
                "label_source": "Vink2021_CRISPRCasdb_repeat_and_proximity_subtype_agreement",
                "label_confidence": "computational_proximity",
                "pam_rule": str(row["PAM"]),
            }
        )

    output = pd.DataFrame(rows, columns=REPEAT_CAS_DATASET_COLUMNS)
    validate_repeat_cas_training_table(output)
    return output


def _normalize_subtype(value: object) -> str:
    text = str(value).strip()
    if not text.startswith("CAS-Type"):
        return ""
    suffix = text.replace("CAS-Type", "", 1).upper()
    for cas_type in ["VIII", "VII", "VI", "IV", "III", "II", "V", "I"]:
        if suffix.startswith(cas_type) and len(suffix) > len(cas_type):
            return f"{cas_type}-{suffix[len(cas_type):]}"
    return ""


def _is_usable_repeat(repeat: str) -> bool:
    return 23 <= len(repeat) <= 47 and set(repeat).issubset(VALID_BASES)


def _has_matching_proximity_subtype(row: pd.Series) -> bool:
    subtype = str(row["subtype"]).strip()
    proximity_values = {
        value.strip()
        for value in str(row["subtypesinproximity"]).replace(";", "_").split("_")
        if value.strip()
    }
    return subtype in proximity_values


def _balanced_sample(table: pd.DataFrame, max_per_subtype: int) -> pd.DataFrame:
    if max_per_subtype <= 0:
        return table
    parts = []
    ranked = table.sort_values(
        ["normalized_subtype", "spacer_count", "accessionnrs", "repeat_sequence"],
        ascending=[True, False, True, True],
    )
    for _, subtype_table in ranked.groupby("normalized_subtype", sort=True):
        parts.append(subtype_table.head(max_per_subtype))
    return pd.concat(parts, ignore_index=True) if parts else table.head(0)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Import conservative candidate repeat/Cas rows from Vink et al. 2021 supplementary data."
    )
    parser.add_argument("input_csv", help="Vink et al. 2021 Additional file 2 CSV.")
    parser.add_argument(
        "--output",
        default="data/training/repeats_cas_types_vink2021_candidate.csv",
        help="Output training CSV path.",
    )
    parser.add_argument(
        "--max-per-subtype",
        type=int,
        default=200,
        help="Maximum rows to keep per subtype. Use 0 for no cap.",
    )
    args = parser.parse_args()

    table = import_vink2021_candidate_repeats(
        args.input_csv,
        max_per_subtype=args.max_per_subtype,
    )
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_path, index=False)
    print(f"Wrote {len(table)} candidate rows to {output_path}")


if __name__ == "__main__":
    main()
