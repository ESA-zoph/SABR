from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from crispr_phage_predictor.ml.dataset import (
    REPEAT_CAS_DATASET_COLUMNS,
    load_cctyper_crisprs_near_cas,
    validate_repeat_cas_training_table,
)


REQUIRED_MANIFEST_COLUMNS = ["cctyper_output_dir", "genome_id"]


def collect_cctyper_training_table(
    manifest_path: str | Path,
    output_path: str | Path,
) -> pd.DataFrame:
    manifest = pd.read_csv(manifest_path)
    _validate_manifest(manifest)

    tables = []
    for _, row in manifest.iterrows():
        output_dir = Path(row["cctyper_output_dir"])
        crisprs_near_cas = output_dir / "crisprs_near_cas.tab"
        if not crisprs_near_cas.exists():
            raise FileNotFoundError(f"Missing CCTyper file: {crisprs_near_cas}")

        tables.append(
            load_cctyper_crisprs_near_cas(
                crisprs_near_cas,
                genome_id=str(row["genome_id"]),
                organism=_optional_string(row, "organism"),
                taxonomy=_optional_string(row, "taxonomy"),
                assembly_level=_optional_string(row, "assembly_level"),
            )
        )

    if tables:
        training_table = pd.concat(tables, ignore_index=True)
    else:
        training_table = pd.DataFrame(columns=REPEAT_CAS_DATASET_COLUMNS)

    validate_repeat_cas_training_table(training_table)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    training_table.to_csv(output_path, index=False)
    return training_table


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Collect CCTyper output directories into repeats_cas_types.csv."
    )
    parser.add_argument("manifest_csv", type=Path)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/training/repeats_cas_types.csv"),
    )
    args = parser.parse_args()

    table = collect_cctyper_training_table(args.manifest_csv, args.output)
    print(f"Wrote {len(table)} rows to {args.output}")
    if not table.empty:
        print("Subtype counts")
        print(table["cas_subtype"].value_counts().to_string())


def _validate_manifest(manifest: pd.DataFrame) -> None:
    missing = [column for column in REQUIRED_MANIFEST_COLUMNS if column not in manifest.columns]
    if missing:
        raise ValueError(f"Manifest is missing required columns: {', '.join(missing)}")


def _optional_string(row: pd.Series, column: str) -> str:
    value = row.get(column, "")
    return "" if pd.isna(value) else str(value)


if __name__ == "__main__":
    main()
