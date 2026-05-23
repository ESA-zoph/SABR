from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from crispr_phage_predictor.io import sequence_hash
from crispr_phage_predictor.ml.audit_crisprcasdb_candidates import _balanced_candidate_subset
from crispr_phage_predictor.ml.dataset import (
    REPEAT_CAS_DATASET_COLUMNS,
    load_repeat_cas_training_table,
    validate_repeat_cas_training_table,
)


def build_crisprcasdb_augmented_dataset(
    current_training_csv: str | Path,
    candidate_csv: str | Path,
    max_per_subtype: int = 500,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return current training rows plus filtered novel CRISPRCasdb candidates."""
    current = load_repeat_cas_training_table(current_training_csv).copy()
    candidates = load_repeat_cas_training_table(candidate_csv).copy()
    current["_repeat_hash"] = current["repeat_sequence"].map(sequence_hash)
    candidates["_repeat_hash"] = candidates["repeat_sequence"].map(sequence_hash)

    current_hashes = set(current["_repeat_hash"])
    conflict_hashes = _conflicting_repeat_hashes(candidates)
    novel = candidates[
        ~candidates["_repeat_hash"].isin(current_hashes | conflict_hashes)
    ].copy()
    balanced = _balanced_candidate_subset(novel, max_per_subtype=max_per_subtype)
    balanced["label_confidence"] = "computational_nearby_cas_cluster"

    augmented = pd.concat(
        [
            current[REPEAT_CAS_DATASET_COLUMNS],
            balanced[REPEAT_CAS_DATASET_COLUMNS],
        ],
        ignore_index=True,
    )
    validate_repeat_cas_training_table(augmented)
    validate_repeat_cas_training_table(balanced[REPEAT_CAS_DATASET_COLUMNS])
    return augmented, balanced[REPEAT_CAS_DATASET_COLUMNS].reset_index(drop=True)


def _conflicting_repeat_hashes(table: pd.DataFrame) -> set[str]:
    subtype_counts = table.groupby("_repeat_hash")["cas_subtype"].nunique()
    return set(subtype_counts[subtype_counts > 1].index)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build a filtered SABR training table augmented with CRISPRCasdb SQL candidates."
    )
    parser.add_argument("current_training_csv", help="Current SABR repeat/Cas training CSV.")
    parser.add_argument("candidate_csv", help="CRISPRCasdb SQL candidate repeat/Cas CSV.")
    parser.add_argument(
        "--output",
        default="data/training/repeats_cas_types_augmented_crisprcasdb_sql_balanced.csv",
        help="Output augmented training CSV path.",
    )
    parser.add_argument(
        "--candidate-output",
        default="data/training/repeats_cas_types_crisprcasdb_sql_balanced_additions.csv",
        help="Output CSV path for the filtered candidate additions only.",
    )
    parser.add_argument(
        "--max-per-subtype",
        type=int,
        default=500,
        help="Maximum novel CRISPRCasdb rows to add per subtype. Use 0 for no cap.",
    )
    args = parser.parse_args()

    augmented, additions = build_crisprcasdb_augmented_dataset(
        args.current_training_csv,
        args.candidate_csv,
        max_per_subtype=args.max_per_subtype,
    )
    output_path = Path(args.output)
    additions_path = Path(args.candidate_output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    additions_path.parent.mkdir(parents=True, exist_ok=True)
    augmented.to_csv(output_path, index=False)
    additions.to_csv(additions_path, index=False)
    print(f"Wrote {len(augmented)} augmented rows to {output_path}")
    print(f"Wrote {len(additions)} candidate addition rows to {additions_path}")


if __name__ == "__main__":
    main()
