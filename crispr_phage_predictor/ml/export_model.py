from __future__ import annotations

import argparse
from pathlib import Path

from crispr_phage_predictor.ml.model_artifact import (
    DEFAULT_MODEL_PATH,
    DEFAULT_TRAINING_TABLE,
    save_artifact,
    train_extra_trees_artifact,
)


def main() -> None:
    parser = argparse.ArgumentParser(description="Export the runtime Cas subtype model artifact.")
    parser.add_argument(
        "--training-table",
        type=Path,
        default=DEFAULT_TRAINING_TABLE,
        help="Repeat/Cas subtype training CSV.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_MODEL_PATH,
        help="Output joblib artifact path.",
    )
    parser.add_argument("--min-class-count", type=int, default=20)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("--n-estimators", type=int, default=400)
    args = parser.parse_args()

    artifact = train_extra_trees_artifact(
        training_table_path=args.training_table,
        min_class_count=args.min_class_count,
        random_state=args.random_state,
        n_estimators=args.n_estimators,
    )
    output = save_artifact(artifact, args.output)
    print(f"Wrote {artifact.metadata['method']} model artifact to {output}")
    print(f"Training rows: {artifact.metadata['training_rows']}")
    print(f"Classes: {', '.join(artifact.classes)}")


if __name__ == "__main__":
    main()
