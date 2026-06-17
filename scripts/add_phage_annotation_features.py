from __future__ import annotations

import argparse
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from crispr_phage_predictor.phage_annotation_features import add_phage_annotation_features


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Append GenBank-derived phage annotation features to interaction features."
    )
    parser.add_argument(
        "--features",
        type=Path,
        default=Path("data/training/phage_host_interaction_features_with_targeting.tsv"),
    )
    parser.add_argument(
        "--genbank-dir",
        type=Path,
        default=Path("data/curation/downloads/phages_genbank"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/training/phage_host_interaction_features_with_annotations.tsv"),
    )
    args = parser.parse_args()

    table = pd.read_csv(args.features, sep="\t")
    augmented = add_phage_annotation_features(table, genbank_dir=args.genbank_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    augmented.to_csv(args.output, sep="\t", index=False)
    print(f"rows\t{len(augmented)}")
    print(f"columns\t{len(augmented.columns)}")
    annotation_columns = [column for column in augmented.columns if column.startswith("phage_") and column.endswith("_count")]
    nonzero = {
        column: int((augmented[column] > 0).sum())
        for column in annotation_columns
        if int((augmented[column] > 0).sum()) > 0
    }
    print("nonzero_annotation_columns")
    for column, count in sorted(nonzero.items()):
        print(f"{column}\t{count}")


if __name__ == "__main__":
    main()
