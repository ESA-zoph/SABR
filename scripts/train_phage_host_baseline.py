from __future__ import annotations

import argparse
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from crispr_phage_predictor.phage_host_baseline import evaluate_baseline_models


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Train baseline susceptibility classifiers on SABR phage-host features."
    )
    parser.add_argument(
        "--features",
        type=Path,
        default=Path("data/training/phage_host_interaction_features.tsv"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/training/phage_host_baseline"),
    )
    parser.add_argument(
        "--splits",
        nargs="+",
        default=["row_random", "group_by_phage", "group_by_source"],
        choices=["row_random", "group_by_phage", "group_by_source"],
    )
    parser.add_argument("--test-size", type=float, default=0.3)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument(
        "--restrict-tier",
        default="",
        help="Optional dataset_tier filter, e.g. tier1_exact_pair.",
    )
    args = parser.parse_args()

    table = pd.read_csv(args.features, sep="\t")
    if args.restrict_tier:
        table = table[table["dataset_tier"] == args.restrict_tier].copy()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    summary_rows = []
    for split_strategy in args.splits:
        results = evaluate_baseline_models(
            table,
            split_strategy=split_strategy,
            test_size=args.test_size,
            random_state=args.random_state,
        )
        split_dir = args.output_dir / split_strategy
        split_dir.mkdir(parents=True, exist_ok=True)
        for result in results:
            summary_rows.append(
                {
                    "split_strategy": result.split_strategy,
                    "method": result.method,
                    "train_size": result.train_size,
                    "test_size": result.test_size,
                    "accuracy": result.accuracy,
                    "macro_f1": result.macro_f1,
                }
            )
            result.confusion.to_csv(split_dir / f"{result.method}_confusion.csv")
            (split_dir / f"{result.method}_report.txt").write_text(
                result.report,
                encoding="utf-8",
            )

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(args.output_dir / "summary.csv", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
