from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from sklearn.metrics import accuracy_score, f1_score, confusion_matrix


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Evaluate CRISPR targeting-score thresholds against infection labels."
    )
    parser.add_argument(
        "--features",
        type=Path,
        default=Path("data/training/phage_host_interaction_features_with_annotations.tsv"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/training/targeting_threshold_evaluation.tsv"),
    )
    parser.add_argument(
        "--restrict-tier",
        default="",
        help="Optional dataset_tier filter, e.g. tier1_exact_pair.",
    )
    args = parser.parse_args()

    table = pd.read_csv(args.features, sep="\t")
    table = table[table["binary_susceptibility"].isin(["susceptible", "resistant"])].copy()
    table = table[table["eop_class"] != "not_reported"].copy()
    if args.restrict_tier:
        table = table[table["dataset_tier"] == args.restrict_tier].copy()

    rows = []
    for score_column in ["crispr_targeting_score", "graded_crispr_interference_score"]:
        if score_column not in table.columns:
            continue
        for threshold in _thresholds(table[score_column]):
            predictions = [
                "resistant" if float(score) >= threshold else "susceptible"
                for score in table[score_column]
            ]
            truth = list(table["binary_susceptibility"])
            matrix = confusion_matrix(
                truth,
                predictions,
                labels=["resistant", "susceptible"],
            )
            rows.append(
                {
                    "score_column": score_column,
                    "threshold": threshold,
                    "rows": len(table),
                    "dataset_tier": args.restrict_tier or "all_modelable",
                    "accuracy": round(float(accuracy_score(truth, predictions)), 6),
                    "macro_f1": round(
                        float(
                            f1_score(
                                truth,
                                predictions,
                                labels=["resistant", "susceptible"],
                                average="macro",
                                zero_division=0,
                            )
                        ),
                        6,
                    ),
                    "true_resistant_pred_resistant": int(matrix[0, 0]),
                    "true_resistant_pred_susceptible": int(matrix[0, 1]),
                    "true_susceptible_pred_resistant": int(matrix[1, 0]),
                    "true_susceptible_pred_susceptible": int(matrix[1, 1]),
                }
            )

    results = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(args.output, sep="\t", index=False)
    if results.empty:
        print("rows\t0")
        return
    best = results.sort_values(["macro_f1", "accuracy"], ascending=False).head(10)
    print(best.to_string(index=False))


def _thresholds(scores: pd.Series) -> list[float]:
    values = sorted({round(float(value), 6) for value in scores if float(value) > 0})
    defaults = [1.0, 10.0, 25.0, 40.0, 50.0, 60.0, 75.0]
    return sorted(set(defaults + values))


if __name__ == "__main__":
    main()
