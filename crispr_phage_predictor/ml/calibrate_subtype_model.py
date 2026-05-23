from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.ensemble import ExtraTreesClassifier

from crispr_phage_predictor.ml.dataset import load_repeat_cas_training_table
from crispr_phage_predictor.ml.features import build_repeat_feature_table, feature_columns
from crispr_phage_predictor.ml.train_classifier import _filter_min_class_count, _split_table


def calibrate_extra_trees_subtype_model(
    training_csv: str | Path,
    split_strategy: str = "group_holdout",
    group_column: str = "genus",
    test_size: float = 0.25,
    random_state: int = 42,
    min_class_count: int = 20,
    n_bins: int = 10,
) -> dict[str, pd.DataFrame]:
    """Evaluate whether ExtraTrees subtype probabilities are calibrated."""
    table = load_repeat_cas_training_table(training_csv)
    table = _filter_min_class_count(table, min_class_count=min_class_count)
    train, test = _split_table(
        table=table,
        test_size=test_size,
        random_state=random_state,
        split_strategy=split_strategy,
        group_column=group_column,
    )

    train_features = build_repeat_feature_table(train)
    test_features = build_repeat_feature_table(test)
    features = feature_columns(train_features)
    aligned_test = test_features.reindex(columns=features, fill_value=0.0)

    model = ExtraTreesClassifier(
        n_estimators=400,
        random_state=random_state,
        class_weight="balanced",
        n_jobs=-1,
    )
    model.fit(train_features[features], train_features["cas_subtype"])
    probabilities = model.predict_proba(aligned_test)
    classes = [str(label) for label in model.classes_]

    rows = []
    for row_index, row_probabilities in enumerate(probabilities):
        best_index = int(row_probabilities.argmax())
        predicted = classes[best_index]
        confidence = float(row_probabilities[best_index])
        true_label = str(test_features.iloc[row_index]["cas_subtype"])
        rows.append(
            {
                "true_subtype": true_label,
                "predicted_subtype": predicted,
                "confidence": confidence,
                "correct": predicted == true_label,
                "genome_id": test.iloc[row_index].get("genome_id", ""),
                "organism": test.iloc[row_index].get("organism", ""),
                "repeat_sequence": test.iloc[row_index].get("repeat_sequence", ""),
            }
        )
    predictions = pd.DataFrame(rows)
    bins = _confidence_bins(predictions, n_bins=n_bins)
    subtype = _subtype_confidence(predictions)
    summary = pd.DataFrame(
        [
            {
                "training_csv": str(training_csv),
                "split_strategy": split_strategy,
                "group_column": group_column,
                "train_rows": len(train),
                "test_rows": len(test),
                "accuracy": float(predictions["correct"].mean()),
                "mean_confidence": float(predictions["confidence"].mean()),
                "expected_calibration_error": _expected_calibration_error(bins),
                "min_class_count": min_class_count,
                "n_bins": n_bins,
            }
        ]
    )
    thresholds = _accuracy_by_threshold(predictions)
    return {
        "summary": summary,
        "predictions": predictions,
        "confidence_bins": bins,
        "subtype_confidence": subtype,
        "accuracy_by_threshold": thresholds,
    }


def write_calibration_outputs(outputs: dict[str, pd.DataFrame], output_dir: str | Path) -> None:
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    for name, table in outputs.items():
        table.to_csv(output_path / f"{name}.csv", index=False)
    _plot_reliability(outputs["confidence_bins"], output_path / "reliability_curve.png")
    _plot_accuracy_by_threshold(
        outputs["accuracy_by_threshold"],
        output_path / "accuracy_by_confidence_threshold.png",
    )
    _plot_subtype_confidence(
        outputs["subtype_confidence"],
        output_path / "subtype_confidence_accuracy.png",
    )


def _confidence_bins(predictions: pd.DataFrame, n_bins: int) -> pd.DataFrame:
    bins = []
    for index in range(n_bins):
        lower = index / n_bins
        upper = (index + 1) / n_bins
        if index == n_bins - 1:
            subset = predictions[(predictions["confidence"] >= lower) & (predictions["confidence"] <= upper)]
        else:
            subset = predictions[(predictions["confidence"] >= lower) & (predictions["confidence"] < upper)]
        if subset.empty:
            accuracy = 0.0
            mean_confidence = (lower + upper) / 2
        else:
            accuracy = float(subset["correct"].mean())
            mean_confidence = float(subset["confidence"].mean())
        bins.append(
            {
                "bin_lower": lower,
                "bin_upper": upper,
                "row_count": len(subset),
                "mean_confidence": mean_confidence,
                "accuracy": accuracy,
                "calibration_gap": accuracy - mean_confidence,
            }
        )
    return pd.DataFrame(bins)


def _expected_calibration_error(bins: pd.DataFrame) -> float:
    total = int(bins["row_count"].sum())
    if total == 0:
        return 0.0
    weighted_gap = (
        bins["row_count"] * (bins["accuracy"] - bins["mean_confidence"]).abs()
    ).sum()
    return float(weighted_gap / total)


def _subtype_confidence(predictions: pd.DataFrame) -> pd.DataFrame:
    return (
        predictions.groupby("true_subtype")
        .agg(
            row_count=("correct", "size"),
            accuracy=("correct", "mean"),
            mean_confidence=("confidence", "mean"),
            median_confidence=("confidence", "median"),
        )
        .reset_index()
        .sort_values(["accuracy", "row_count", "true_subtype"], ascending=[True, False, True])
    )


def _accuracy_by_threshold(predictions: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for threshold in [round(value / 100, 2) for value in range(0, 101, 5)]:
        subset = predictions[predictions["confidence"] >= threshold]
        rows.append(
            {
                "confidence_threshold": threshold,
                "row_count": len(subset),
                "coverage_fraction": len(subset) / len(predictions) if len(predictions) else 0.0,
                "accuracy": float(subset["correct"].mean()) if not subset.empty else 0.0,
            }
        )
    return pd.DataFrame(rows)


def _plot_reliability(bins: pd.DataFrame, output_path: Path) -> None:
    plt.figure(figsize=(6, 6))
    plt.plot([0, 1], [0, 1], color="#555555", linestyle="--", label="perfect calibration")
    nonempty = bins[bins["row_count"] > 0]
    plt.plot(nonempty["mean_confidence"], nonempty["accuracy"], marker="o", color="#2F6F73")
    plt.xlabel("Mean predicted confidence")
    plt.ylabel("Observed accuracy")
    plt.title("Subtype Confidence Reliability")
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(output_path, dpi=180)
    plt.close()


def _plot_accuracy_by_threshold(thresholds: pd.DataFrame, output_path: Path) -> None:
    plt.figure(figsize=(7, 5))
    plt.plot(thresholds["confidence_threshold"], thresholds["accuracy"], marker="o", color="#2F6F73")
    plt.plot(thresholds["confidence_threshold"], thresholds["coverage_fraction"], marker="s", color="#9A6B2F")
    plt.xlabel("Minimum confidence threshold")
    plt.ylabel("Fraction")
    plt.title("Accuracy and Coverage by Confidence Threshold")
    plt.ylim(0, 1.05)
    plt.legend(["accuracy", "coverage"], frameon=False)
    plt.tight_layout()
    plt.savefig(output_path, dpi=180)
    plt.close()


def _plot_subtype_confidence(subtype: pd.DataFrame, output_path: Path) -> None:
    ordered = subtype.sort_values("accuracy")
    labels = ordered["true_subtype"]
    x = range(len(ordered))
    plt.figure(figsize=(10, 5))
    plt.bar(x, ordered["accuracy"], color="#2F6F73", alpha=0.82, label="accuracy")
    plt.scatter(x, ordered["mean_confidence"], color="#9A6B2F", label="mean confidence", zorder=3)
    plt.xticks(x, labels, rotation=45, ha="right")
    plt.ylabel("Fraction")
    plt.ylim(0, 1.05)
    plt.title("Observed Accuracy and Mean Confidence by True Subtype")
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(output_path, dpi=180)
    plt.close()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Calibrate ExtraTrees Cas subtype probabilities on a held-out split."
    )
    parser.add_argument("training_csv", help="Repeat/Cas subtype training CSV.")
    parser.add_argument("--output-dir", required=True, help="Directory for calibration outputs.")
    parser.add_argument(
        "--split-strategy",
        choices=["row_random", "group_holdout"],
        default="group_holdout",
    )
    parser.add_argument("--group-column", default="genus")
    parser.add_argument("--test-size", type=float, default=0.25)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("--min-class-count", type=int, default=20)
    parser.add_argument("--n-bins", type=int, default=10)
    args = parser.parse_args()

    outputs = calibrate_extra_trees_subtype_model(
        training_csv=args.training_csv,
        split_strategy=args.split_strategy,
        group_column=args.group_column,
        test_size=args.test_size,
        random_state=args.random_state,
        min_class_count=args.min_class_count,
        n_bins=args.n_bins,
    )
    write_calibration_outputs(outputs, args.output_dir)
    print(outputs["summary"].to_string(index=False))
    print(f"Wrote calibration outputs to {args.output_dir}")


if __name__ == "__main__":
    main()
