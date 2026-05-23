from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix

from crispr_phage_predictor.ml.dataset import load_repeat_cas_training_table
from crispr_phage_predictor.ml.features import build_repeat_feature_table, feature_columns
from crispr_phage_predictor.ml.train_classifier import _filter_min_class_count


def evaluate_external_dataset(
    train_csv: str | Path,
    test_csv: str | Path,
    min_class_count: int = 20,
    random_state: int = 42,
    n_estimators: int = 400,
) -> dict[str, object]:
    """Train on one repeat/Cas table and evaluate on an independent table."""
    train = _filter_min_class_count(
        load_repeat_cas_training_table(train_csv),
        min_class_count=min_class_count,
    )
    raw_test = load_repeat_cas_training_table(test_csv)
    train_labels = set(train["cas_subtype"].astype(str))
    test = raw_test[raw_test["cas_subtype"].astype(str).isin(train_labels)].copy()
    excluded = raw_test[~raw_test["cas_subtype"].astype(str).isin(train_labels)].copy()

    train_features = build_repeat_feature_table(train)
    test_features = build_repeat_feature_table(test)
    features = feature_columns(train_features)
    aligned_test = test_features.reindex(columns=features, fill_value=0.0)

    model = ExtraTreesClassifier(
        n_estimators=n_estimators,
        random_state=random_state,
        class_weight="balanced",
        n_jobs=-1,
    )
    model.fit(train_features[features], train_features["cas_subtype"])
    predictions = list(model.predict(aligned_test))
    true_labels = list(test_features["cas_subtype"])
    labels = sorted(set(true_labels) | set(predictions))

    return {
        "train_rows": len(train),
        "raw_test_rows": len(raw_test),
        "evaluated_test_rows": len(test),
        "excluded_test_rows": len(excluded),
        "excluded_test_subtypes": sorted(excluded["cas_subtype"].astype(str).unique()),
        "accuracy": accuracy_score(true_labels, predictions),
        "report": classification_report(true_labels, predictions, labels=labels, zero_division=0),
        "confusion": pd.DataFrame(
            confusion_matrix(true_labels, predictions, labels=labels),
            index=[f"true_{label}" for label in labels],
            columns=[f"pred_{label}" for label in labels],
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Train on one repeat/Cas table and evaluate on another table."
    )
    parser.add_argument("train_csv", help="Training repeat/Cas CSV.")
    parser.add_argument("test_csv", help="External test repeat/Cas CSV.")
    parser.add_argument("--min-class-count", type=int, default=20)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("--n-estimators", type=int, default=400)
    args = parser.parse_args()

    result = evaluate_external_dataset(
        train_csv=args.train_csv,
        test_csv=args.test_csv,
        min_class_count=args.min_class_count,
        random_state=args.random_state,
        n_estimators=args.n_estimators,
    )
    print(f"Train rows: {result['train_rows']}")
    print(f"Raw test rows: {result['raw_test_rows']}")
    print(f"Evaluated test rows: {result['evaluated_test_rows']}")
    print(f"Excluded test rows: {result['excluded_test_rows']}")
    print(f"Excluded test subtypes: {', '.join(result['excluded_test_subtypes'])}")
    print(f"Accuracy: {result['accuracy']:.4f}")
    print()
    print("Classification report")
    print(result["report"])
    print("Confusion matrix")
    print(result["confusion"].to_string())


if __name__ == "__main__":
    main()
