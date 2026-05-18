from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
from sklearn.model_selection import train_test_split

from crispr_phage_predictor.ml.classifier import NearestRepeatClassifier, RepeatCasSubtypeClassifier
from crispr_phage_predictor.ml.dataset import (
    filter_high_confidence_labels,
    load_repeat_cas_training_table,
)


@dataclass(frozen=True)
class EvaluationResult:
    method: str
    accuracy: float
    labels: list[str]
    confusion: pd.DataFrame
    report: str
    train_size: int
    test_size: int


def evaluate_classifier(
    training_table: pd.DataFrame,
    test_size: float = 0.25,
    random_state: int = 42,
    high_confidence_only: bool = True,
) -> EvaluationResult:
    table = filter_high_confidence_labels(training_table) if high_confidence_only else training_table.copy()
    _validate_evaluation_table(table)

    stratify = table["cas_subtype"] if _can_stratify(table["cas_subtype"]) else None
    train, test = train_test_split(
        table,
        test_size=test_size,
        random_state=random_state,
        stratify=stratify,
    )

    classifier = RepeatCasSubtypeClassifier(random_state=random_state)
    classifier.fit(train)
    predicted_labels = [prediction.cas_subtype for prediction in classifier.predict_table(test)]
    return _build_evaluation_result(
        method="random_forest",
        true_labels=list(test["cas_subtype"]),
        predicted_labels=predicted_labels,
        train_size=len(train),
        test_size=len(test),
    )


def evaluate_methods(
    training_table: pd.DataFrame,
    test_size: float = 0.25,
    random_state: int = 42,
    high_confidence_only: bool = True,
) -> list[EvaluationResult]:
    table = filter_high_confidence_labels(training_table) if high_confidence_only else training_table.copy()
    _validate_evaluation_table(table)

    stratify = table["cas_subtype"] if _can_stratify(table["cas_subtype"]) else None
    train, test = train_test_split(
        table,
        test_size=test_size,
        random_state=random_state,
        stratify=stratify,
    )

    true_labels = list(test["cas_subtype"])

    nearest = NearestRepeatClassifier().fit(train)
    nearest_predictions = [prediction.cas_subtype for prediction in nearest.predict_table(test)]

    random_forest = RepeatCasSubtypeClassifier(random_state=random_state).fit(train)
    random_forest_predictions = [
        prediction.cas_subtype for prediction in random_forest.predict_table(test)
    ]

    return [
        _build_evaluation_result(
            method="nearest_repeat",
            true_labels=true_labels,
            predicted_labels=nearest_predictions,
            train_size=len(train),
            test_size=len(test),
        ),
        _build_evaluation_result(
            method="random_forest",
            true_labels=true_labels,
            predicted_labels=random_forest_predictions,
            train_size=len(train),
            test_size=len(test),
        ),
    ]


def _build_evaluation_result(
    method: str,
    true_labels: list[str],
    predicted_labels: list[str],
    train_size: int,
    test_size: int,
) -> EvaluationResult:
    labels = sorted(set(true_labels) | set(predicted_labels))

    return EvaluationResult(
        method=method,
        accuracy=accuracy_score(true_labels, predicted_labels),
        labels=labels,
        confusion=pd.DataFrame(
            confusion_matrix(true_labels, predicted_labels, labels=labels),
            index=[f"true_{label}" for label in labels],
            columns=[f"pred_{label}" for label in labels],
        ),
        report=classification_report(true_labels, predicted_labels, labels=labels, zero_division=0),
        train_size=train_size,
        test_size=test_size,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Train and evaluate the baseline CRISPR repeat Cas-subtype classifier."
    )
    parser.add_argument("training_csv", type=Path, help="Path to repeats_cas_types.csv")
    parser.add_argument("--test-size", type=float, default=0.25)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument(
        "--include-medium-confidence",
        action="store_true",
        help="Use all validated labels instead of filtering to high-confidence labels.",
    )
    args = parser.parse_args()

    table = load_repeat_cas_training_table(args.training_csv)
    results = evaluate_methods(
        table,
        test_size=args.test_size,
        random_state=args.random_state,
        high_confidence_only=not args.include_medium_confidence,
    )

    for result in results:
        print(f"Method: {result.method}")
        print(f"Rows used: {result.train_size + result.test_size}")
        print(f"Train rows: {result.train_size}")
        print(f"Test rows: {result.test_size}")
        print(f"Accuracy: {result.accuracy:.4f}")
        print()
        print("Classification report")
        print(result.report)
        print("Confusion matrix")
        print(result.confusion.to_string())
        print()


def _validate_evaluation_table(table: pd.DataFrame) -> None:
    if table.empty:
        raise ValueError("No rows available for evaluation after filtering")
    subtype_counts = table["cas_subtype"].value_counts()
    if len(subtype_counts) < 2:
        raise ValueError("At least two Cas subtypes are required for evaluation")
    if len(table) < 4:
        raise ValueError("At least four rows are required for a train/test evaluation")


def _can_stratify(labels: pd.Series) -> bool:
    return bool((labels.value_counts() >= 2).all())


if __name__ == "__main__":
    main()
