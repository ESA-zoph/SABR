from __future__ import annotations

import argparse
from dataclasses import dataclass
from math import log2
from pathlib import Path

import pandas as pd
from sklearn.ensemble import ExtraTreesClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
from sklearn.model_selection import train_test_split
from sklearn.multiclass import OneVsRestClassifier
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC

from crispr_phage_predictor.ml.classifier import NearestRepeatClassifier, RepeatCasSubtypeClassifier
from crispr_phage_predictor.ml.dataset import (
    filter_high_confidence_labels,
    load_repeat_cas_training_table,
)
from crispr_phage_predictor.ml.features import build_repeat_feature_table, feature_columns


@dataclass(frozen=True)
class EvaluationResult:
    method: str
    accuracy: float
    labels: list[str]
    confusion: pd.DataFrame
    report: str
    train_size: int
    test_size: int
    split_strategy: str = "row_random"


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
    split_strategy: str = "row_random",
    group_column: str = "genome_id",
    min_class_count: int = 1,
    methods: tuple[str, ...] = (
        "nearest_repeat",
        "logistic_regression",
        "linear_svm",
        "gradient_boosting",
        "extra_trees",
        "hybrid_extra_trees",
        "random_forest",
    ),
) -> list[EvaluationResult]:
    table = filter_high_confidence_labels(training_table) if high_confidence_only else training_table.copy()
    table = _filter_min_class_count(table, min_class_count=min_class_count)
    _validate_evaluation_table(table)

    train, test = _split_table(
        table=table,
        test_size=test_size,
        random_state=random_state,
        split_strategy=split_strategy,
        group_column=group_column,
    )

    true_labels = list(test["cas_subtype"])

    results = []
    if "nearest_repeat" in methods:
        nearest = NearestRepeatClassifier().fit(train)
        nearest_predictions = [prediction.cas_subtype for prediction in nearest.predict_table(test)]
        results.append(
            _build_evaluation_result(
                method="nearest_repeat",
                true_labels=true_labels,
                predicted_labels=nearest_predictions,
                train_size=len(train),
                test_size=len(test),
                split_strategy=split_strategy,
            )
        )

    results.extend(_evaluate_feature_models(
        train=train,
        test=test,
        random_state=random_state,
        split_strategy=split_strategy,
        methods=methods,
    ))

    if "random_forest" in methods:
        random_forest = RepeatCasSubtypeClassifier(random_state=random_state).fit(train)
        random_forest_predictions = [
            prediction.cas_subtype for prediction in random_forest.predict_table(test)
        ]
        results.append(
            _build_evaluation_result(
                method="random_forest",
                true_labels=true_labels,
                predicted_labels=random_forest_predictions,
                train_size=len(train),
                test_size=len(test),
                split_strategy=split_strategy,
            )
        )
    return results


def _evaluate_feature_models(
    train: pd.DataFrame,
    test: pd.DataFrame,
    random_state: int,
    split_strategy: str,
    methods: tuple[str, ...],
) -> list[EvaluationResult]:
    train_features = build_repeat_feature_table(train)
    test_features = build_repeat_feature_table(test)
    features = feature_columns(train_features)
    x_train = train_features[features]
    y_train = train_features["cas_subtype"]
    x_test = test_features.reindex(columns=features, fill_value=0.0)
    y_test = list(test_features["cas_subtype"])

    models = [
        (
            "logistic_regression",
            make_pipeline(
                StandardScaler(),
                OneVsRestClassifier(
                    LogisticRegression(
                        max_iter=1000,
                        class_weight="balanced",
                        random_state=random_state,
                        solver="liblinear",
                    )
                ),
            ),
        ),
        (
            "linear_svm",
            make_pipeline(
                StandardScaler(),
                LinearSVC(
                    class_weight="balanced",
                    random_state=random_state,
                    dual="auto",
                    max_iter=5000,
                ),
            ),
        ),
        (
            "gradient_boosting",
            GradientBoostingClassifier(random_state=random_state),
        ),
        (
            "extra_trees",
            ExtraTreesClassifier(
                n_estimators=400,
                random_state=random_state,
                class_weight="balanced",
                n_jobs=-1,
            ),
        ),
    ]

    results = []
    for method, model in models:
        if method not in methods:
            continue
        model.fit(x_train, y_train)
        predictions = list(model.predict(x_test))
        results.append(
            _build_evaluation_result(
                method=method,
                true_labels=y_test,
                predicted_labels=predictions,
                train_size=len(train),
                test_size=len(test),
                split_strategy=split_strategy,
            )
        )
    if "hybrid_extra_trees" in methods:
        hybrid_x_train = _feature_table_with_neighbor_features(train, train, exclude_self=True)
        hybrid_x_test = _feature_table_with_neighbor_features(test, train, exclude_self=False)
        model = ExtraTreesClassifier(
            n_estimators=400,
            random_state=random_state,
            class_weight="balanced",
            n_jobs=-1,
        )
        model.fit(hybrid_x_train, y_train)
        predictions = list(model.predict(hybrid_x_test))
        results.append(
            _build_evaluation_result(
                method="hybrid_extra_trees",
                true_labels=y_test,
                predicted_labels=predictions,
                train_size=len(train),
                test_size=len(test),
                split_strategy=split_strategy,
            )
        )
    if "hierarchical_extra_trees" in methods:
        predictions = _predict_with_hierarchical_extra_trees(
            train_features=train_features,
            test_features=test_features,
            feature_names=features,
            random_state=random_state,
        )
        results.append(
            _build_evaluation_result(
                method="hierarchical_extra_trees",
                true_labels=y_test,
                predicted_labels=predictions,
                train_size=len(train),
                test_size=len(test),
                split_strategy=split_strategy,
            )
        )
    return results


def _predict_with_hierarchical_extra_trees(
    train_features: pd.DataFrame,
    test_features: pd.DataFrame,
    feature_names: list[str],
    random_state: int,
) -> list[str]:
    x_train = train_features[feature_names]
    x_test = test_features.reindex(columns=feature_names, fill_value=0.0)
    type_model = ExtraTreesClassifier(
        n_estimators=400,
        random_state=random_state,
        class_weight="balanced",
        n_jobs=-1,
    )
    type_model.fit(x_train, train_features["cas_type"])

    subtype_models = {}
    constant_subtypes = {}
    global_model = ExtraTreesClassifier(
        n_estimators=400,
        random_state=random_state,
        class_weight="balanced",
        n_jobs=-1,
    )
    global_model.fit(x_train, train_features["cas_subtype"])

    for cas_type, type_rows in train_features.groupby("cas_type"):
        subtype_count = type_rows["cas_subtype"].nunique()
        if subtype_count == 1:
            constant_subtypes[cas_type] = str(type_rows["cas_subtype"].iloc[0])
            continue
        model = ExtraTreesClassifier(
            n_estimators=400,
            random_state=random_state,
            class_weight="balanced",
            n_jobs=-1,
        )
        model.fit(type_rows[feature_names], type_rows["cas_subtype"])
        subtype_models[cas_type] = model

    predicted_types = list(type_model.predict(x_test))
    predictions = []
    for row_index, predicted_type in enumerate(predicted_types):
        row = x_test.iloc[[row_index]]
        if predicted_type in constant_subtypes:
            predictions.append(constant_subtypes[predicted_type])
        elif predicted_type in subtype_models:
            predictions.append(str(subtype_models[predicted_type].predict(row)[0]))
        else:
            predictions.append(str(global_model.predict(row)[0]))
    return predictions


def _build_evaluation_result(
    method: str,
    true_labels: list[str],
    predicted_labels: list[str],
    train_size: int,
    test_size: int,
    split_strategy: str = "row_random",
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
        split_strategy=split_strategy,
    )


def _feature_table_with_neighbor_features(
    query_table: pd.DataFrame,
    reference_table: pd.DataFrame,
    exclude_self: bool,
) -> pd.DataFrame:
    query_features = build_repeat_feature_table(query_table)
    feature_names = feature_columns(query_features)
    numeric = query_features[feature_names].reset_index(drop=True)
    reference_examples = _reference_neighbor_examples(reference_table)
    neighbor_rows = [
        _nearest_neighbor_feature_row(
            query_row,
            reference_examples,
            exclude_self=exclude_self,
        )
        for _, query_row in query_table.reset_index(drop=True).iterrows()
    ]
    neighbor_features = pd.DataFrame(neighbor_rows).reset_index(drop=True)
    return pd.concat([numeric, neighbor_features], axis=1)


def _nearest_neighbor_feature_row(
    query_row: pd.Series,
    reference_examples: list[tuple[str, str, str, str, str]],
    exclude_self: bool,
    close_threshold: float = 0.9,
) -> dict[str, float]:
    query_repeat = str(query_row["repeat_sequence"]).upper()
    query_genome = str(query_row.get("genome_id", ""))
    query_contig = str(query_row.get("contig_id", ""))
    query_start = str(query_row.get("array_start", ""))
    scored: list[tuple[float, str]] = []
    for reference_repeat, subtype, reference_genome, reference_contig, reference_start in reference_examples:
        if exclude_self and _same_array_identity(
            query_genome,
            query_contig,
            query_start,
            reference_genome,
            reference_contig,
            reference_start,
        ):
            continue
        scored.append((_global_identity(query_repeat, reference_repeat), subtype))

    if not scored:
        return {
            "neighbor_best_identity": 0.0,
            "neighbor_second_identity": 0.0,
            "neighbor_identity_margin": 0.0,
            "neighbor_close_count_90": 0.0,
            "neighbor_top_vote_fraction_90": 0.0,
            "neighbor_subtype_entropy_90": 0.0,
        }

    scored.sort(key=lambda item: item[0], reverse=True)
    best_identity = scored[0][0]
    second_identity = scored[1][0] if len(scored) > 1 else 0.0
    close_subtypes = [subtype for identity, subtype in scored if identity >= close_threshold]
    close_counts = pd.Series(close_subtypes).value_counts() if close_subtypes else pd.Series(dtype=int)
    top_vote_fraction = (
        float(close_counts.iloc[0] / close_counts.sum()) if not close_counts.empty else 0.0
    )
    entropy = _entropy([int(value) for value in close_counts.tolist()])
    return {
        "neighbor_best_identity": round(best_identity, 6),
        "neighbor_second_identity": round(second_identity, 6),
        "neighbor_identity_margin": round(best_identity - second_identity, 6),
        "neighbor_close_count_90": float(len(close_subtypes)),
        "neighbor_top_vote_fraction_90": round(top_vote_fraction, 6),
        "neighbor_subtype_entropy_90": round(entropy, 6),
    }


def _reference_neighbor_examples(reference_table: pd.DataFrame) -> list[tuple[str, str, str, str, str]]:
    return [
        (
            str(row["repeat_sequence"]).upper(),
            str(row["cas_subtype"]),
            str(row.get("genome_id", "")),
            str(row.get("contig_id", "")),
            str(row.get("array_start", "")),
        )
        for _, row in reference_table.iterrows()
    ]


def _same_array_identity(
    query_genome: str,
    query_contig: str,
    query_start: str,
    reference_genome: str,
    reference_contig: str,
    reference_start: str,
) -> bool:
    return (
        query_genome == reference_genome
        and query_contig == reference_contig
        and query_start == reference_start
    )


def _entropy(counts: list[int]) -> float:
    total = sum(counts)
    if total <= 0:
        return 0.0
    return -sum((count / total) * log2(count / total) for count in counts if count)


def _global_identity(left: str, right: str) -> float:
    if not left and not right:
        return 1.0
    max_length = max(len(left), len(right))
    if max_length == 0:
        return 0.0
    matches = sum(1 for left_base, right_base in zip(left, right) if left_base == right_base)
    return matches / max_length


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Train and evaluate the baseline CRISPR repeat Cas-subtype classifier."
    )
    parser.add_argument("training_csv", type=Path, help="Path to repeats_cas_types.csv")
    parser.add_argument("--test-size", type=float, default=0.25)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument(
        "--split-strategy",
        choices=["row_random", "group_holdout"],
        default="row_random",
        help="Use row_random for smoke tests or group_holdout for accession/genome holdout.",
    )
    parser.add_argument(
        "--group-column",
        default="genome_id",
        help="Column to hold out as groups when --split-strategy group_holdout is used.",
    )
    parser.add_argument(
        "--methods",
        default="nearest_repeat,logistic_regression,linear_svm,random_forest",
        help=(
            "Comma-separated methods to evaluate. Available: nearest_repeat, "
            "logistic_regression, linear_svm, gradient_boosting, extra_trees, "
            "hybrid_extra_trees, hierarchical_extra_trees, random_forest."
        ),
    )
    parser.add_argument(
        "--min-class-count",
        type=int,
        default=1,
        help="Drop Cas subtypes with fewer than this many rows before evaluation.",
    )
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
        split_strategy=args.split_strategy,
        group_column=args.group_column,
        min_class_count=args.min_class_count,
        methods=tuple(method.strip() for method in args.methods.split(",") if method.strip()),
    )

    for result in results:
        print(f"Method: {result.method}")
        print(f"Rows used: {result.train_size + result.test_size}")
        print(f"Split strategy: {result.split_strategy}")
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


def _filter_min_class_count(table: pd.DataFrame, min_class_count: int) -> pd.DataFrame:
    if min_class_count <= 1:
        return table.copy()
    counts = table["cas_subtype"].value_counts()
    kept = set(counts[counts >= min_class_count].index)
    return table[table["cas_subtype"].isin(kept)].copy()


def _split_table(
    table: pd.DataFrame,
    test_size: float,
    random_state: int,
    split_strategy: str,
    group_column: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if split_strategy == "row_random":
        stratify = table["cas_subtype"] if _can_stratify(table["cas_subtype"]) else None
        return train_test_split(
            table,
            test_size=test_size,
            random_state=random_state,
            stratify=stratify,
        )

    if split_strategy != "group_holdout":
        raise ValueError("split_strategy must be row_random or group_holdout")

    groups = _group_values(table, group_column)
    unique_groups = list(groups.drop_duplicates())
    if len(unique_groups) < 2:
        raise ValueError("At least two groups are required for group_holdout evaluation")

    shuffled = pd.Series(unique_groups).sample(frac=1.0, random_state=random_state).tolist()
    target_test_rows = max(1, int(round(len(table) * test_size)))
    test_groups: set[str] = set()
    test_rows = 0
    group_sizes = groups.value_counts().to_dict()
    for group in shuffled:
        test_groups.add(group)
        test_rows += int(group_sizes[group])
        if test_rows >= target_test_rows:
            break

    test_mask = groups.isin(test_groups)
    train = table[~test_mask].copy()
    test = table[test_mask].copy()
    _validate_split(train, test)
    return train, test


def _group_values(table: pd.DataFrame, group_column: str) -> pd.Series:
    if group_column in table.columns:
        values = table[group_column].astype(str).fillna("")
    elif group_column == "genus":
        values = table.apply(_derive_genus, axis=1)
    elif group_column == "species":
        values = table.apply(_derive_species, axis=1)
    else:
        raise ValueError(f"Group column not found: {group_column}")

    values = values.replace({"": "unknown"})
    if (values == "unknown").all():
        raise ValueError(f"Group column {group_column} could not be derived")
    return values


def _derive_genus(row: pd.Series) -> str:
    organism = str(row.get("organism", "")).strip()
    if organism:
        return organism.split()[0]
    taxonomy = str(row.get("taxonomy", "")).strip()
    if taxonomy:
        return taxonomy.split(";")[-1].strip() or "unknown"
    return "unknown"


def _derive_species(row: pd.Series) -> str:
    organism = str(row.get("organism", "")).strip()
    parts = organism.split()
    if len(parts) >= 2:
        return " ".join(parts[:2])
    return _derive_genus(row)


def _validate_split(train: pd.DataFrame, test: pd.DataFrame) -> None:
    if train.empty or test.empty:
        raise ValueError("Train and test splits must both contain rows")
    if train["cas_subtype"].nunique() < 2:
        raise ValueError("Training split must contain at least two Cas subtypes")
    unseen_labels = sorted(set(test["cas_subtype"]) - set(train["cas_subtype"]))
    if unseen_labels:
        raise ValueError(
            "Test split contains labels absent from training split: " + ", ".join(unseen_labels)
        )


def _can_stratify(labels: pd.Series) -> bool:
    return bool((labels.value_counts() >= 2).all())


if __name__ == "__main__":
    main()
