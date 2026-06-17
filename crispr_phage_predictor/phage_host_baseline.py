from __future__ import annotations

from dataclasses import dataclass

import pandas as pd
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, f1_score
from sklearn.model_selection import GroupShuffleSplit, train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


NON_FEATURE_COLUMNS = {
    "interaction_id",
    "source_key",
    "bacterium",
    "strain",
    "phage",
    "eop_class",
    "susceptibility_label",
    "binary_susceptibility",
    "eop_value",
    "host_accession",
    "host_local_path",
    "phage_accession",
    "phage_local_path",
    "dataset_tier",
}


@dataclass(frozen=True)
class PhageHostEvaluation:
    method: str
    split_strategy: str
    train_size: int
    test_size: int
    accuracy: float
    macro_f1: float
    labels: list[str]
    confusion: pd.DataFrame
    report: str


def evaluate_baseline_models(
    feature_table: pd.DataFrame,
    split_strategy: str = "group_by_phage",
    test_size: float = 0.3,
    random_state: int = 42,
) -> list[PhageHostEvaluation]:
    table = _modeling_table(feature_table)
    train, test = _split_table(table, split_strategy, test_size, random_state)
    feature_names = numeric_feature_columns(train)
    x_train = train[feature_names]
    y_train = train["binary_susceptibility"]
    x_test = test.reindex(columns=train.columns)[feature_names]
    y_test = test["binary_susceptibility"]

    models = [
        ("majority_baseline", DummyClassifier(strategy="most_frequent")),
        (
            "logistic_regression",
            make_pipeline(
                StandardScaler(),
                LogisticRegression(
                    max_iter=2000,
                    class_weight="balanced",
                    random_state=random_state,
                ),
            ),
        ),
        (
            "random_forest",
            RandomForestClassifier(
                n_estimators=400,
                random_state=random_state,
                class_weight="balanced",
                n_jobs=-1,
            ),
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
        model.fit(x_train, y_train)
        predictions = model.predict(x_test)
        results.append(
            _evaluation_result(
                method=method,
                split_strategy=split_strategy,
                y_true=list(y_test),
                y_pred=list(predictions),
                train_size=len(train),
                test_size=len(test),
            )
        )
    return results


def numeric_feature_columns(table: pd.DataFrame) -> list[str]:
    excluded = NON_FEATURE_COLUMNS | {
        "host_linkage_status",
        "phage_linkage_status",
    }
    features = []
    for column in table.columns:
        if column in excluded:
            continue
        if pd.api.types.is_bool_dtype(table[column]) or pd.api.types.is_numeric_dtype(table[column]):
            features.append(column)
    return features


def _modeling_table(feature_table: pd.DataFrame) -> pd.DataFrame:
    table = feature_table.copy()
    table = table[table["binary_susceptibility"].isin(["susceptible", "resistant"])].copy()
    table = table[table["eop_class"] != "not_reported"].copy()
    if table["binary_susceptibility"].nunique() < 2:
        raise ValueError("Need at least two binary susceptibility classes")
    for column in table.columns:
        if pd.api.types.is_bool_dtype(table[column]):
            table[column] = table[column].astype(int)
    return table


def _split_table(
    table: pd.DataFrame,
    split_strategy: str,
    test_size: float,
    random_state: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if split_strategy == "row_random":
        stratify = table["binary_susceptibility"] if _can_stratify(table["binary_susceptibility"]) else None
        return train_test_split(
            table,
            test_size=test_size,
            random_state=random_state,
            stratify=stratify,
        )
    group_column = {
        "group_by_phage": "phage",
        "group_by_source": "source_key",
    }.get(split_strategy)
    if group_column is None:
        raise ValueError(f"Unknown split strategy: {split_strategy}")
    splitter = GroupShuffleSplit(
        n_splits=1,
        test_size=test_size,
        random_state=random_state,
    )
    train_index, test_index = next(
        splitter.split(table, table["binary_susceptibility"], groups=table[group_column])
    )
    return table.iloc[train_index].copy(), table.iloc[test_index].copy()


def _evaluation_result(
    method: str,
    split_strategy: str,
    y_true: list[str],
    y_pred: list[str],
    train_size: int,
    test_size: int,
) -> PhageHostEvaluation:
    labels = sorted(set(y_true) | set(y_pred))
    return PhageHostEvaluation(
        method=method,
        split_strategy=split_strategy,
        train_size=train_size,
        test_size=test_size,
        accuracy=round(float(accuracy_score(y_true, y_pred)), 6),
        macro_f1=round(float(f1_score(y_true, y_pred, labels=labels, average="macro")), 6),
        labels=labels,
        confusion=pd.DataFrame(
            confusion_matrix(y_true, y_pred, labels=labels),
            index=[f"true_{label}" for label in labels],
            columns=[f"pred_{label}" for label in labels],
        ),
        report=classification_report(y_true, y_pred, labels=labels, zero_division=0),
    )


def _can_stratify(labels: pd.Series) -> bool:
    return bool((labels.value_counts() >= 2).all())
