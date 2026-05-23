from __future__ import annotations

from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.metrics import accuracy_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from crispr_phage_predictor.ml.dataset import load_repeat_cas_training_table
from crispr_phage_predictor.ml.features import build_repeat_feature_table, feature_columns
from crispr_phage_predictor.ml.train_classifier import _filter_min_class_count, _split_table


DATASET = ROOT / "data" / "training" / "repeats_cas_types_crisprcasdb_sql_candidate.csv"
OUTPUT_DIR = ROOT / "docs" / "interpretability" / "crisprcasdb_extratrees"
ASSET_DIR = ROOT / "docs" / "manuscript_assets"


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    ASSET_DIR.mkdir(parents=True, exist_ok=True)
    train_features, test_features, feature_names, model = _fit_model()

    builtin = _builtin_importance(model, feature_names)
    builtin.to_csv(OUTPUT_DIR / "builtin_feature_importance.csv", index=False)
    _plot_top_features(
        builtin,
        "ExtraTrees Built-In Feature Importance",
        ASSET_DIR / "crisprcasdb_builtin_feature_importance",
    )

    permutation_features = _permutation_feature_subset(builtin, feature_names, top_n=40)
    permutation = _permutation_importance(model, test_features, feature_names, permutation_features)
    permutation.to_csv(OUTPUT_DIR / "permutation_feature_importance.csv", index=False)
    _plot_top_features(
        permutation.rename(columns={"importance_mean": "importance"}),
        "Held-Out Permutation Feature Importance",
        ASSET_DIR / "crisprcasdb_permutation_feature_importance",
    )

    category = _category_summary(builtin, permutation)
    category.to_csv(OUTPUT_DIR / "feature_category_importance.csv", index=False)
    _plot_category_summary(category, ASSET_DIR / "crisprcasdb_feature_category_importance")

    error_focus = _typeiii_error_feature_summary(test_features)
    error_focus.to_csv(OUTPUT_DIR / "typeiii_error_feature_summary.csv", index=False)
    _plot_typeiii_error_summary(error_focus, ASSET_DIR / "crisprcasdb_typeiii_error_feature_summary")
    print(f"Wrote interpretability outputs to {OUTPUT_DIR}")


def _fit_model():
    table = load_repeat_cas_training_table(DATASET)
    table = _filter_min_class_count(table, min_class_count=20)
    train, test = _split_table(
        table=table,
        test_size=0.25,
        random_state=42,
        split_strategy="group_holdout",
        group_column="genome_id",
    )
    train_features = build_repeat_feature_table(train)
    test_features = build_repeat_feature_table(test)
    feature_names = feature_columns(train_features)
    aligned_test = test_features.reindex(columns=feature_names, fill_value=0.0)
    model = ExtraTreesClassifier(
        n_estimators=400,
        random_state=42,
        class_weight="balanced",
        n_jobs=-1,
    )
    model.fit(train_features[feature_names], train_features["cas_subtype"])
    aligned_test["cas_subtype"] = list(test_features["cas_subtype"])
    predictions = model.predict(aligned_test[feature_names])
    aligned_test["predicted_subtype"] = predictions
    aligned_test["correct"] = aligned_test["predicted_subtype"] == aligned_test["cas_subtype"]
    return train_features, aligned_test, feature_names, model


def _builtin_importance(model: ExtraTreesClassifier, feature_names: list[str]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "feature": feature_names,
            "importance": model.feature_importances_,
        }
    ).assign(category=lambda table: table["feature"].map(_feature_category)).sort_values(
        "importance", ascending=False
    )


def _permutation_importance(
    model: ExtraTreesClassifier,
    test_features: pd.DataFrame,
    all_feature_names: list[str],
    permutation_feature_names: list[str],
) -> pd.DataFrame:
    rng = np.random.default_rng(42)
    x = test_features[all_feature_names].copy()
    y = test_features["cas_subtype"]
    baseline = accuracy_score(y, model.predict(x))
    rows = []
    for feature in permutation_feature_names:
        drops = []
        original = x[feature].to_numpy(copy=True)
        for _ in range(5):
            shuffled = original.copy()
            rng.shuffle(shuffled)
            x[feature] = shuffled
            drops.append(baseline - accuracy_score(y, model.predict(x)))
        x[feature] = original
        rows.append(
            {
                "feature": feature,
                "importance_mean": float(np.mean(drops)),
                "importance_std": float(np.std(drops, ddof=0)),
                "baseline_accuracy": float(baseline),
            }
        )
    return pd.DataFrame(rows).assign(category=lambda table: table["feature"].map(_feature_category)).sort_values(
        "importance_mean", ascending=False
    )


def _permutation_feature_subset(
    builtin: pd.DataFrame,
    feature_names: list[str],
    top_n: int,
) -> list[str]:
    priority = [
        "repeat_length",
        "repeat_gc_percent",
        "repeat_at_percent",
        "repeat_gc_skew",
        "spacer_count",
        "mean_spacer_length",
        "repeat_self_rc_identity",
        "repeat_longest_inverted_stem",
        "repeat_hairpin_score",
        "array_length_estimate",
        "spacer_repeat_length_ratio",
    ]
    selected = list(builtin.head(top_n)["feature"])
    for feature in priority:
        if feature in feature_names and feature not in selected:
            selected.append(feature)
    return selected


def _category_summary(builtin: pd.DataFrame, permutation: pd.DataFrame) -> pd.DataFrame:
    built = builtin.groupby("category", as_index=False).agg(
        builtin_importance=("importance", "sum"),
        feature_count=("feature", "size"),
    )
    perm = permutation.groupby("category", as_index=False).agg(
        permutation_importance=("importance_mean", "sum")
    )
    merged = built.merge(perm, on="category", how="outer").fillna(0.0)
    return merged.sort_values("builtin_importance", ascending=False)


def _typeiii_error_feature_summary(test_features: pd.DataFrame) -> pd.DataFrame:
    features = [
        "repeat_length",
        "repeat_gc_percent",
        "repeat_at_percent",
        "repeat_gc_skew",
        "spacer_count",
        "mean_spacer_length",
        "repeat_self_rc_identity",
        "repeat_longest_inverted_stem",
        "repeat_hairpin_score",
    ]
    groups = {
        "III-B correct": test_features[(test_features["cas_subtype"] == "III-B") & test_features["correct"]],
        "III-B wrong to I-B": test_features[
            (test_features["cas_subtype"] == "III-B")
            & (test_features["predicted_subtype"] == "I-B")
        ],
        "III-D correct": test_features[(test_features["cas_subtype"] == "III-D") & test_features["correct"]],
        "III-D wrong to III-A/B": test_features[
            (test_features["cas_subtype"] == "III-D")
            & (test_features["predicted_subtype"].isin(["III-A", "III-B"]))
        ],
    }
    rows = []
    for group_name, subset in groups.items():
        for feature in features:
            rows.append(
                {
                    "group": group_name,
                    "feature": feature,
                    "row_count": len(subset),
                    "mean": float(subset[feature].mean()) if not subset.empty else 0.0,
                    "median": float(subset[feature].median()) if not subset.empty else 0.0,
                }
            )
    return pd.DataFrame(rows)


def _feature_category(feature: str) -> str:
    if feature.startswith("kmer_"):
        return "whole-repeat k-mer"
    if feature.startswith("terminal_kmer_"):
        return "terminal k-mer"
    if feature.startswith("repeat_start_") or feature.startswith("repeat_end_"):
        return "terminal composition"
    if feature in {"spacer_count", "mean_spacer_length", "repeat_count", "array_length_estimate", "spacer_repeat_length_ratio"}:
        return "array statistics"
    if feature in {"repeat_self_rc_identity", "repeat_longest_inverted_stem", "repeat_hairpin_score"}:
        return "repeat structure"
    if feature.startswith("repeat_"):
        return "repeat composition"
    return "other"


def _plot_top_features(table: pd.DataFrame, title: str, output_base: Path, top_n: int = 25) -> None:
    top = table.head(top_n).sort_values("importance")
    plt.figure(figsize=(8, 7))
    plt.barh(top["feature"], top["importance"], color="#2F6F73")
    plt.xlabel("Importance")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output_base.with_suffix(".png"), dpi=220)
    plt.savefig(output_base.with_suffix(".svg"))
    plt.close()


def _plot_category_summary(table: pd.DataFrame, output_base: Path) -> None:
    ordered = table.sort_values("builtin_importance")
    x = range(len(ordered))
    plt.figure(figsize=(9, 5))
    plt.barh(x, ordered["builtin_importance"], color="#2F6F73", label="built-in")
    plt.scatter(ordered["permutation_importance"], x, color="#9A6B2F", label="permutation")
    plt.yticks(x, ordered["category"])
    plt.xlabel("Summed importance")
    plt.title("Feature Importance by Category")
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(output_base.with_suffix(".png"), dpi=220)
    plt.savefig(output_base.with_suffix(".svg"))
    plt.close()


def _plot_typeiii_error_summary(table: pd.DataFrame, output_base: Path) -> None:
    selected = table[table["feature"].isin(["repeat_gc_percent", "mean_spacer_length", "repeat_hairpin_score", "repeat_self_rc_identity"])]
    pivot = selected.pivot(index="group", columns="feature", values="mean")
    normalized = (pivot - pivot.min()) / (pivot.max() - pivot.min()).replace(0, 1)
    plt.figure(figsize=(8, 4))
    plt.imshow(normalized, cmap="viridis", aspect="auto")
    plt.xticks(range(len(normalized.columns)), normalized.columns, rotation=35, ha="right")
    plt.yticks(range(len(normalized.index)), normalized.index)
    plt.colorbar(label="normalized mean")
    plt.title("Type III Correct vs Error Feature Pattern")
    plt.tight_layout()
    plt.savefig(output_base.with_suffix(".png"), dpi=220)
    plt.savefig(output_base.with_suffix(".svg"))
    plt.close()


if __name__ == "__main__":
    main()
