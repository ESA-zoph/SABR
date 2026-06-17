from __future__ import annotations

from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, auc, confusion_matrix, roc_auc_score, roc_curve
from sklearn.multiclass import OneVsRestClassifier
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler, label_binarize

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from crispr_phage_predictor.ml.dataset import load_repeat_cas_training_table
from crispr_phage_predictor.ml.features import build_repeat_feature_table, feature_columns
from crispr_phage_predictor.ml.train_classifier import _filter_min_class_count, _split_table


ASSETS = ROOT / "docs" / "manuscript_assets"
DEVELOPMENT_DATASET = ROOT / "data" / "training" / "repeats_cas_types_augmented_vink_genbank_targeted.csv"
SELECTED_MODEL_DATASET = ROOT / "data" / "training" / "repeats_cas_types_crisprcasdb_sql_candidate.csv"
EVIDENCE_MATRIX = ROOT / "outputs" / "runs" / "20260520_210025_rescored" / "evidence_matrix.csv"
PROJECTION_COORDINATES = ASSETS / "selected_model_error_projection_coordinates.csv"


def main() -> None:
    ASSETS.mkdir(parents=True, exist_ok=True)
    _plot_architecture_selection_roc()
    selected_test, model, features = _fit_selected_model()
    _plot_selected_model_roc(selected_test, model, features)
    _plot_selected_confusion_heatmap(selected_test, model, features)
    _plot_tsne_panels()
    _plot_targeting_score_heatmap()
    print(f"Wrote advanced manuscript figures to {ASSETS}")


def _save(name: str) -> None:
    plt.tight_layout()
    plt.savefig(ASSETS / f"{name}.png", dpi=240)
    plt.savefig(ASSETS / f"{name}.svg")
    plt.close()


def _development_split() -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    table = _filter_min_class_count(load_repeat_cas_training_table(DEVELOPMENT_DATASET), min_class_count=20)
    train, test = _split_table(
        table=table,
        test_size=0.25,
        random_state=42,
        split_strategy="group_holdout",
        group_column="genus",
    )
    train_features = build_repeat_feature_table(train)
    test_features = build_repeat_feature_table(test)
    features = feature_columns(train_features)
    return train_features, test_features.reindex(columns=features + ["cas_subtype"], fill_value=0.0), features


def _plot_architecture_selection_roc() -> None:
    train, test, features = _development_split()
    x_train = train[features]
    y_train = train["cas_subtype"]
    x_test = test[features]
    y_test = test["cas_subtype"]
    models = [
        (
            "Logistic regression",
            make_pipeline(
                StandardScaler(),
                OneVsRestClassifier(
                    LogisticRegression(
                        max_iter=1000,
                        class_weight="balanced",
                        random_state=42,
                        solver="liblinear",
                    )
                ),
            ),
            "#8996A3",
        ),
        (
            "Random forest",
            RandomForestClassifier(n_estimators=200, random_state=42, class_weight="balanced", n_jobs=-1),
            "#B07B33",
        ),
        (
            "ExtraTrees",
            ExtraTreesClassifier(n_estimators=400, random_state=42, class_weight="balanced", n_jobs=-1),
            "#21717A",
        ),
    ]
    classes = sorted(y_train.unique())
    y_binary = label_binarize(y_test, classes=classes)
    plt.figure(figsize=(8, 6))
    rows = []
    for label, model, color in models:
        model.fit(x_train, y_train)
        probabilities = model.predict_proba(x_test)
        if not np.array_equal(getattr(model, "classes_", classes), classes):
            model_classes = list(getattr(model, "classes_", classes))
            probabilities = probabilities[:, [model_classes.index(value) for value in classes]]
        fpr, tpr, _ = roc_curve(y_binary.ravel(), probabilities.ravel())
        micro_auc = auc(fpr, tpr)
        evaluable = [
            index
            for index in range(len(classes))
            if 0 < y_binary[:, index].sum() < len(y_binary)
        ]
        class_aucs = [
            roc_auc_score(y_binary[:, index], probabilities[:, index])
            for index in evaluable
        ]
        macro_auc = float(np.mean(class_aucs))
        accuracy = accuracy_score(y_test, model.predict(x_test))
        rows.append({"model": label, "accuracy": accuracy, "micro_auc": micro_auc, "macro_auc": macro_auc})
        plt.plot(fpr, tpr, linewidth=2.2, color=color, label=f"{label} (macro AUC={macro_auc:.3f})")
    plt.plot([0, 1], [0, 1], linestyle="--", color="#777777", linewidth=1)
    plt.xlabel("False positive rate")
    plt.ylabel("True positive rate")
    plt.title("Architecture Selection: One-vs-Rest ROC on Common Genus Holdout")
    plt.legend(frameon=False, loc="lower right")
    plt.text(
        0.04,
        0.10,
        "Same held-out genera; macro AUC over evaluable test subtypes",
        fontsize=9,
        color="#555555",
    )
    _save("architecture_selection_roc")
    pd.DataFrame(rows).to_csv(ASSETS / "architecture_selection_roc_metrics.csv", index=False)


def _fit_selected_model() -> tuple[pd.DataFrame, ExtraTreesClassifier, list[str]]:
    table = _filter_min_class_count(load_repeat_cas_training_table(SELECTED_MODEL_DATASET), min_class_count=20)
    train, test = _split_table(
        table=table,
        test_size=0.25,
        random_state=42,
        split_strategy="group_holdout",
        group_column="genome_id",
    )
    train_features = build_repeat_feature_table(train)
    test_features = build_repeat_feature_table(test)
    features = feature_columns(train_features)
    aligned_test = test_features.reindex(columns=features, fill_value=0.0)
    aligned_test["cas_subtype"] = list(test_features["cas_subtype"])
    model = ExtraTreesClassifier(n_estimators=400, random_state=42, class_weight="balanced", n_jobs=-1)
    model.fit(train_features[features], train_features["cas_subtype"])
    return aligned_test, model, features


def _plot_selected_model_roc(test: pd.DataFrame, model: ExtraTreesClassifier, features: list[str]) -> None:
    classes = list(model.classes_)
    probabilities = model.predict_proba(test[features])
    y_binary = label_binarize(test["cas_subtype"], classes=classes)
    plt.figure(figsize=(8.5, 6.5))
    micro_fpr, micro_tpr, _ = roc_curve(y_binary.ravel(), probabilities.ravel())
    plt.plot(
        micro_fpr,
        micro_tpr,
        color="#1F4D55",
        linewidth=3,
        label=f"Micro-average (AUC={auc(micro_fpr, micro_tpr):.3f})",
    )
    palette = {"III-A": "#C25B47", "III-B": "#9B3C46", "III-D": "#D08738", "I-B": "#617C9B", "I-C": "#3B8C87"}
    for subtype, color in palette.items():
        if subtype not in classes:
            continue
        index = classes.index(subtype)
        if y_binary[:, index].sum() == 0:
            continue
        fpr, tpr, _ = roc_curve(y_binary[:, index], probabilities[:, index])
        plt.plot(fpr, tpr, linewidth=1.8, color=color, label=f"{subtype} (AUC={auc(fpr, tpr):.3f})")
    plt.plot([0, 1], [0, 1], linestyle="--", color="#777777", linewidth=1)
    plt.xlabel("False positive rate")
    plt.ylabel("True positive rate")
    plt.title("Selected SABR Subtype Model: Genome-Held-Out ROC Curves")
    plt.legend(frameon=False, loc="lower right", fontsize=9)
    _save("selected_model_roc_focus_subtypes")


def _plot_selected_confusion_heatmap(test: pd.DataFrame, model: ExtraTreesClassifier, features: list[str]) -> None:
    true = test["cas_subtype"]
    predicted = model.predict(test[features])
    labels = sorted(set(true) | set(predicted))
    matrix = confusion_matrix(true, predicted, labels=labels, normalize="true")
    plt.figure(figsize=(10, 8))
    plt.imshow(matrix, cmap="YlGnBu", vmin=0, vmax=1)
    plt.colorbar(label="Row-normalized proportion")
    plt.xticks(range(len(labels)), labels, rotation=45, ha="right")
    plt.yticks(range(len(labels)), labels)
    for row in range(len(labels)):
        for column in range(len(labels)):
            value = matrix[row, column]
            if value >= 0.08 or row == column:
                plt.text(column, row, f"{value:.2f}", ha="center", va="center", fontsize=7, color="white" if value > 0.55 else "black")
    plt.xlabel("Predicted subtype")
    plt.ylabel("True subtype")
    plt.title("Selected SABR Model: Normalized Confusion Heatmap")
    _save("selected_model_confusion_heatmap")


def _plot_tsne_panels() -> None:
    table = pd.read_csv(PROJECTION_COORDINATES).dropna(subset=["TSNE1", "TSNE2"]).copy()
    table["cas_group"] = table["cas_type"].str.replace("Type ", "", regex=False)
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    group_colors = {"I": "#21717A", "II": "#D08738", "III": "#9B3C46", "V": "#617C9B", "VI": "#8B7BA8"}
    for group in sorted(table["cas_group"].dropna().unique()):
        subset = table[table["cas_group"] == group]
        axes[0].scatter(subset["TSNE1"], subset["TSNE2"], s=10, alpha=0.6, color=group_colors.get(group, "#777777"), label=f"Type {group}")
    axes[0].set_title("A. Broad Cas Type Grouping")
    axes[0].set_xlabel("t-SNE 1")
    axes[0].set_ylabel("t-SNE 2")
    axes[0].legend(frameon=False, markerscale=1.6)
    focus = {"III-A", "III-B", "III-D", "I-A", "I-B", "I-C"}
    focal = table[(table["true_subtype"].isin(focus)) | (table["predicted_subtype"].isin(focus))].copy()
    colors = {
        "III-A": "#C25B47", "III-B": "#9B3C46", "III-D": "#D08738",
        "I-A": "#8BA7A5", "I-B": "#21717A", "I-C": "#617C9B",
    }
    for subtype in sorted(focus):
        subset = focal[focal["true_subtype"] == subtype]
        axes[1].scatter(subset["TSNE1"], subset["TSNE2"], s=13, alpha=0.65, color=colors[subtype], label=subtype)
    errors = focal[~focal["correct"].astype(bool)]
    axes[1].scatter(errors["TSNE1"], errors["TSNE2"], s=28, facecolors="none", edgecolors="#222222", linewidths=0.65, label="Wrong call")
    axes[1].set_title("B. Type III / Adjacent Type I Error Region")
    axes[1].set_xlabel("t-SNE 1")
    axes[1].set_ylabel("t-SNE 2")
    axes[1].legend(frameon=False, fontsize=8, markerscale=1.4, ncol=2)
    fig.suptitle("Selected SABR Model: Held-Out Repeat/Array Feature Space", fontsize=15)
    plt.tight_layout()
    plt.savefig(ASSETS / "selected_model_tsne_group_and_subtype_panels.png", dpi=240)
    plt.savefig(ASSETS / "selected_model_tsne_group_and_subtype_panels.svg")
    plt.close()


def _plot_targeting_score_heatmap() -> None:
    evidence = pd.read_csv(EVIDENCE_MATRIX)
    matrix = evidence.pivot(index="bacterium", columns="phage", values="crispr_targeting_score")
    matrix = matrix.loc[:, matrix.max(axis=0) > 0]
    display = matrix.copy()
    display.index = display.index.str.replace(".fasta", "", regex=False)
    display.columns = display.columns.str.replace(".fasta", "", regex=False)
    plt.figure(figsize=(9, 7))
    plt.imshow(display.values, cmap="YlGnBu", vmin=0, vmax=100, aspect="auto")
    plt.colorbar(label="CRISPR targeting score (0-100)")
    plt.xticks(range(len(display.columns)), display.columns, rotation=45, ha="right")
    plt.yticks(range(len(display.index)), display.index)
    for row in range(display.shape[0]):
        for column in range(display.shape[1]):
            value = display.iloc[row, column]
            if value > 0:
                plt.text(column, row, f"{value:.0f}", ha="center", va="center", fontsize=8, color="white" if value >= 60 else "black")
    plt.xlabel("Phage genome")
    plt.ylabel("Bacterial genome")
    plt.title("SABR Targeting-Evidence Matrix: Non-Zero Phage Columns")
    _save("targeting_score_heatmap")


if __name__ == "__main__":
    main()
