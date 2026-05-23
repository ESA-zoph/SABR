from __future__ import annotations

from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.metrics import confusion_matrix, precision_recall_fscore_support, roc_curve, auc
from sklearn.preprocessing import label_binarize

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from crispr_phage_predictor.ml.dataset import load_repeat_cas_training_table
from crispr_phage_predictor.ml.features import build_repeat_feature_table, feature_columns
from crispr_phage_predictor.ml.train_classifier import _filter_min_class_count, _split_table


FIGURES_DIR = ROOT / "docs" / "figures"
MODEL_COMPARISON = ROOT / "docs" / "model_comparison_current.csv"
BEST_DATASET = ROOT / "data" / "training" / "repeats_cas_types_augmented_vink_genbank_targeted.csv"
PREDICTIONS_CSV = ROOT / "docs" / "best_model_predictions.csv"


def main() -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    _plot_model_accuracy_comparison()
    test_features, feature_names, model = _fit_best_model()
    _plot_confusion_matrix(test_features, model, feature_names)
    _plot_per_class_f1(test_features, model, feature_names)
    _plot_multiclass_roc(test_features, model, feature_names)
    _plot_feature_importance(model, feature_names)
    predictions = _export_predictions(test_features, model, feature_names)
    _plot_error_by_subtype(predictions)
    _plot_confidence_correct_vs_wrong(predictions)
    _plot_top_errors(predictions)
    print(f"Wrote model figures to {FIGURES_DIR}")


def _plot_model_accuracy_comparison() -> None:
    table = pd.read_csv(MODEL_COMPARISON)
    subset = table[table["split_strategy"].isin(["genus_holdout_min20", "genome_holdout"])]
    labels = subset["split_strategy"] + "\n" + subset["method"]
    plt.figure(figsize=(10, 5))
    colors = ["#2F6F73" if "genus" in split else "#7A8C99" for split in subset["split_strategy"]]
    plt.bar(labels, subset["accuracy"], color=colors)
    plt.ylim(0.75, 0.95)
    plt.ylabel("Accuracy")
    plt.title("Model Accuracy Comparison")
    plt.xticks(rotation=35, ha="right")
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "model_accuracy_comparison.png", dpi=180)
    plt.close()


def _fit_best_model():
    table = load_repeat_cas_training_table(BEST_DATASET)
    table = _filter_min_class_count(table, min_class_count=20)
    train, test = _split_table(
        table,
        test_size=0.25,
        random_state=42,
        split_strategy="group_holdout",
        group_column="genus",
    )
    train_features = build_repeat_feature_table(train)
    test_features = build_repeat_feature_table(test)
    feature_names = feature_columns(train_features)
    x_train = train_features[feature_names]
    y_train = train_features["cas_subtype"]
    model = ExtraTreesClassifier(
        n_estimators=400,
        random_state=42,
        class_weight="balanced",
        n_jobs=-1,
    )
    model.fit(x_train, y_train)
    aligned_test = test_features.reindex(columns=feature_names, fill_value=0.0)
    metadata_columns = [
        "genome_id",
        "contig_id",
        "organism",
        "repeat_sequence",
        "repeat_length",
        "spacer_count",
        "cas_subtype",
    ]
    for column in metadata_columns:
        aligned_test[column] = list(test[column] if column in test.columns else test_features[column])
    return aligned_test, feature_names, model


def _plot_confusion_matrix(test_features: pd.DataFrame, model: ExtraTreesClassifier, feature_names: list[str]) -> None:
    true_labels = list(test_features["cas_subtype"])
    predictions = list(model.predict(test_features[feature_names]))
    labels = sorted(set(true_labels) | set(predictions))
    matrix = confusion_matrix(true_labels, predictions, labels=labels)
    plt.figure(figsize=(9, 8))
    plt.imshow(matrix, cmap="Blues")
    plt.title("Best Model Confusion Matrix")
    plt.colorbar(label="Rows")
    plt.xticks(range(len(labels)), labels, rotation=45, ha="right")
    plt.yticks(range(len(labels)), labels)
    plt.xlabel("Predicted subtype")
    plt.ylabel("True subtype")
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "best_model_confusion_matrix.png", dpi=180)
    plt.close()


def _plot_per_class_f1(test_features: pd.DataFrame, model: ExtraTreesClassifier, feature_names: list[str]) -> None:
    true_labels = list(test_features["cas_subtype"])
    predictions = list(model.predict(test_features[feature_names]))
    labels = sorted(set(true_labels) | set(predictions))
    _, _, f1_scores, supports = precision_recall_fscore_support(
        true_labels,
        predictions,
        labels=labels,
        zero_division=0,
    )
    order = sorted(range(len(labels)), key=lambda index: f1_scores[index])
    ordered_labels = [labels[index] for index in order]
    ordered_f1 = [f1_scores[index] for index in order]
    ordered_supports = [supports[index] for index in order]
    plt.figure(figsize=(8, max(4, 0.32 * len(labels))))
    plt.barh(ordered_labels, ordered_f1, color="#2F6F73")
    for y_index, support in enumerate(ordered_supports):
        plt.text(0.02, y_index, f"n={support}", va="center", color="white", fontsize=8)
    plt.xlim(0, 1)
    plt.xlabel("F1 score")
    plt.ylabel("Subtype")
    plt.title("Best Model Per-Class F1")
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "best_model_per_class_f1.png", dpi=180)
    plt.close()


def _plot_multiclass_roc(test_features: pd.DataFrame, model: ExtraTreesClassifier, feature_names: list[str]) -> None:
    classes = list(model.classes_)
    probabilities = model.predict_proba(test_features[feature_names])
    y_true = label_binarize(test_features["cas_subtype"], classes=classes)
    plt.figure(figsize=(8, 6))
    plotted = 0
    for index, subtype in enumerate(classes):
        positives = y_true[:, index].sum()
        if positives == 0:
            continue
        false_positive_rate, true_positive_rate, _ = roc_curve(y_true[:, index], probabilities[:, index])
        roc_auc = auc(false_positive_rate, true_positive_rate)
        if positives >= 5 or subtype in {"V-A", "VI-B1"}:
            plt.plot(false_positive_rate, true_positive_rate, label=f"{subtype} AUC={roc_auc:.2f}")
            plotted += 1
    plt.plot([0, 1], [0, 1], linestyle="--", color="#777777", linewidth=1)
    plt.xlabel("False positive rate")
    plt.ylabel("True positive rate")
    plt.title("Best Model One-vs-Rest ROC Curves")
    if plotted:
        plt.legend(fontsize=8, loc="lower right")
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "best_model_roc_curve.png", dpi=180)
    plt.close()


def _plot_feature_importance(model: ExtraTreesClassifier, feature_names: list[str], top_n: int = 25) -> None:
    importances = pd.Series(model.feature_importances_, index=feature_names).sort_values(ascending=False)
    top = importances.head(top_n).sort_values()
    plt.figure(figsize=(8, 7))
    plt.barh(top.index, top.values, color="#2F6F73")
    plt.xlabel("Feature importance")
    plt.ylabel("Feature")
    plt.title("Best Model Top Feature Importances")
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "best_model_feature_importance.png", dpi=180)
    plt.close()


def _export_predictions(
    test_features: pd.DataFrame,
    model: ExtraTreesClassifier,
    feature_names: list[str],
) -> pd.DataFrame:
    probabilities = model.predict_proba(test_features[feature_names])
    predictions = model.classes_[probabilities.argmax(axis=1)]
    confidences = probabilities.max(axis=1)
    rows = []
    for index, (_, row) in enumerate(test_features.iterrows()):
        true_label = row["cas_subtype"]
        predicted = predictions[index]
        rows.append(
            {
                "genome_id": row.get("genome_id", ""),
                "contig_id": row.get("contig_id", ""),
                "organism": row.get("organism", ""),
                "cas_subtype": true_label,
                "predicted_subtype": predicted,
                "prediction_confidence": round(float(confidences[index]), 6),
                "correct": bool(true_label == predicted),
                "repeat_sequence": row.get("repeat_sequence", ""),
                "repeat_length": row.get("repeat_length", ""),
                "spacer_count": row.get("spacer_count", ""),
            }
        )
    output = pd.DataFrame(rows)
    output.to_csv(PREDICTIONS_CSV, index=False)
    return output


def _plot_error_by_subtype(predictions: pd.DataFrame) -> None:
    summary = (
        predictions.groupby("cas_subtype")
        .agg(total=("correct", "size"), errors=("correct", lambda values: int((~values).sum())))
        .reset_index()
    )
    summary["error_rate"] = summary["errors"] / summary["total"]
    summary = summary.sort_values("error_rate")
    plt.figure(figsize=(8, max(4, 0.32 * len(summary))))
    plt.barh(summary["cas_subtype"], summary["error_rate"], color="#8F3D3D")
    for y_index, row in summary.reset_index(drop=True).iterrows():
        plt.text(
            0.01,
            y_index,
            f"{int(row['errors'])}/{int(row['total'])}",
            va="center",
            color="white",
            fontsize=8,
        )
    plt.xlim(0, 1)
    plt.xlabel("Error rate")
    plt.ylabel("Subtype")
    plt.title("Best Model Error Rate by Subtype")
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "best_model_error_by_subtype.png", dpi=180)
    plt.close()


def _plot_confidence_correct_vs_wrong(predictions: pd.DataFrame) -> None:
    correct = predictions[predictions["correct"]]["prediction_confidence"]
    wrong = predictions[~predictions["correct"]]["prediction_confidence"]
    plt.figure(figsize=(7, 4))
    plt.hist(correct, bins=20, alpha=0.75, label="Correct", color="#2F6F73")
    plt.hist(wrong, bins=20, alpha=0.75, label="Wrong", color="#A65F2B")
    plt.xlabel("Prediction confidence")
    plt.ylabel("Rows")
    plt.title("Prediction Confidence: Correct vs Wrong")
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "best_model_confidence_correct_vs_wrong.png", dpi=180)
    plt.close()


def _plot_top_errors(predictions: pd.DataFrame, top_n: int = 12) -> None:
    errors = predictions[~predictions["correct"]]
    if errors.empty:
        return
    counts = (
        errors.groupby(["cas_subtype", "predicted_subtype"])
        .size()
        .sort_values(ascending=False)
        .head(top_n)
    )
    labels = [f"{true} -> {predicted}" for true, predicted in counts.index]
    plt.figure(figsize=(8, max(4, 0.35 * len(labels))))
    plt.barh(labels[::-1], counts.values[::-1], color="#8F3D3D")
    plt.xlabel("Error count")
    plt.ylabel("True -> predicted")
    plt.title("Most Common Best-Model Errors")
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "best_model_top_errors.png", dpi=180)
    plt.close()


if __name__ == "__main__":
    main()
