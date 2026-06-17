from __future__ import annotations

from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from crispr_phage_predictor.ml.dataset import load_repeat_cas_training_table
from crispr_phage_predictor.ml.features import build_repeat_feature_table, feature_columns
from crispr_phage_predictor.ml.train_classifier import _filter_min_class_count, _split_table


DATASET = ROOT / "data" / "training" / "repeats_cas_types_crisprcasdb_sql_candidate.csv"
OUTPUT_DIR = ROOT / "docs" / "manuscript_assets"


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    predictions = _heldout_predictions()
    predictions.to_csv(OUTPUT_DIR / "selected_model_error_projection_coordinates.csv", index=False)
    _plot_projection(
        predictions,
        "PCA1",
        "PCA2",
        "correct",
        OUTPUT_DIR / "selected_model_pca_correct_vs_wrong",
        "Selected SABR subtype model: PCA correct vs wrong calls",
    )
    _plot_projection(
        predictions,
        "TSNE1",
        "TSNE2",
        "correct",
        OUTPUT_DIR / "selected_model_tsne_correct_vs_wrong",
        "Selected SABR subtype model: t-SNE correct vs wrong calls",
    )
    errors = predictions[~predictions["correct"]].copy()
    errors.to_csv(OUTPUT_DIR / "selected_model_error_rows.csv", index=False)
    _plot_projection(
        errors,
        "TSNE1",
        "TSNE2",
        "true_subtype",
        OUTPUT_DIR / "selected_model_tsne_errors_by_true_subtype",
        "Wrong calls clustered by true subtype",
    )
    _plot_error_pair_counts(errors)
    print(f"Wrote error projection figures to {OUTPUT_DIR}")


def _heldout_predictions() -> pd.DataFrame:
    table = load_repeat_cas_training_table(DATASET)
    table = _filter_min_class_count(table, min_class_count=20)
    train, test = _split_table(
        table,
        test_size=0.25,
        random_state=42,
        split_strategy="group_holdout",
        group_column="genome_id",
    )
    train_features = build_repeat_feature_table(train)
    test_features = build_repeat_feature_table(test)
    features = feature_columns(train_features)
    aligned_test = test_features.reindex(columns=features, fill_value=0.0)
    model = ExtraTreesClassifier(
        n_estimators=400,
        random_state=42,
        class_weight="balanced",
        n_jobs=-1,
    )
    model.fit(train_features[features], train_features["cas_subtype"])
    predicted = model.predict(aligned_test)
    probabilities = model.predict_proba(aligned_test)
    confidence = probabilities.max(axis=1)

    scaled = StandardScaler().fit_transform(aligned_test[features])
    pca_coords = PCA(n_components=2, random_state=42).fit_transform(scaled)
    tsne_sample_index = _sample_indices(test, max_rows=3000)
    sampled_scaled = scaled[tsne_sample_index]
    tsne_coords_sampled = TSNE(
        n_components=2,
        random_state=42,
        init="pca",
        learning_rate="auto",
        perplexity=35,
        max_iter=1000,
    ).fit_transform(sampled_scaled)

    output = test.reset_index(drop=True)[
        ["source", "genome_id", "organism", "repeat_sequence", "cas_type", "cas_subtype"]
    ].copy()
    output = output.rename(columns={"cas_subtype": "true_subtype"})
    output["predicted_subtype"] = predicted
    output["confidence"] = confidence
    output["correct"] = output["true_subtype"] == output["predicted_subtype"]
    output["PCA1"] = pca_coords[:, 0]
    output["PCA2"] = pca_coords[:, 1]
    output["TSNE1"] = pd.NA
    output["TSNE2"] = pd.NA
    output.loc[tsne_sample_index, "TSNE1"] = tsne_coords_sampled[:, 0]
    output.loc[tsne_sample_index, "TSNE2"] = tsne_coords_sampled[:, 1]
    output["error_pair"] = output["true_subtype"] + " -> " + output["predicted_subtype"]
    return output


def _sample_indices(table: pd.DataFrame, max_rows: int) -> list[int]:
    if len(table) <= max_rows:
        return list(range(len(table)))
    sampled = []
    reset = table.reset_index(drop=True)
    for _, subtype_table in reset.groupby("cas_subtype", sort=True):
        fraction = len(subtype_table) / len(reset)
        subtype_n = max(10, int(round(max_rows * fraction)))
        sampled.extend(
            subtype_table.sample(n=min(len(subtype_table), subtype_n), random_state=42).index.tolist()
        )
    if len(sampled) > max_rows:
        sampled = list(pd.Series(sampled).sample(n=max_rows, random_state=42))
    return sorted(sampled)


def _plot_projection(table: pd.DataFrame, x_col: str, y_col: str, color_col: str, output_base: Path, title: str) -> None:
    plot_table = table.dropna(subset=[x_col, y_col]).copy()
    if plot_table.empty:
        return
    groups = list(plot_table[color_col].astype(str).sort_values().unique())
    cmap = plt.get_cmap("tab20")
    plt.figure(figsize=(9, 7))
    for index, group in enumerate(groups):
        subset = plot_table[plot_table[color_col].astype(str) == group]
        size = 26 if group == "False" else 11
        alpha = 0.9 if group == "False" else 0.45
        plt.scatter(
            subset[x_col],
            subset[y_col],
            s=size,
            alpha=alpha,
            label=group,
            color="#C0443E" if group == "False" else cmap(index % 20),
            edgecolors="none",
        )
    plt.title(title)
    plt.xlabel(x_col)
    plt.ylabel(y_col)
    plt.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=8, frameon=False)
    plt.tight_layout()
    plt.savefig(output_base.with_suffix(".png"), dpi=220)
    plt.savefig(output_base.with_suffix(".svg"))
    plt.close()


def _plot_error_pair_counts(errors: pd.DataFrame) -> None:
    counts = errors["error_pair"].value_counts().head(15).sort_values()
    plt.figure(figsize=(8, 5))
    plt.barh(counts.index, counts.values, color="#C0443E")
    plt.xlabel("Wrong calls")
    plt.title("Selected SABR Model: Most Frequent Confusions")
    plt.tight_layout()
    plt.savefig((OUTPUT_DIR / "selected_model_top_error_pairs").with_suffix(".png"), dpi=220)
    plt.savefig((OUTPUT_DIR / "selected_model_top_error_pairs").with_suffix(".svg"))
    plt.close()


if __name__ == "__main__":
    main()
