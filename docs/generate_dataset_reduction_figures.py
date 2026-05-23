from __future__ import annotations

from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from crispr_phage_predictor.io import sequence_hash
from crispr_phage_predictor.ml.dataset import load_repeat_cas_training_table
from crispr_phage_predictor.ml.features import build_repeat_feature_table, feature_columns


CURRENT_DATASET = ROOT / "data" / "training" / "repeats_cas_types_augmented_vink_genbank_targeted.csv"
ADDITIONS_DATASET = ROOT / "data" / "training" / "repeats_cas_types_crisprcasdb_sql_balanced_additions.csv"
TYPEIII_ADDITIONS_DATASET = ROOT / "data" / "training" / "repeats_cas_types_crisprcasdb_typeiii_balanced_additions.csv"
OUTPUT_DIR = ROOT / "docs" / "manuscript_assets"


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    current = _load_with_group(CURRENT_DATASET, "current_training")
    additions = _load_with_group(ADDITIONS_DATASET, "crisprcasdb_balanced_addition")
    typeiii_additions = _load_with_group(TYPEIII_ADDITIONS_DATASET, "crisprcasdb_typeiii_addition")

    combined = pd.concat([current, additions], ignore_index=True)
    typeiii_combined = pd.concat([current, typeiii_additions], ignore_index=True)

    _write_projection_bundle(combined, "current_plus_crisprcasdb_balanced")
    _write_projection_bundle(typeiii_combined, "current_plus_crisprcasdb_typeiii")
    print(f"Wrote dimensionality-reduction figures and tables to {OUTPUT_DIR}")


def _load_with_group(path: Path, group: str) -> pd.DataFrame:
    table = load_repeat_cas_training_table(path)
    table = table.copy()
    table["dataset_group"] = group
    table["repeat_hash"] = table["repeat_sequence"].map(sequence_hash)
    return table


def _write_projection_bundle(table: pd.DataFrame, prefix: str) -> None:
    features = build_repeat_feature_table(table)
    columns = feature_columns(features)
    scaled = StandardScaler().fit_transform(features[columns])

    pca = PCA(n_components=2, random_state=42)
    pca_coords = pca.fit_transform(scaled)
    pca_table = _projection_table(table, pca_coords, "PCA1", "PCA2")
    pca_table.to_csv(OUTPUT_DIR / f"{prefix}_pca_coordinates.csv", index=False)
    _plot_projection(
        pca_table,
        x_col="PCA1",
        y_col="PCA2",
        title=f"{prefix}: PCA of repeat/array features",
        output_path=OUTPUT_DIR / f"{prefix}_pca_by_subtype.png",
        color_col="cas_subtype",
    )
    _plot_projection(
        pca_table,
        x_col="PCA1",
        y_col="PCA2",
        title=f"{prefix}: PCA by dataset source",
        output_path=OUTPUT_DIR / f"{prefix}_pca_by_dataset_group.png",
        color_col="dataset_group",
    )

    sampled = _sample_for_tsne(table, max_rows=2500)
    sampled_features = build_repeat_feature_table(sampled)
    sampled_columns = feature_columns(sampled_features)
    sampled_scaled = StandardScaler().fit_transform(sampled_features[sampled_columns])
    tsne = TSNE(
        n_components=2,
        random_state=42,
        init="pca",
        learning_rate="auto",
        perplexity=35,
        max_iter=1000,
    )
    tsne_coords = tsne.fit_transform(sampled_scaled)
    tsne_table = _projection_table(sampled, tsne_coords, "TSNE1", "TSNE2")
    tsne_table.to_csv(OUTPUT_DIR / f"{prefix}_tsne_coordinates_sampled.csv", index=False)
    _plot_projection(
        tsne_table,
        x_col="TSNE1",
        y_col="TSNE2",
        title=f"{prefix}: sampled t-SNE of repeat/array features",
        output_path=OUTPUT_DIR / f"{prefix}_tsne_by_subtype.png",
        color_col="cas_subtype",
    )
    _plot_projection(
        tsne_table,
        x_col="TSNE1",
        y_col="TSNE2",
        title=f"{prefix}: sampled t-SNE by dataset source",
        output_path=OUTPUT_DIR / f"{prefix}_tsne_by_dataset_group.png",
        color_col="dataset_group",
    )


def _projection_table(table: pd.DataFrame, coords, x_col: str, y_col: str) -> pd.DataFrame:
    metadata = table[
        [
            "dataset_group",
            "source",
            "genome_id",
            "organism",
            "repeat_hash",
            "repeat_sequence",
            "cas_type",
            "cas_subtype",
            "label_confidence",
        ]
    ].reset_index(drop=True)
    metadata[x_col] = coords[:, 0]
    metadata[y_col] = coords[:, 1]
    return metadata


def _sample_for_tsne(table: pd.DataFrame, max_rows: int) -> pd.DataFrame:
    if len(table) <= max_rows:
        return table.copy()
    parts = []
    for _, subtype_table in table.groupby("cas_subtype", sort=True):
        fraction = len(subtype_table) / len(table)
        subtype_n = max(10, int(round(max_rows * fraction)))
        parts.append(subtype_table.sample(n=min(len(subtype_table), subtype_n), random_state=42))
    sampled = pd.concat(parts, ignore_index=True)
    if len(sampled) > max_rows:
        sampled = sampled.sample(n=max_rows, random_state=42)
    return sampled.reset_index(drop=True)


def _plot_projection(
    table: pd.DataFrame,
    x_col: str,
    y_col: str,
    title: str,
    output_path: Path,
    color_col: str,
) -> None:
    groups = list(table[color_col].dropna().astype(str).sort_values().unique())
    cmap = plt.get_cmap("tab20")
    plt.figure(figsize=(9, 7))
    for index, group in enumerate(groups):
        subset = table[table[color_col].astype(str) == group]
        plt.scatter(
            subset[x_col],
            subset[y_col],
            s=12,
            alpha=0.72,
            label=group,
            color=cmap(index % 20),
            edgecolors="none",
        )
    plt.title(title)
    plt.xlabel(x_col)
    plt.ylabel(y_col)
    plt.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=8, frameon=False)
    plt.tight_layout()
    plt.savefig(output_path, dpi=180)
    plt.savefig(output_path.with_suffix(".svg"))
    plt.close()


if __name__ == "__main__":
    main()
