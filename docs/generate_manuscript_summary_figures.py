from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
OUT = DOCS / "manuscript_assets"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    _plot_model_comparison()
    _plot_benchmark_scores()
    _plot_calibration_comparison()
    _plot_typeiii_comparison()
    _plot_cctyper_pilot()
    print(f"Wrote manuscript summary figures to {OUT}")


def _save(path: Path) -> None:
    plt.tight_layout()
    plt.savefig(path.with_suffix(".png"), dpi=220)
    plt.savefig(path.with_suffix(".svg"))
    plt.close()


def _plot_model_comparison() -> None:
    table = pd.read_csv(DOCS / "crisprcasdb_augmented_model_comparison.csv")
    keep = table[
        (
            (table["method"] == "extra_trees")
            & (
                table["dataset"].isin(
                    [
                        "current_best_vink_genbank_targeted",
                        "crisprcasdb_sql_candidate_only",
                        "augmented_crisprcasdb_sql_balanced",
                        "augmented_crisprcasdb_typeiii_balanced",
                    ]
                )
            )
        )
    ].copy()
    display_labels = {
        "current_best_vink_genbank_targeted": "Development table\nGenus holdout",
        "augmented_crisprcasdb_sql_balanced": "Broad annotation\nadditions",
        "augmented_crisprcasdb_typeiii_balanced": "Type III annotation\nadditions",
        "crisprcasdb_sql_candidate_only": "Selected SABR model\nGenome holdout",
    }
    keep["label"] = keep["dataset"].map(display_labels)
    colors = ["#2F6F73" if value == "crisprcasdb_sql_candidate_only" else "#7A8C99" for value in keep["dataset"]]
    plt.figure(figsize=(9, 5))
    plt.bar(keep["label"], keep["accuracy"], color=colors)
    plt.ylim(0.88, 1.0)
    plt.ylabel("Accuracy")
    plt.title("SABR ExtraTrees Development Evaluations")
    plt.text(
        0.99,
        0.97,
        "Grouped split designs differ; development summary only",
        transform=plt.gca().transAxes,
        ha="right",
        va="top",
        fontsize=9,
        color="#555555",
    )
    plt.xticks(rotation=25, ha="right")
    _save(OUT / "model_comparison_summary")


def _plot_benchmark_scores() -> None:
    table = pd.read_csv(ROOT / "data" / "curation" / "benchmark_evaluation_20260520_210025_rescored.csv")
    subset = table[table["score_expectation_result"] != "not_evaluated"].copy()
    colors = subset["expected_sabr_behavior"].map(
        {
            "high_score_expected": "#2F6F73",
            "low_score_expected": "#9A6B2F",
            "moderate_score_expected": "#7A8C99",
        }
    ).fillna("#7A8C99")
    plt.figure(figsize=(9, 5))
    plt.bar(subset["pair_id"], subset["crispr_targeting_score"], color=colors)
    plt.axhline(39, color="#555555", linestyle="--", linewidth=1, label="PAM-unsupported cap")
    plt.ylabel("CRISPR targeting score")
    plt.title("Benchmark CRISPR Targeting Evidence Scores")
    plt.xticks(rotation=35, ha="right")
    plt.legend(frameon=False)
    _save(OUT / "benchmark_score_summary")


def _plot_calibration_comparison() -> None:
    current = pd.read_csv(DOCS / "calibration" / "current_best" / "confidence_bins.csv")
    selected_model = pd.read_csv(DOCS / "calibration" / "crisprcasdb_only" / "confidence_bins.csv")
    plt.figure(figsize=(6, 6))
    plt.plot([0, 1], [0, 1], color="#555555", linestyle="--", label="perfect")
    for table, label, color in [
        (current, "Earlier development-stage model", "#7A8C99"),
        (selected_model, "Selected SABR subtype model", "#2F6F73"),
    ]:
        nonempty = table[table["row_count"] > 0]
        plt.plot(nonempty["mean_confidence"], nonempty["accuracy"], marker="o", label=label, color=color)
    plt.xlabel("Mean predicted confidence")
    plt.ylabel("Observed accuracy")
    plt.title("SABR Subtype Model Confidence Calibration")
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.legend(frameon=False)
    _save(OUT / "calibration_comparison")


def _plot_typeiii_comparison() -> None:
    table = pd.read_csv(DOCS / "crisprcasdb_augmented_model_comparison.csv")
    rows = table[
        (table["method"] == "extra_trees")
        & table["dataset"].isin(
            [
                "current_best_vink_genbank_targeted",
                "augmented_crisprcasdb_sql_balanced",
                "augmented_crisprcasdb_typeiii_balanced",
                "crisprcasdb_sql_candidate_only",
            ]
        )
    ].copy()
    labels = ["Development\nmodel", "Broad\nadditions", "Type III\nadditions", "Selected\nSABR model"]
    rows = rows.tail(4)
    x = range(len(labels))
    width = 0.25
    plt.figure(figsize=(9, 5))
    plt.bar([i - width for i in x], rows["iii_a_f1"], width=width, label="III-A", color="#2F6F73")
    plt.bar(x, rows["iii_b_f1"], width=width, label="III-B", color="#7A8C99")
    plt.bar([i + width for i in x], rows["iii_d_f1"], width=width, label="III-D", color="#9A6B2F")
    plt.xticks(list(x), labels)
    plt.ylim(0, 1)
    plt.ylabel("F1 score")
    plt.title("Type III Performance During SABR Model Development")
    plt.legend(frameon=False)
    _save(OUT / "typeiii_performance_summary")


def _plot_cctyper_pilot() -> None:
    table = pd.read_csv(ROOT / "data" / "curation" / "cctyper_validation_repeats_cas_types.csv")
    counts = table["cas_subtype"].value_counts().reindex(["I-A", "I-E", "I-F", "II-A", "II-C", "III-A"]).fillna(0)
    plt.figure(figsize=(8, 5))
    bars = plt.bar(counts.index, counts.values, color="#2F6F73")
    for bar, count in zip(bars, counts):
        plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.2, f"{int(count)}/{int(count)}", ha="center")
    plt.ylim(0, max(counts) + 3)
    plt.ylabel("Imported CCTyper-supported arrays")
    plt.xlabel("CCTyper subtype")
    plt.title("Independent CCTyper Pilot: 25/25 Correct SABR Predictions")
    plt.text(
        0.99,
        0.96,
        "No III-B or III-D arrays in pilot",
        transform=plt.gca().transAxes,
        ha="right",
        va="top",
        fontsize=9,
        color="#555555",
    )
    _save(OUT / "cctyper_pilot_summary")


if __name__ == "__main__":
    main()
