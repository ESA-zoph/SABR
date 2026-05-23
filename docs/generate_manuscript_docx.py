from __future__ import annotations

from pathlib import Path

import pandas as pd
from docx import Document
from docx.shared import Inches


ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
ASSETS = DOCS / "manuscript_assets"
OUTPUT = DOCS / "SABR_manuscript.docx"


def main() -> None:
    document = Document()
    document.add_heading("SABR: Spacer Alignment-Based Recognition for CRISPR-Phage Targeting Evidence", 0)
    document.add_paragraph(
        "A transparent workflow for CRISPR spacer-targeting evidence, PAM/PFS diagnostics, "
        "and repeat-derived CRISPR-Cas subtype prediction."
    )

    _section(
        document,
        "Abstract",
        "SABR is an early-stage bioinformatics workflow that accepts bacterial and phage FASTA files, "
        "detects candidate CRISPR arrays, extracts spacers and repeats, searches spacers against phage "
        "genomes, evaluates PAM/PFS support, predicts likely CRISPR-Cas subtype from repeat/array "
        "features, and reports bacteria-by-phage CRISPR targeting evidence. SABR is intentionally "
        "framed as an evidence mapper rather than a direct resistance caller. The best documented "
        "current production-candidate model was a flat ExtraTrees classifier with 0.9152 genus-holdout "
        "accuracy. A CRISPRCasdb-trained experimental ExtraTrees model achieved 0.9455 internal "
        "genome-holdout accuracy and transferred well to the current SABR table, but it remains a "
        "computational-label model requiring independent curated validation before final replacement."
    )

    _section(
        document,
        "Scientific Framing",
        "Spacer-protospacer matches are evidence of possible CRISPR targeting or prior exposure, not "
        "proof of biological resistance. Resistance also depends on subtype, PAM/PFS compatibility, "
        "seed conservation, functional cas genes, expression, phage escape, anti-CRISPR genes, and "
        "assembly quality. SABR therefore reports CRISPR targeting evidence, not confirmed resistance."
    )

    _section(
        document,
        "Workflow",
        "The pipeline parses bacterial and phage FASTA files, reports sequence diagnostics, detects "
        "candidate CRISPR arrays, extracts repeats and spacers, matches spacers against phage genomes, "
        "predicts Cas subtype from repeat/array features, evaluates PAM/PFS support, summarizes seed "
        "mismatches, and exports evidence matrices and detailed tables."
    )

    _add_figure(document, ASSETS / "benchmark_score_summary.png", "Figure 1. Benchmark CRISPR targeting evidence scores.")

    _section(
        document,
        "Cas Subtype Model Development",
        "SABR uses features derivable from uploaded FASTA files, including repeat length, base "
        "composition, k-mer composition, terminal features, spacer count, mean spacer length, and "
        "hairpin-like repeat features. Runtime prediction does not use organism name, taxonomy, "
        "accession, source database, or cas-gene annotation as model features."
    )

    _add_figure(document, ASSETS / "model_comparison_summary.png", "Figure 2. Model comparison across current and CRISPRCasdb-derived datasets.")

    document.add_heading("Model Results", level=1)
    comparison = pd.read_csv(DOCS / "crisprcasdb_augmented_model_comparison.csv")
    _add_table(
        document,
        comparison[["dataset", "method", "group_column", "accuracy", "notes"]].tail(8),
        "Table 1. Main model comparisons and transfer experiments.",
    )

    _section(
        document,
        "CRISPRCasdb Experiments",
        "CRISPRCasdb release 34 was imported as both direct-repeat inventories and SQL-derived "
        "candidate repeat/Cas labels. The SQL importer linked CRISPR loci to nearest same-sequence "
        "Cas clusters and produced 23,507 computational candidate rows. CRISPRCasdb-only training "
        "performed strongly, but the labels are computational candidates rather than manually curated "
        "gold-standard truth."
    )

    _add_figure(document, ASSETS / "typeiii_performance_summary.png", "Figure 3. Type III subtype F1 scores across experiments.")
    _add_figure(document, ASSETS / "calibration_comparison.png", "Figure 4. Probability calibration comparison.")
    _add_figure(document, ASSETS / "current_plus_crisprcasdb_balanced_pca_by_subtype.png", "Figure 5. PCA projection colored by subtype.")
    _add_figure(document, ASSETS / "current_plus_crisprcasdb_balanced_tsne_by_subtype.png", "Figure 6. Sampled t-SNE projection colored by subtype.")
    _add_figure(document, ASSETS / "current_plus_crisprcasdb_balanced_pca_by_dataset_group.png", "Figure 7. PCA projection showing current rows and CRISPRCasdb additions.")
    _add_figure(document, ASSETS / "current_plus_crisprcasdb_typeiii_tsne_by_dataset_group.png", "Figure 8. Type III-focused sampled t-SNE by dataset source.")
    _add_figure(document, ASSETS / "crisprcasdb_model_tsne_correct_vs_wrong.png", "Figure 9. Held-out CRISPRCasdb model t-SNE projection highlighting wrong calls.")
    _add_figure(document, ASSETS / "crisprcasdb_model_tsne_errors_by_true_subtype.png", "Figure 10. Wrong-call t-SNE projection colored by true subtype.")
    _add_figure(document, ASSETS / "crisprcasdb_model_top_error_pairs.png", "Figure 11. Most frequent wrong-call pairs for the CRISPRCasdb-trained model.")

    _section(
        document,
        "Model Interpretability",
        "The CRISPRCasdb-trained ExtraTrees model is not a neural-network black box: it can be "
        "interrogated through built-in tree impurity importance, held-out permutation tests, feature "
        "category summaries, and targeted error analyses. Built-in importance is distributed across "
        "repeat length, spacer/repeat length ratio, spacer-length statistics, terminal repeat k-mers, "
        "GC/AT composition, and hairpin-like repeat features. Category-level importance is dominated "
        "by whole-repeat k-mers and terminal k-mers, with smaller but biologically interpretable "
        "contributions from terminal composition, array statistics, repeat composition, and repeat "
        "structure. Held-out permutation drops are small because many repeat-derived features are "
        "correlated, so these results should be read as supporting evidence for feature families rather "
        "than as proof that any single nucleotide motif is causal."
    )
    _add_figure(document, ASSETS / "crisprcasdb_builtin_feature_importance.png", "Figure 12. Built-in feature importance for the CRISPRCasdb-trained ExtraTrees model.")
    _add_figure(document, ASSETS / "crisprcasdb_permutation_feature_importance.png", "Figure 13. Held-out permutation importance for selected high-priority features.")
    _add_figure(document, ASSETS / "crisprcasdb_feature_category_importance.png", "Figure 14. Feature importance summarized by biological feature category.")
    _add_figure(document, ASSETS / "crisprcasdb_typeiii_error_feature_summary.png", "Figure 15. Type III correct and wrong-call feature pattern summary.")

    _section(
        document,
        "Current Interpretation",
        "The strongest current result is that CRISPRCasdb-derived training data are highly useful for "
        "repeat-to-subtype learning. However, because CRISPRCasdb labels are computationally derived "
        "and overlap conceptually with some existing sources, final claims require independent "
        "curated/literature/CCTyper-supported validation. Error projection shows that wrong calls "
        "are not randomly distributed: they concentrate around Type III and adjacent Type I-B/I-C/I-A "
        "regions, especially III-B to I-B and III-D to III-A/III-B confusions. SABR's core contribution "
        "remains the integrated, cautious, reproducible evidence framework."
    )

    _section(
        document,
        "Limitations",
        "The internal CRISPR detector is an exact-repeat baseline. Spacer matching is exact by default. "
        "PAM/PFS rules are incomplete. CRISPRCasdb candidate labels are not gold-standard labels. "
        "Resistance/sensitivity claims remain limited by the small curated benchmark panel. Type III-B "
        "and Type III-D remain challenging in several experiments."
    )

    _section(
        document,
        "Future Work",
        "Next steps include independent CCTyper/literature validation, larger bacteria-phage benchmark "
        "panels, confidence calibration in the GUI, anti-CRISPR and PAM-failure controls, and a careful "
        "decision about replacing the runtime model with the CRISPRCasdb-trained artifact."
    )

    document.add_heading("SVG Assets", level=1)
    document.add_paragraph(
        "Manuscript figures are embedded as PNG for Word compatibility. Matching SVG files are stored "
        "in docs/manuscript_assets/ for publication-quality editing."
    )
    document.save(OUTPUT)
    print(f"Wrote {OUTPUT}")


def _section(document: Document, heading: str, text: str) -> None:
    document.add_heading(heading, level=1)
    document.add_paragraph(text)


def _add_figure(document: Document, path: Path, caption: str) -> None:
    if path.exists():
        document.add_picture(str(path), width=Inches(6.2))
        document.add_paragraph(caption)


def _add_table(document: Document, table: pd.DataFrame, caption: str) -> None:
    document.add_paragraph(caption)
    doc_table = document.add_table(rows=1, cols=len(table.columns))
    for index, column in enumerate(table.columns):
        doc_table.rows[0].cells[index].text = str(column)
    for _, row in table.iterrows():
        cells = doc_table.add_row().cells
        for index, column in enumerate(table.columns):
            value = row[column]
            cells[index].text = f"{value:.4f}" if isinstance(value, float) else str(value)


if __name__ == "__main__":
    main()
