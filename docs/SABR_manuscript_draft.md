# SABR: A Transparent CRISPR Spacer-Targeting Evidence and Repeat-Based Cas-Typing Tool

## Working Title

SABR: Spacer Alignment-Based Recognition for cautious CRISPR-phage targeting evidence and repeat-derived CRISPR-Cas subtype prediction

## Abstract

CRISPR spacer-protospacer matches provide evidence of prior phage exposure or potential CRISPR targeting, but they do not by themselves prove biological resistance. Resistance depends on additional factors including CRISPR-Cas subtype, PAM or PFS compatibility, seed integrity, cas gene functionality, array expression, phage escape mutations, and anti-CRISPR mechanisms. We developed SABR, an early-stage, transparent bioinformatics workflow that accepts bacterial and phage FASTA inputs, detects candidate CRISPR arrays, extracts spacers and repeats, searches spacers against phage genomes, evaluates candidate PAM/PFS support, predicts likely CRISPR-Cas subtype from repeat and array features, and reports a bacteria-by-phage CRISPR targeting evidence matrix.

SABR is intentionally framed as an evidence mapper rather than a direct resistance caller. The current implementation includes internal exact-repeat CRISPR detection, internal exact spacer matching, optional external MinCED-compatible and BLASTN backends, PAM/PFS diagnostics, seed-region mismatch summaries, evidence scoring, saved run artifacts, and a Streamlit graphical interface. For repeat-derived Cas subtype prediction, we compared nearest-repeat baselines and classical machine-learning models using FASTA-derived repeat and array features. The best current model is a flat ExtraTrees classifier trained on the augmented Vink/CRISPRCasdb plus targeted GenBank/MinCED dataset, achieving 0.9152 accuracy under genus-holdout validation after filtering subtypes with fewer than 20 rows. Exploratory CRISPRCasdb SQL augmentation improved Type III-A recall but reduced overall accuracy, indicating that more computational candidate rows do not automatically improve the model. These results support SABR as a reproducible framework for CRISPR targeting evidence, while emphasizing that broader curated benchmarks and experimentally grounded labels are still required before making resistance claims.

## Keywords

CRISPR, bacteriophage, spacer, protospacer, PAM, PFS, CRISPR-Cas typing, machine learning, bioinformatics, evidence scoring

## 1. Introduction

Bacteria and archaea can acquire CRISPR spacers from invading mobile genetic elements, including bacteriophages. A spacer that matches a phage protospacer can indicate previous exposure or potential targeting by a CRISPR-Cas immune system. However, spacer matches alone are insufficient to conclude resistance. Effective targeting may depend on CRISPR-Cas subtype, the presence and orientation of an appropriate PAM or PFS, seed-region conservation, functional cas genes, transcription and processing of the CRISPR array, phage escape mutations, and anti-CRISPR proteins.

Many practical workflows still treat spacer-protospacer matches as a strong proxy for host-phage interaction or resistance. SABR was developed to make this inference more cautious and transparent. The tool reports CRISPR targeting evidence and repeat-derived Cas subtype predictions, while explicitly avoiding direct claims of confirmed biological resistance unless supported by external phenotype, literature, or curated benchmark evidence.

The goals of SABR are:

1. Provide an easy-to-run FASTA-based workflow for CRISPR spacer-targeting evidence.
2. Produce bacteria-by-phage evidence matrices suitable for inspection and export.
3. Predict likely CRISPR-Cas type/subtype from repeat and array features available from uploaded FASTA files.
4. Evaluate subtype-aware PAM/PFS support where possible.
5. Keep data provenance, benchmark labels, and model limitations visible.

## 2. System Overview

SABR accepts one or more bacterial FASTA files and one or more phage FASTA files. The workflow parses and validates input records, detects candidate CRISPR arrays in bacterial genomes, extracts repeat and spacer sequences, compares spacers against phage genomes, evaluates protospacer flanking sequences against candidate PAM/PFS rules, predicts Cas subtype from repeat/array features, and produces a bacteria-by-phage evidence matrix.

The primary output is a CRISPR targeting evidence matrix. Each bacterium-phage pair can include the targeting score, number of matching spacers, best spacer-protospacer identity, PAM/PFS support level, predicted Cas subtype, classifier confidence, and an evidence interpretation. Detailed tables for arrays, spacers, spacer hits, matrix rows, heatmap values, and metadata are saved for reproducibility.

## 3. Methods

### 3.1 FASTA Parsing and Input Diagnostics

Uploaded FASTA records are parsed into normalized sequence records with source file, record identifier, sequence length, GC content, accession extraction, and SHA-256 sequence hashes. SABR reports duplicate sequences and accession conflicts so users can identify redundant inputs or records with the same accession but different sequence content.

### 3.2 CRISPR Array Detection

The internal detector is a transparent baseline that searches for exact direct repeats of plausible CRISPR repeat length separated by plausible spacer lengths. This choice favors reproducibility and deployability over maximal sensitivity. It is not intended to replace established CRISPR callers. Optional MinCED-compatible detection is available when external dependencies are installed.

### 3.3 Spacer-Phage Matching

The internal matcher searches extracted spacer sequences against uploaded phage genomes on both strands. Optional BLASTN support is available when `blastn` and `makeblastdb` are present. BLAST-derived hits can include identity, alignment length, coverage, e-value, and bitscore. Exact matching remains a reproducible baseline.

### 3.4 PAM/PFS Evaluation

SABR evaluates candidate protospacer flanking sequences using subtype-aware PAM/PFS rules when a subtype prediction is available. PAM/PFS support is reported as evidence rather than as a definitive functional claim. Current production scoring caps PAM/PFS-unsupported rows at 39.0 when PAM/PFS was evaluated and no hit supports the expected rule.

### 3.5 Evidence Scoring

The preferred score name is `crispr_targeting_score`. Historical `hypothetical_resistance_score` aliases are retained only for backward compatibility with older artifacts. The score combines spacer-hit evidence, match quality, PAM/PFS support, and related diagnostics. SABR avoids wording that implies confirmed resistance from spacer matches alone.

### 3.6 Repeat-Derived Cas Subtype Prediction

The Cas subtype model uses features derivable from uploaded FASTA files, including repeat length, GC/AT composition, GC skew, spacer count, mean spacer length, terminal base composition, k-mer composition, terminal k-mer composition, reverse-complement self-identity, and simple hairpin-like repeat features. Taxonomy, organism name, accession, source database, and cas gene annotations are not used as runtime model features. These fields are retained for provenance and validation split design.

Evaluated models include nearest-repeat classification, logistic regression, linear SVM, random forest, ExtraTrees, hybrid ExtraTrees with nearest-neighbor features, and a broad type-to-subtype hierarchical ExtraTrees experiment.

### 3.7 Training Data Sources

The current best training table is `data/training/repeats_cas_types_augmented_vink_genbank_targeted.csv`. It combines computationally filtered CRISPRCasdb-derived Vink et al. rows with targeted GenBank/MinCED additions for selected underrepresented subtypes. The data are treated as a development set rather than a publication-grade gold standard.

Additional CRISPRCasdb release 34 files were organized under `data/training/external_sources/crisprcasdb_34/`. Direct-repeat FASTA exports were imported as unlabeled inventories. The PostgreSQL dump was parsed to build computational candidate labels by linking CRISPR loci to nearest same-sequence Cas clusters. This yielded 23,507 candidate rows, but audit showed substantial overlap with current data and 125 repeat hashes with conflicting subtype labels.

### 3.8 Validation Strategy

Random row splits are treated as smoke tests because they may overestimate performance through repeat, accession, or lineage leakage. Stronger validation uses genome, species, or genus holdout splits. The current main comparison uses genus holdout with subtypes below 20 rows removed.

## 4. Results

### 4.1 Best Current Predictor

The best current predictor is the flat ExtraTrees subtype classifier trained on `repeats_cas_types_augmented_vink_genbank_targeted.csv`. Under genus-holdout validation with subtypes below 20 rows removed, it achieved:

- Accuracy: 0.9152
- Rows used: 4,848
- Train rows: 3,633
- Test rows: 1,215

This model outperformed nearest-repeat classification, random forest, hybrid ExtraTrees, and hierarchical ExtraTrees in the current evaluation artifacts.

### 4.2 Model Comparison Summary

On the current best dataset:

- nearest repeat: 0.8922
- ExtraTrees: 0.9152
- hybrid ExtraTrees: 0.9119
- hierarchical ExtraTrees: 0.9012

The flat ExtraTrees model remains the runtime candidate. Hybrid and hierarchical variants have not yet justified replacing it.

### 4.3 CRISPRCasdb Augmentation Experiments

The CRISPRCasdb SQL candidate importer produced 23,507 computational candidate rows. After removing repeats already present in the current dataset, excluding repeat hashes with conflicting subtype labels, and capping additions at 500 rows per subtype, 3,535 additions were retained for a balanced augmentation experiment.

All-subtype CRISPRCasdb augmentation:

- Rows used after filtering: 8,378
- Accuracy: 0.9050
- III-A f1/recall: 0.68 / 0.66
- III-B f1/recall: 0.55 / 0.48
- III-D f1/recall: 0.38 / 0.32

Type III-only CRISPRCasdb augmentation:

- Added rows: 714
- Rows used after filtering: 5,578
- Accuracy: 0.9004
- III-A f1/recall: 0.76 / 0.86
- III-B f1/recall: 0.49 / 0.57
- III-D f1/recall: 0.33 / 0.26

These experiments show that CRISPRCasdb candidates are useful for improving III-A recall, but they do not broadly solve Type III subtype prediction and reduce overall accuracy.

### 4.4 Dimensionality Reduction and Dataset Presentation

We generated PCA and sampled t-SNE projections from the same repeat/array feature space used by the classifier. These plots are intended to show whether subtype labels form separable clusters, whether CRISPRCasdb additions occupy the same feature space as current training rows, and whether problematic classes such as Type III-B and Type III-D overlap with other subtypes.

Generated assets are stored under `docs/manuscript_assets/`. These visualizations should be interpreted descriptively; they are not validation metrics.

## 5. Discussion

SABR currently performs best as a transparent CRISPR targeting evidence mapper with repeat-derived Cas subtype support. The best model shows strong overall genus-holdout performance, but weak subtypes remain. Type III-A can be improved by targeted computational candidate additions, but Type III-B and Type III-D remain difficult. This suggests that repeat/array features alone may be insufficient for some Type III distinctions, or that current labels are noisy, overlapping, or too sparse.

The CRISPRCasdb experiments highlight a central point: adding more rows does not automatically improve prediction. Data additions improve model reliability only when labels are correct, rows are diverse, duplicate leakage is controlled, and subtype conflicts are removed. More entries can reduce performance if they introduce computational-label noise, class imbalance, repeated lineage-specific motifs, or ambiguous repeat-to-subtype mappings.

## 6. Probability and Confidence Interpretation

The ExtraTrees classifier can produce class probabilities. SABR can therefore report a predicted subtype and confidence-like value, such as `I-E` with probability 0.94. However, raw model probabilities should not be interpreted as exact biological probabilities unless calibrated on held-out validation data. A reported probability is best interpreted as model confidence given the training distribution and feature representation.

The probability of a correct Cas type or subtype call generally improves with more high-quality, diverse, correctly labeled training entries. It does not necessarily improve with raw dataset size. Rows that are duplicated, biased, mislabeled, or derived from ambiguous computational sources can lower real-world accuracy.

## 7. Limitations

SABR has several current limitations:

- The internal CRISPR detector is an exact-repeat baseline and may miss divergent arrays.
- Spacer matching is exact by default, although BLASTN can be used when available.
- PAM/PFS rules are incomplete and subtype-dependent.
- Repeat-derived Cas typing cannot confirm cas gene functionality.
- Current training labels include computational candidate rows and are not a final curated gold standard.
- Type III-B and Type III-D prediction remains weak.
- Current benchmark panels are small and not yet sufficient for calibrated resistance thresholds.
- SABR reports targeting evidence, not confirmed resistance.

## 8. Future Work

Priority next steps include:

1. Calibrate model probabilities using held-out validation.
2. Improve Type III labels using CCTyper or other locus-level cas-gene-supported sources.
3. Expand curated bacteria-phage benchmark pairs.
4. Add anti-CRISPR and PAM-failure controls.
5. Compare internal CRISPR detection against established tools.
6. Evaluate approximate spacer-protospacer matching.
7. Build manuscript-quality PCA/UMAP/t-SNE figures with clear legends and caveats.
8. Package SABR for easier reproducible deployment.

## 9. Current Conclusion

SABR provides a cautious, reproducible framework for mapping CRISPR spacer-targeting evidence and predicting likely CRISPR-Cas subtype from FASTA-derived repeat and array features. The current best predictor is a flat ExtraTrees model with 0.9152 genus-holdout accuracy. CRISPRCasdb expansion improves III-A but does not yet justify replacing the production model. The tool is scientifically useful when interpreted as an evidence mapper, not as a direct resistance caller.
