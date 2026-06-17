# SABR: A Transparent Workflow for CRISPR Spacer-Targeting Evidence and Repeat-Derived Cas Subtype Prediction

## Working Title

SABR: Spacer Alignment-Based Recognition for cautious CRISPR-phage targeting evidence and FASTA-derived CRISPR-Cas subtype prediction

## Abstract

CRISPR spacer matches to bacteriophage genomes can support hypotheses of prior exposure or potential immune targeting, but they cannot alone demonstrate bacterial resistance. Functional interference depends on additional evidence, including the CRISPR-Cas system subtype, appropriate protospacer adjacent motif or protospacer flanking sequence (PAM/PFS), seed-region conservation, intact cas machinery, expression, phage escape, and anti-CRISPR activity. We developed SABR, a transparent bioinformatics workflow that accepts bacterial and phage FASTA files, detects candidate CRISPR arrays, extracts spacers, searches for phage protospacer matches, evaluates subtype-aware PAM/PFS support, predicts likely Cas subtype from repeat and array features, and exports a bacteria-by-phage targeting evidence matrix. SABR is deliberately an evidence mapper rather than a resistance classifier.

For subtype inference, we evaluated interpretable FASTA-derived models and selected a flat ExtraTrees classifier as the current SABR subtype model. The selected model was trained on 23,478 computationally annotated repeat-array rows after minimum-class filtering and achieved 0.9455 accuracy under genome-held-out internal evaluation. It was better calibrated on this split than an earlier development-table model (expected calibration error 0.0362 versus 0.2254). A first independent CCTyper-supported pilot comprised 25 array rows from 11 genomes and was classified correctly in all represented subtypes (25/25), although this pilot was dominated by Type I arrays. A subsequent enriched strict CCTyper screen retained 20 strict-confirmed genomes and 25 expected-subtype array rows across I-B, I-C, III-A, and III-D; SABR correctly classified 23/25 rows (accuracy 0.9200), with the two errors being I-C arrays predicted as III-B. The strict screen still did not yield an expected-subtype III-B array row, so Type III-B remains insufficiently tested. SABR provides a reproducible framework for layered CRISPR targeting evidence while preserving the distinction between computational evidence and experimentally confirmed phage resistance.

## Keywords

CRISPR, bacteriophage, spacer, protospacer, PAM, PFS, CRISPR-Cas subtype, ExtraTrees, evidence mapping, CCTyper

## 1. Introduction

CRISPR-Cas adaptive immunity stores fragments of previous mobile genetic element encounters as spacers within bacterial or archaeal CRISPR arrays. When a spacer matches a protospacer in a phage genome, the match is biologically meaningful evidence: it may represent a historical encounter, a potential interference target, or a feature associated with reduced infection. Nevertheless, a sequence match is not equivalent to experimentally demonstrated resistance. An apparently matching phage can remain infective if the required PAM/PFS is absent, if critical seed positions are altered, if the bacterial Cas machinery is inactive, if the array is not expressed, or if the phage encodes evasion mechanisms.

SABR was developed to make spacer-based analysis explicit, inspectable, and cautious. The workflow combines FASTA ingestion, candidate array detection, spacer-to-phage matching, subtype inference, PAM/PFS evaluation, and evidence scoring within one reproducible interface. Rather than returning a binary resistant/susceptible call, SABR returns evidence components and an interpretable targeting score that can be audited against the underlying sequences and model outputs.

A central part of this work is repeat-derived Cas subtype prediction. Cas subtype is relevant because interference biology and PAM/PFS rules differ across CRISPR-Cas systems. SABR must operate on uploaded sequence data, so the runtime model is designed to rely on repeat and array features available from FASTA input rather than on taxonomy, organism identity, database source, or pre-existing cas annotations. Annotation-bearing resources and independent software are used to obtain development or validation labels, not as competing predictive models in the SABR workflow.

The objectives of this study were to:

1. Implement an end-to-end, auditable workflow for CRISPR spacer-targeting evidence against uploaded phage genomes.
2. Develop a FASTA-derived Cas subtype classifier suitable for integration in that workflow.
3. Evaluate candidate subtype models using leakage-aware grouped validation and confidence diagnostics.
4. Assess the selected model on independent CCTyper-supported pilot and enriched strict-confirmation datasets.
5. Identify remaining error regions and define the validation needed before biological resistance claims are considered.

## 2. Software and Analytical Workflow

### 2.1 SABR Inputs and Outputs

SABR accepts one or more bacterial genome FASTA files and one or more phage genome FASTA files through a Streamlit interface. Input diagnostics include parsed record count, record identifier, sequence length, GC content, accession extraction when available, and SHA-256 sequence hashes used to detect duplicate or conflicting uploads.

For each bacterial-phage pair, SABR can report candidate arrays, extracted spacers, spacer-protospacer hits, orientation, matching statistics, predicted Cas subtype and model confidence, PAM/PFS support, seed-related diagnostics, and a composite `crispr_targeting_score`. Run metadata and detailed tables are written to local output folders so that a graphical result can be traced to source sequences and intermediate calculations.

### 2.2 Candidate CRISPR Array Detection

The current internal detector searches bacterial sequences for exact direct repeats within plausible CRISPR-repeat length ranges and with intervening spacer lengths compatible with candidate arrays. This intentionally transparent detector provides a reproducible baseline and supports immediate deployment without external binaries. It is not presented as a replacement for specialized CRISPR discovery tools. SABR can also use a MinCED-compatible external detector when available, enabling future sensitivity comparisons.

### 2.3 Spacer-to-Phage Matching

Extracted spacers are compared with uploaded phage genomes on both strands. The internal matching mode records exact spacer-protospacer matches as a conservative evidence layer. When external dependencies are installed, BLASTN matching can report approximate alignments with identity, aligned length, spacer coverage, E-value, and bitscore. This layered design preserves a deterministic baseline while supporting more permissive analysis where explicitly selected.

### 2.4 PAM/PFS and Seed Evidence

When a likely Cas subtype is available, SABR evaluates flanking nucleotides against the subtype-associated PAM or PFS rules presently encoded in the system. PAM/PFS observations are evidence components, not definitive functional proof. In the current evidence score, a row evaluated for PAM/PFS but without a supporting hit is capped at 39.0, preventing an unsupported exact sequence match from being presented as high-confidence targeting evidence. Seed-region summaries are retained as further diagnostics for approximate alignments and possible escape variants.

### 2.5 Evidence Score Terminology

The current result field is named `crispr_targeting_score`. Earlier internal artifacts may contain the historical name `hypothetical_resistance_score`; that alias is retained for compatibility only. The manuscript uses targeting-evidence terminology throughout because SABR has not yet been validated as a phenotype-level resistance predictor.

### 2.6 Mathematical Definition of the Targeting-Evidence Score

The implemented targeting score is an explicit bounded evidence function, not a fitted probability of resistance. For a bacterial-phage pair with one or more detected spacer hits, let `n_s` denote the number of unique matching spacers, `I_max` the best alignment identity expressed on a 0-1 scale, `C_max` the best alignment coverage, and `q_max` the highest predicted subtype confidence associated with evaluated hits.

Equation 1. `S_spacer = min(35, 18 * log2(n_s + 1))`

Equation 2. `S_identity = min(25, 25 * I_max * C_max)`

The PAM/PFS component is:

Equation 3. `S_PAM = 20` if at least one evaluated hit is PAM/PFS-supported; `S_PAM = -25` if flanks were evaluated but none supported the rule; otherwise `S_PAM = 0`.

The seed component is applied only when at least one hit has supported PAM/PFS. If `m_seed` is the minimum PAM-proximal seed mismatch count among evaluated hits:

Equation 4. `S_seed = 10` when `m_seed = 0`; `S_seed = 5` when `1 <= m_seed <= 2`; `S_seed = -5` when `m_seed > 2`; otherwise `S_seed = 0`.

Subtype confidence contributes:

Equation 5. `S_subtype = min(10, 10 * q_max)` when confidence is available, and `0` otherwise.

The final targeting-evidence score is:

Equation 6. `S_target = min(C_PAM, max(0, S_spacer + S_identity + S_PAM + S_seed + S_subtype))`

where `C_PAM = 39` if PAM/PFS was evaluated but unsupported by every hit, and `C_PAM = 100` otherwise. The score is interpreted as weak evidence when `0 < S_target < 50`, moderate evidence when `50 <= S_target < 75`, and strong candidate targeting evidence when `S_target >= 75`.

## 3. Cas Subtype Model Development

### 3.1 Model Scope

The subtype classifier predicts likely CRISPR-Cas subtype from an observed repeat and associated array summary. It is part of SABR and is referred to hereafter as the **selected SABR ExtraTrees subtype model**. Database names identify label provenance only; they do not identify another model or a competing tool.

### 3.2 FASTA-Derived Feature Representation

The evaluated feature space includes:

- repeat length, GC percentage, AT percentage, and GC skew;
- whole-repeat k-mer frequencies and terminal k-mer features;
- terminal base composition;
- spacer count and spacer-length summaries;
- spacer-to-repeat length ratio and estimated array geometry;
- reverse-complement self-identity, longest inverted-stem proxy, and hairpin-like repeat features.

The runtime feature matrix excludes organism name, taxonomy, accession, source database, and cas-gene annotation. Those fields may be retained as provenance variables and used to define held-out evaluation groups, but they do not enter prediction.

### 3.3 Development Label Sources and Provenance

Several development tables were assembled during model construction. Early tables combined candidate repeat/subtype labels derived from published CRISPR annotation resources with targeted GenBank and MinCED-supported rows. A larger candidate-label table was subsequently built from CRISPRCasdb release 34 by linking array repeats to nearby same-sequence Cas clusters represented in the source export. After filtering subtypes with fewer than 20 rows for evaluation, this table supplied 23,478 rows for the selected SABR model experiment.

CRISPRCasdb is therefore a source of computationally derived training annotations, not a model against which SABR is compared. These annotations are suitable for model development and internal grouped evaluation but are not treated as independently curated biological truth. An audit identified substantial overlap with earlier development data and 125 repeat hashes with conflicting subtype assignments, which is why independent validation is required.

### 3.4 Evaluated Models

Model development considered:

- nearest-repeat classification as an interpretable sequence-similarity baseline;
- logistic regression and linear support-vector classification;
- random forest;
- flat ExtraTrees;
- hybrid ExtraTrees incorporating nearest-neighbor confidence features;
- an experimental hierarchical ExtraTrees structure predicting broad type before subtype.

The selected architecture is flat ExtraTrees because it achieved the strongest documented grouped-evaluation result while remaining compatible with interpretable FASTA-derived tabular features.

### 3.5 Validation Design and Metrics

Random row splitting was used only for early smoke tests because rows from related genomes or shared repeats can inflate apparent performance. More stringent experiments held out genome/accession groups or genera so that closely related development records did not freely cross the train-test boundary. Reported metrics include overall accuracy, subtype precision, recall and F1 for difficult Type III classes, confidence calibration, and error-pair frequency. For future independent panels, macro-F1 and full confusion matrices will be added alongside accuracy.

For subtype `k`, performance was evaluated using:

Equation 7. `Precision_k = TP_k / (TP_k + FP_k)`

Equation 8. `Recall_k = TP_k / (TP_k + FN_k)`

Equation 9. `F1_k = 2 * Precision_k * Recall_k / (Precision_k + Recall_k)`

For confidence calibration, expected calibration error (ECE) was computed over confidence bins `b = 1,...,B`:

Equation 10. `ECE = sum_b (n_b / N) * abs(accuracy_b - confidence_b)`

Receiver operating characteristic (ROC) analyses used one-vs-rest subtype discrimination. During architecture selection, probability-capable models were compared on the identical genus-held-out evaluation rows; macro AUC was averaged only over subtypes represented by both positive and negative test examples.

### 3.6 Dimensionality Reduction and Visual Error Analysis

Principal component analysis (PCA) and t-distributed stochastic neighbor embedding (t-SNE) were applied to standardized repeat/array feature vectors for descriptive visualization. These projections were not used to train the subtype classifier and are not independent metrics of accuracy. They were used to visualize broad Cas-type grouping, overlap in Type III and adjacent Type I subtypes, and the spatial distribution of wrong calls.

## 4. Results

### 4.1 End-to-End Evidence Workflow

The implemented SABR application parses multi-FASTA inputs, detects candidate arrays, extracts repeats and spacers, searches spacers against phage sequences, performs subtype-aware PAM/PFS checks, exports run-level metadata and detailed result tables, and renders the bacteria-by-phage evidence matrix. Internal exact methods are available as reproducible baselines; optional external detection and BLASTN matching are supported when installed.

A curated benchmark scoring exercise demonstrated the intended interpretation of the evidence layer: high evidence requires supporting components, while PAM/PFS-unsupported matches are prevented from reaching the high-score range by the current 39.0 cap. This benchmark assesses score behavior and reporting logic; it is not a phenotype-level sensitivity or resistance trial.

The full evidence-matrix view highlights the reporting purpose of SABR: in the current demonstration run, non-zero scores are sparse and biologically interpretable as candidate targeting signals, including strong evidence for selected *Pseudomonas aeruginosa* PA14/JBD phage combinations and an archaeal SIRV2-associated signal. Zero-score cells represent no detected evidence under the selected pipeline, not proof of susceptibility.

### 4.2 Architecture Selection on the Earlier Development Table

The earlier targeted development table included 4,848 eligible rows under genus-held-out evaluation after applying the minimum class-count threshold. Flat ExtraTrees achieved the best result in that evaluation, exceeding the similarity baseline and the tested hybrid and hierarchical variants.

| Model on development table | Split design | Rows evaluated | Accuracy |
| --- | --- | ---: | ---: |
| Nearest-repeat baseline | Genus holdout | 4,848 | 0.8922 |
| Flat ExtraTrees | Genus holdout | 4,848 | 0.9152 |
| Hybrid ExtraTrees | Genus holdout | 4,848 | 0.9119 |
| Hierarchical ExtraTrees | Genus holdout | 4,848 | 0.9012 |

These experiments supported the choice of a flat ExtraTrees architecture. The hierarchical formulation did not improve overall prediction and reduced performance in difficult Type III classes.

ROC analysis was performed for the probability-capable architectures evaluated on the same genus-held-out development split. ExtraTrees produced the highest accuracy and the highest macro one-vs-rest AUC, marginally above random forest and clearly above logistic regression.

| Architecture on common genus holdout | Accuracy | Micro AUC | Macro AUC |
| --- | ---: | ---: | ---: |
| Logistic regression | 0.8337 | 0.9515 | 0.9144 |
| Random forest | 0.9128 | 0.9916 | 0.9719 |
| ExtraTrees | 0.9152 | 0.9927 | 0.9725 |

### 4.3 Training-Annotation Expansion Experiments

We next tested whether additional computational candidate annotations improved the same ExtraTrees architecture. Adding 3,535 non-conflicting, subtype-capped candidate rows to the earlier development table reduced genus-held-out accuracy from 0.9152 to 0.9050. Restricting additions to 714 Type III candidate rows raised III-A recall but yielded overall accuracy of 0.9004 and did not resolve III-B or III-D performance.

| ExtraTrees training-table experiment | Split design | Rows evaluated | Accuracy | III-A F1 | III-B F1 | III-D F1 |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Earlier development table | Genus holdout | 4,848 | 0.9152 | 0.58 | 0.59 | 0.40 |
| Broad candidate-label additions | Genus holdout | 8,378 | 0.9050 | 0.68 | 0.55 | 0.38 |
| Type III-focused additions | Genus holdout | 5,578 | 0.9004 | 0.76 | 0.49 | 0.33 |

These results show that increasing annotation volume alone is insufficient: computational label quality, conflict filtering, subtype coverage, and lineage diversity materially affect generalization.

### 4.4 Selected SABR ExtraTrees Subtype Model

The current selected SABR model uses the flat ExtraTrees architecture trained on the larger computational annotation-derived table. Under genome-held-out evaluation, after removing subtypes represented by fewer than 20 rows, the model used 23,478 rows and achieved 0.9455 accuracy. The nearest-repeat baseline on the same evaluation achieved 0.9389.

| Model | Training rows available for evaluation | Split design | Accuracy | III-A F1 | III-B F1 | III-D F1 |
| --- | ---: | --- | ---: | ---: | ---: | ---: |
| Nearest-repeat baseline | 23,478 | Genome holdout | 0.9389 | 0.86 | 0.52 | 0.54 |
| Selected SABR ExtraTrees subtype model | 23,478 | Genome holdout | 0.9455 | 0.89 | 0.58 | 0.54 |

This internal result supports use of the selected model in SABR development, but it is not an independent external accuracy estimate because the held-out rows derive from the same computational annotation source used for training-table construction.

The selected model's one-vs-rest ROC analysis was additionally inspected for the Type III and adjacent Type I classes implicated by the error analysis. The corresponding normalized confusion heatmap makes clear that strong overall accuracy coexists with subtype-specific difficulty: the principal remaining reductions in recall occur in III-B, III-C, and III-D rather than across all subtypes uniformly.

### 4.5 Compatibility Transfer Checks

As a development consistency check, the selected SABR architecture trained on the larger annotation-derived table was evaluated on compatible rows in the earlier SABR development table. It achieved 0.9811 accuracy across 4,806 rows with represented labels. In the reciprocal experiment, the earlier-table-trained ExtraTrees model achieved 0.9559 across 22,981 compatible rows from the larger table.

These transfer checks determine whether two development-label constructions lead to compatible repeat-to-subtype mappings. They must not be interpreted as comparison of SABR with a separate database model, and they are not substitutes for an independent validation dataset.

### 4.6 Confidence Calibration

Confidence diagnostics differed substantially between model-development stages. On the genus-held-out earlier development table, the ExtraTrees model achieved accuracy of 0.9152 with mean maximum predicted confidence of 0.6898 and expected calibration error (ECE) of 0.2254, indicating under-confidence. The selected SABR model, evaluated by genome holdout within its training-annotation source, achieved accuracy of 0.9455, mean confidence of 0.9227, and ECE of 0.0362.

| Model evaluation | Accuracy | Mean predicted confidence | ECE |
| --- | ---: | ---: | ---: |
| Earlier development-table ExtraTrees, genus holdout | 0.9152 | 0.6898 | 0.2254 |
| Selected SABR ExtraTrees model, genome holdout | 0.9455 | 0.9227 | 0.0362 |

The better internal calibration is encouraging, but probability values should remain confidence indicators rather than biological probabilities until calibration is confirmed on broader independent data.

### 4.7 Model Interpretability and Error Analysis

Built-in ExtraTrees importance, held-out permutation analyses, feature-category aggregation, and Type III-focused error summaries were generated for the selected SABR model. The largest built-in contributions involved repeat length, spacer/repeat length ratio, spacer-length statistics, terminal repeat k-mers, repeat composition, and hairpin-like proxies. When grouped biologically, whole-repeat k-mers and terminal k-mers contributed the largest total built-in importance, followed by terminal composition, array statistics, repeat composition, and repeat-structure features.

Permutation importance values were comparatively small, as expected when numerous sequence-derived features are correlated. Consequently, these analyses support interpretation at the level of feature families rather than causal claims about individual motifs.

Residual wrong calls were concentrated rather than uniformly distributed. The most frequent internal confusion pairs included:

| True subtype | Predicted subtype | Wrong calls |
| --- | --- | ---: |
| III-B | I-B | 37 |
| III-A | I-B | 25 |
| III-D | III-A | 18 |
| III-B | I-C | 17 |
| III-B | I-A | 16 |
| I-B | III-B | 14 |
| III-D | III-B | 14 |

This pattern identifies Type III-B, Type III-D, Type III-A, and adjacent Type I subtypes as the priority region for independent enrichment and future model improvement.

### 4.8 Dimensionality-Reduction View of the Error Region

t-SNE visualization of held-out repeat/array feature vectors demonstrates broad clustering by major Cas type while also revealing overlap between Type III and neighboring Type I subtypes. In the focused boundary-region view, wrong calls are concentrated within overlapping feature-space neighborhoods rather than appearing uniformly across the projection. This visual pattern is consistent with the confusion analysis and supports targeted external sampling of III-B, III-D, III-A, I-B, I-C, and I-A rather than undirected expansion.

### 4.9 Independent CCTyper-Supported Validation

An Ubuntu/VirtualBox installation of CRISPRCasTyper/CCTyper version 1.8.0 was used to produce an initial independent subtype-labelled pilot. Nineteen bacterial FASTA inputs were processed; 11 genomes generated importer-compatible `crisprs_near_cas.tab` outputs. Import into SABR produced 25 CCTyper-supported array rows.

| CCTyper subtype | Imported array rows | Correct SABR predictions |
| --- | ---: | ---: |
| I-E | 11 | 11 |
| I-A | 8 | 8 |
| I-F | 3 | 3 |
| II-A | 1 | 1 |
| II-C | 1 | 1 |
| III-A | 1 | 1 |
| Total | 25 | 25 |

The selected SABR model correctly classified all 25 imported arrays in this pilot (accuracy 1.0000). This is a useful first independent technical validation, but its scope is limited: it is small, concentrated in common Type I subtypes, contains one Type III-A row, and contains no Type III-B or III-D rows. It therefore did not test the internal model's principal error region, motivating the enriched strict CCTyper screen below.

We therefore assembled a second candidate pool enriched for Type III and adjacent Type I subtypes and screened it with CCTyper. Only genomes with a strict subtype-consistent CCTyper Cas-operon call were retained. This conservative filter produced 20 strict-confirmed genomes: 6 I-B, 5 I-C, 5 III-A, 1 III-B, and 3 III-D. No I-A genome passed strict confirmation, although several had hybrid or partial subtype-compatible predictions.

Because CCTyper operon confirmation and array-level subtype linkage are not identical, the imported validation table was filtered to expected-subtype array rows. This yielded 25 evaluable rows: 8 I-B, 6 I-C, 6 III-A, and 5 III-D. The selected SABR model correctly classified 23 of 25 rows (accuracy 0.9200).

| Strict CCTyper validation subtype | Strict-confirmed genomes | Expected-subtype array rows | Correct SABR predictions |
| --- | ---: | ---: | ---: |
| I-B | 6 | 8 | 8 |
| I-C | 5 | 6 | 4 |
| III-A | 5 | 6 | 6 |
| III-B | 1 | 0 | 0 |
| III-D | 3 | 5 | 5 |
| Total | 20 | 25 | 23 |

The two errors were I-C rows predicted as III-B. The strict-confirmed III-B genome did not contribute an expected-subtype III-B array row; its imported nearby array was labelled I-G by CCTyper. Two strict-confirmed III-D genomes similarly had nearby rows labelled III-B rather than III-D, leaving the Nostoc III-D genome as the source of the five expected III-D rows. These results strengthen external support for I-B, III-A, and III-D repeat-derived predictions but show that III-B remains unresolved and that strict genome-level confirmation should not be overinterpreted as balanced array-level validation.

## 5. Discussion

SABR integrates several evidence layers required for careful interpretation of spacer-phage relationships. The workflow is not designed to convert an exact spacer match into an unsupported resistance statement. Instead, it exposes candidate array structure, protospacer evidence, likely Cas subtype, PAM/PFS compatibility, confidence, and scoring outputs in a form that can be reviewed and exported.

The subtype-model results establish a defensible current model choice. Flat ExtraTrees consistently provided the strongest architecture among the evaluated FASTA-derived approaches, and the selected SABR model achieved high genome-held-out internal accuracy with substantially improved internal calibration. Importantly, this is one SABR model developed using annotation-derived training material; it is not a competition between SABR and CRISPRCasdb. CRISPRCasdb supplied computational candidate labels, while CCTyper is being used to build independent validation material.

The CCTyper validation results demonstrate that the trained classifier can agree with independent locus-supported calls across multiple subtype groups. However, neither the perfect first pilot nor the enriched strict screen should be overinterpreted. The first pilot did not challenge the known Type III error region, while the enriched strict screen was reduced by conservative CCTyper filtering and still did not yield an evaluable III-B array-level example. Additional III-B and III-D validation remains necessary to determine whether strong internal performance extends to the subtypes most likely to be confused.

The results also caution against assuming that more computational labels automatically yield a better classifier. Broad additions and targeted Type III additions to the earlier development table altered subtype metrics without improving overall held-out accuracy. Future development should prioritize independently supported, subtype-balanced examples and maintain separation between training expansion data and external evaluation data.

## 6. Limitations

The current study has several limitations:

1. The internal CRISPR detector is an exact-repeat baseline and may miss divergent or irregular arrays.
2. Exact spacer matching is conservative; approximate matching and alignment interpretation require further benchmarking.
3. PAM/PFS rules are incomplete and depend on accurate subtype assignment.
4. Repeat/array-derived subtype inference cannot confirm cas locus integrity or expression.
5. The selected model's large training table uses computational candidate annotations rather than a fully manually curated gold standard.
6. Internal errors remain concentrated in Type III and neighboring Type I classes.
7. The independent CCTyper validation remains small after strict filtering, and Type III-B still lacks an evaluable expected-subtype array row.
8. The bacteria-phage evidence benchmark is not yet large enough for phenotype-level resistance calibration.
9. Anti-CRISPR activity and other phage escape mechanisms are not yet incorporated as complete evidence layers.

## 7. Ongoing and Planned Work

The immediate ongoing experiment is continued enrichment of independent CCTyper validation material, especially array-level III-B and additional III-D examples. Future evaluations will use the selected SABR model without retraining and will report accuracy, macro-F1, subtype recall, confidence calibration, and confusion structure.

Subsequent work will:

1. Expand independently supported subtype labels while preserving a locked external evaluation panel.
2. Recover or assemble phenotype-labelled bacteria-phage pairs for targeting-evidence validation separate from subtype validation.
3. Benchmark internal CRISPR detection against established callers.
4. Benchmark exact and approximate protospacer matching settings.
5. Extend PAM/PFS support and anti-CRISPR/evasion evidence.
6. Prepare a reproducible release with frozen model metadata, data provenance tables, and manuscript-quality figures.

## 8. Conclusion

SABR is a transparent FASTA-driven framework for CRISPR spacer-targeting evidence and repeat-derived Cas subtype prediction. Its selected ExtraTrees subtype model achieved 0.9455 genome-held-out internal accuracy, correctly classified all 25 arrays in a first CCTyper-supported pilot, and correctly classified 23 of 25 expected-subtype rows in an enriched strict CCTyper validation screen. Clear limitations remain for Type III-B-focused independent validation and phenotype-level interpretation. The present evidence supports SABR as an auditable targeting-evidence tool and subtype-prediction framework, not yet as a direct predictor of bacteriophage resistance.

## Figure Plan for the First Draft

1. Targeting-evidence heatmap from the demonstration run, showing sparse bacteria-phage signal patterns.
2. Benchmark targeting-evidence scores, including the PAM/PFS-unsupported cap.
3. Common-split ROC comparison of probability-capable model architectures tested before model selection.
4. Model-development accuracy and Type III subtype performance summaries.
5. Selected SABR model ROC curves and normalized confusion heatmap.
6. Held-out t-SNE visualization by broad Cas type and focused Type III/Type I boundary subtypes.
7. Confidence calibration for the earlier development-stage and selected SABR models.
8. Feature-importance categories and leading FASTA-derived features for the selected model.
9. Error-pair summary highlighting Type III and adjacent Type I confusions.
10. Independent CCTyper validation subtype distributions and correct classifications, including the pilot and enriched strict-confirmation screen.

## Reporting Note

All statements about the selected subtype model distinguish internal annotation-derived evaluation, independent CCTyper-supported validation, unresolved Type III-B validation, and untested phenotype-level resistance claims.
