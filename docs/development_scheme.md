# SABR Development Scheme

SABR is being built as a cautious CRISPR-phage targeting evidence tool, not as a tool that directly claims confirmed resistance.

## Track 1: Core Analysis Pipeline

Goal: turn bacterial and phage FASTA inputs into a transparent bacteria-by-phage evidence matrix.

Steps:

1. FASTA parsing and validation.
2. Sequence hashing, accession extraction, and duplicate/conflict diagnostics.
3. CRISPR array detection.
4. Spacer extraction.
5. Spacer-phage matching with internal exact matching or BLASTN.
6. Protospacer flank extraction.
7. Cas subtype prediction from repeat/array features.
8. Candidate PAM/PFS rule selection from curated subtype rules.
9. PAM/PFS evaluation.
10. PAM-proximal seed mismatch summary.
11. Hypothetical CRISPR targeting evidence score.
12. Heatmap, evidence matrix, and detailed exports.

## Track 2: Data Expansion

Goal: increase the amount and quality of data behind SABR.

Steps:

1. Expand repeat/Cas subtype training rows toward 20,000-50,000 rows.
2. Prioritize underrepresented Type III, Type V, and Type VI systems.
3. Use locus-level labels when genomes contain multiple CRISPR-Cas systems.
4. Expand curated bacteria-phage benchmark pairs toward 100-300 pairs.
5. Include positives, negatives, anti-CRISPR cases, host-range controls, and PAM-failure controls.
6. Track accession, source, confidence, caveats, and local file status for every curated row.

## Track 3: Model Development

Goal: make Cas subtype prediction useful without shortcut learning.

Steps:

1. Keep nearest-repeat prediction as an interpretable baseline.
2. Compare tabular ML models on FASTA-only repeat/array features.
3. Validate with genome, species, and genus holdout splits.
4. Calibrate prediction confidence.
5. Use subtype confidence to decide whether PAM/PFS support should be applied.
6. Add neural models only after enough high-confidence, diverse labels exist.

## Track 4: Benchmarking and Calibration

Goal: evaluate whether SABR scores match known biology without overclaiming.

Steps:

1. Build a formal bacteria-phage benchmark label table.
2. Separate CRISPR-mediated resistance from receptor/host-range resistance.
3. Flag spacer-only evidence, PAM-failure cases, and anti-CRISPR cases.
4. Run SABR on the benchmark panel.
5. Compare scores against curated labels.
6. Calibrate score thresholds only after the benchmark panel is large enough.

## Track 5: Productization

Goal: make SABR easy to run, inspect, and reproduce.

Steps:

1. Streamlit GUI for interactive use.
2. CLI pipeline for batch and publication workflows.
3. Saved timestamped run artifacts.
4. Documentation for inputs, outputs, caveats, and methods.
5. Docker or environment packaging.
6. Example datasets and expected-output walkthroughs.

## Current Near-Term Order

1. Structured benchmark label schema.
2. Benchmark validation tooling.
3. Larger curated benchmark panel.
4. Larger repeat/Cas subtype data acquisition.
5. Score calibration against benchmark labels.
