# Repeat/Cas-Type Training Data

This directory is reserved for curated CRISPR repeat to Cas type/subtype training tables.

The first target file is:

```text
data/training/repeats_cas_types.csv
```

Required columns are defined in `crispr_phage_predictor.ml.dataset`.

Initial dataset policy:

- Bootstrap from CRISPRCasTyper/CCTyper output because it exposes consensus repeats, spacer statistics, subtype predictions, probability, trusted status, and CRISPR-near-Cas tables.
- Treat nearby Cas-operon-supported subtype as the training label.
- Use repeat-based subtype predictions as a baseline/comparison, not as ground truth.
- Exclude orphan, ambiguous, hybrid, low-confidence, and incomplete systems from the first training set.
- Later, build a reproducible RefSeq/GenBank complete-genome pipeline and rerun annotation with fixed tool versions.

Validation policy:

- Split by genome, species, or genus rather than randomly by array.
- Keep taxonomy and organism metadata for auditing, but avoid using them as first-model features.
