# Repeat/Cas-Type Training Data

This directory is reserved for curated CRISPR repeat to Cas type/subtype training tables.

The first target file is:

```text
data/training/repeats_cas_types.csv
```

To collect completed CCTyper output folders into this file, create a manifest CSV:

```csv
cctyper_output_dir,genome_id,organism,taxonomy,assembly_level
data/training/cctyper_outputs/GCF_000001,GCF_000001,Example bacterium,Bacteria,complete genome
```

Then run:

```bash
python -m crispr_phage_predictor.ml.collect_cctyper_training data/training/cctyper_manifest.csv --output data/training/repeats_cas_types.csv
```

Before running CCTyper locally, check dependencies:

```bash
python -m crispr_phage_predictor.ml.check_cctyper_environment
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
