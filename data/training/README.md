# Repeat/Cas-Type Training Data

This directory is reserved for curated CRISPR repeat-to-Cas type/subtype
training tables.

The first target file is:

```text
data/training/repeats_cas_types.csv
```

## Current Direction

The first training dataset should not use phage-host outcomes as labels and
should not use MinCED alone as the source of Cas type labels.

Use this split:

- MinCED-compatible detection through Diced, or the internal detector as a
  fallback, extracts CRISPR arrays, repeat consensus sequences, spacer counts,
  and spacer length statistics.
- Trusted literature, RefSeq/GenBank annotation, CRISPRCasdb/CRISPRCasFinder
  records, or other documented expert sources provide the Cas type/subtype
  labels.
- CCTyper/CRISPRCasTyper remains optional future scaffolding for Linux/WSL
  scale-up or label cross-checking, but it is not a blocking dependency.

Runtime model constraint:

- Training labels may come from annotation, literature, CRISPRCasdb, GenBank,
  or other trusted sources.
- Model input features should be derivable from a plain uploaded FASTA file.
- Do not use organism name, taxonomy, accession, source database, or Cas-gene
  annotation as first-pass runtime model features.
- Use those metadata fields for provenance, auditing, and genome/species/genus
  validation splits.

## Curated Manifest Workflow

Create a CSV or TSV manifest with one row per trusted labeled genome, contig, or
array/locus:

```csv
fasta_path,genome_id,organism,taxonomy,assembly_level,cas_type,cas_subtype,label_source,label_confidence,label_scope,contig_id,locus_start,locus_end,pam_rule,source,notes
data/curation/downloads/bacteria/example.fasta,GCF_example,Example bacterium,Bacteria,complete genome,Type I,I-F,curated_literature,curated,genome,,,,CC,curated_minced,replace with real source note
```

`label_scope` controls how labels are attached to detected arrays:

- `genome`: attach the label to all detected arrays in the FASTA. Use only when
  the genome/locus is known to contain one unambiguous CRISPR-Cas system.
- `contig`: attach the label to arrays on the listed `contig_id`.
- `array_coordinates`: attach the label only to arrays overlapping
  `contig_id`, `locus_start`, and `locus_end`.

Then run:

```bash
python -m crispr_phage_predictor.ml.collect_curated_minced_training data/training/curated_cas_type_manifest.tsv --output data/training/repeats_cas_types.csv
```

Use `--detector internal` for a pure-Python reproducible baseline, or
`--detector minced` to require Diced/MinCED.

## Optional CCTyper Import

The repository still includes CCTyper import helpers for future Linux/WSL use:

```bash
python -m crispr_phage_predictor.ml.collect_cctyper_training data/training/cctyper_manifest.csv --output data/training/repeats_cas_types.csv
```

Before using that route, check dependencies:

```bash
python -m crispr_phage_predictor.ml.check_cctyper_environment
```

## Scaled Candidate Dataset

For scaling experiments, the repository includes an importer for Vink et al.
2021 supplementary CRISPRCasdb-derived spacer/repeat data:

```bash
python -m crispr_phage_predictor.ml.import_vink2021_repeats data/training/external_sources/vink_2021_additional_file_2.csv --output data/training/repeats_cas_types_vink2021_candidate.csv --max-per-subtype 200
```

For a larger development table:

```bash
python -m crispr_phage_predictor.ml.import_vink2021_repeats data/training/external_sources/vink_2021_additional_file_2.csv --output data/training/repeats_cas_types_vink2021_candidate_1k.csv --max-per-subtype 1000
```

This creates a computational candidate table, not a manually curated gold
dataset. The importer keeps rows where the repeat subtype agrees with nearby
Cas-subtype metadata, collapses spacer-level records into accession + repeat +
subtype rows, and caps each subtype for balance.

Use this table for early scaling tests and feature/model development. For
publication claims, keep it separate from the manually curated seed set and
report it as a computationally filtered CRISPRCasdb-derived candidate dataset.

## CRISPRCasdb Raw Direct-Repeat Inventory

Raw CRISPRCasdb release 34 files are documented under:

```text
data/training/external_sources/crisprcasdb_34/
```

To create an unlabeled direct-repeat inventory from the local FASTA export:

```bash
python -m crispr_phage_predictor.ml.import_crisprcasdb_repeats data/training/external_sources/crisprcasdb_34/dr_34.zip --output data/training/crisprcasdb_34_direct_repeats_inventory.csv
```

This inventory is useful for coverage audits and future joins against the SQL
dump, but it is not a Cas type/subtype training table by itself.

To build a computational candidate repeat/Cas table from the extracted
CRISPRCasdb SQL dump:

```bash
python -m crispr_phage_predictor.ml.import_crisprcasdb_sql data/training/external_sources/crisprcasdb_34/home/pa.charbit/20220414_ccpp_recette_chromo_complete.sql --output data/training/repeats_cas_types_crisprcasdb_sql_candidate.csv
```

The importer links each CRISPR locus to the nearest same-sequence Cas cluster
within a configurable distance threshold and keeps only unambiguous
`CAS-Type...` subtype labels. Treat the output as computational candidates,
not curated gold labels.

Before using the candidate table for training, audit it against the current
best training table:

```bash
python -m crispr_phage_predictor.ml.audit_crisprcasdb_candidates data/training/repeats_cas_types_augmented_vink_genbank_targeted.csv data/training/repeats_cas_types_crisprcasdb_sql_candidate.csv --output-dir data/training/audits/crisprcasdb_sql_candidate
```

The audit reports overlap, repeat-to-subtype conflicts, subtype gains, and a
proposed balanced candidate subset size.

To build the filtered augmented table used for model experiments:

```bash
python -m crispr_phage_predictor.ml.build_crisprcasdb_augmented_dataset data/training/repeats_cas_types_augmented_vink_genbank_targeted.csv data/training/repeats_cas_types_crisprcasdb_sql_candidate.csv --output data/training/repeats_cas_types_augmented_crisprcasdb_sql_balanced.csv --candidate-output data/training/repeats_cas_types_crisprcasdb_sql_balanced_additions.csv
```

This keeps only novel, non-conflicting candidate repeats and caps additions at
500 rows per subtype by default.

Initial evaluation showed that the balanced CRISPRCasdb augmentation improved
`III-A` recall but reduced overall genus-holdout accuracy. Treat it as an
experiment source, not the production training table.

For a Type III-targeted experiment:

```bash
python -m crispr_phage_predictor.ml.build_crisprcasdb_augmented_dataset data/training/repeats_cas_types_augmented_vink_genbank_targeted.csv data/training/repeats_cas_types_crisprcasdb_sql_candidate.csv --include-subtypes III-A,III-B,III-C,III-D --output data/training/repeats_cas_types_augmented_crisprcasdb_typeiii_balanced.csv --candidate-output data/training/repeats_cas_types_crisprcasdb_typeiii_balanced_additions.csv
```

The Type III-targeted experiment improved `III-A` but did not improve `III-B`
or `III-D`; keep it as an experiment source for now.

## Dataset Policy

- Keep source provenance for every row.
- Exclude orphan, ambiguous, hybrid, low-confidence, incomplete, or multi-system
  records from the first training set unless the array-to-Cas association is
  unambiguous.
- Store organism and taxonomy for auditing and split design, not as first-model
  features.

## Validation Policy

- Split by genome, species, or genus rather than randomly by array.
- Compare a nearest-repeat similarity baseline against classical ML models.
- Treat subtype prediction confidence and calibration as core evaluation
  outputs, not optional polish.

## Phage-Host Interaction Features

`phage_host_interaction_features.tsv` is the first model-ready table for the
new sensitivity/resistance direction. It is generated from:

- `data/curation/phage_host_interactions.tsv`
- `data/curation/accession_linkage.tsv`
- `data/curation/accession_linkage_coverage.tsv`
- FASTA files under `data/curation/downloads/`

Regenerate it with:

```bash
python scripts/build_interaction_feature_table.py
```

The current table is hybrid-mode. Strict rows use exact or strain-alias host
genomes and exact phage genomes. Hybrid rows use exact phage genomes plus
downloaded reference-proxy host genomes. Keep `uses_reference_proxy_host` as a
feature or filtering column so proxy-based modeling is not confused with
exact-isolate validation.

Run the first smoke-test baseline with:

```bash
python scripts/train_phage_host_baseline.py
```

Outputs are written to `data/training/phage_host_baseline/`. Treat these as
engineering baselines only. Row-random splits are optimistic because related
rows share phages, host proxies, sources, and assay protocols. Grouped splits
by phage or source are more honest for generalization and currently show that
the model is not yet strong enough for biological claims.

Add first-pass CRISPR spacer-targeting features with:

```bash
python scripts/add_targeting_features.py
```

This creates `phage_host_interaction_features_with_targeting.tsv`. The current
targeting layer uses Diced/MinCED-compatible array detection when available,
exact spacer matching, and SABR's targeting-evidence score. Internal exact
repeat fallback is bounded to small FASTA records so large reference genomes do
not block feature generation.

The targeting layer also includes experimental fuzzy spacer/protospacer
features. These allow imperfect matches up to a bounded mismatch count and
summarize possible seed-edge conservation, distal mismatch burden, and a graded
CRISPR interference score. Evaluate targeting thresholds directly with:

```bash
python scripts/evaluate_targeting_thresholds.py
python scripts/evaluate_targeting_thresholds.py --restrict-tier tier1_exact_pair --output data/training/targeting_threshold_evaluation_tier1_exact_pair.tsv
```

Current threshold behavior is high precision but low recall: on tier1 exact-pair
rows, CRISPR targeting thresholds identify 10 resistant rows with no susceptible
false positives, but miss 21 resistant rows.

Add GenBank-derived phage annotation features with:

```bash
python scripts/download_phage_genbank_records.py
python scripts/add_phage_annotation_features.py
python scripts/train_phage_host_baseline.py --features data/training/phage_host_interaction_features_with_annotations.tsv --output-dir data/training/phage_host_baseline_with_annotations
```

This layer parses CDS product/gene/note/function text for keyword families such
as integrase, repressor, tail fiber, baseplate, depolymerase, holin/endolysin,
DNA methyltransferase, and anti-CRISPR-like terms. These are weak annotation
features and should be interpreted as coarse mechanistic hints, not curated
functional calls.

The current linked feature table uses explicit tiers from
`accession_linkage_coverage.tsv`:

- `tier1_exact_pair`: exact/strain-alias host genome and exact phage genome.
- `tier2_proxy_host_exact_phage`: reference-proxy host genome and exact phage
  genome.
- `tier4_host_only_or_proxy`: phenotype row with host exact/proxy linkage but
  no exact phage genome yet.

`dataset_tier` is metadata and is excluded from model features.

After the first NCBI Assembly host-promotion pass, strict exact-pair rows
increased from 13 to 74. The modelable table remains 213 rows because those are
the rows with exact phage genomes; the improvement is that many Mirzaei/Nilsson
ECOR rows now use exact host assemblies instead of the E. coli reference proxy.
The current combined targeting + annotation baseline is stored in
`phage_host_baseline_with_annotations/`.

To evaluate only strict exact host-phage genome rows:

```bash
python scripts/train_phage_host_baseline.py --features data/training/phage_host_interaction_features_with_annotations.tsv --output-dir data/training/phage_host_baseline_with_fuzzy_targeting_tier1_exact_pair --restrict-tier tier1_exact_pair
```
