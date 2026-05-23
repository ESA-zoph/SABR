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
