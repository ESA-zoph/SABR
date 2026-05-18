# CRISPR-Phage Predictor Technical Handout

## Active Direction

We are extending the MVP CRISPR-phage matching app with a repeat-based CRISPR-Cas type/subtype classifier.

The purpose of the classifier is to infer likely Cas type/subtype from CRISPR array sequence features when explicit genome annotation is absent. The predicted subtype will then choose candidate PAM/PFS rules for checking phage protospacer flanks, improving the resistance evidence score.

## Dataset Strategy

Preferred sequence:

1. Bootstrap from CRISPRCasTyper/CCTyper-derived outputs.
2. Use Cas-operon-supported subtype calls as labels.
3. Train first models on repeat/array-only features.
4. Build a reproducible RefSeq/GenBank complete-genome annotation pipeline after the local model workflow is working.

High-confidence training examples should have:

- trusted CRISPR array
- nearby Cas operon
- confident Cas subtype prediction
- non-ambiguous subtype label
- complete or high-quality source genome where possible

Excluded from the first training set:

- orphan arrays without confident support
- ambiguous or hybrid subtype calls
- incomplete systems
- low-confidence repeat predictions
- labels with no traceable source

## Dataset Schema

The first table is planned as `data/training/repeats_cas_types.csv`.

Columns:

```text
source
genome_id
organism
taxonomy
assembly_level
contig_id
array_start
array_end
repeat_sequence
repeat_length
spacer_count
mean_spacer_length
cas_type
cas_subtype
label_source
label_confidence
pam_rule
```

## CCTyper Importer

The local importer in `crispr_phage_predictor/ml/dataset.py` converts CCTyper `crisprs_near_cas.tab`-style tables into the local schema.

Expected CCTyper fields:

- `Contig`
- `Start`
- `End`
- `Consensus_repeat`
- `N_repeats`
- `Repeat_len`
- `Spacer_len_avg`
- `Trusted`
- `Subtype` or `Prediction`
- `Subtype_probability` when available

Current import rules:

- one output row per CRISPR array
- `spacer_count = N_repeats - 1`
- `cas_type` is inferred from subtype prefix, for example `I-E` becomes `Type I`
- `label_source = nearby_cas_operon`
- trusted rows with subtype probability at least 0.8 are marked high confidence
- ambiguous labels containing `/` or `,` are skipped in the first dataset
- invalid repeat sequences are skipped

## First Model Features

Start with:

- repeat sequence
- repeat length
- repeat GC percent
- repeat 3-mer counts
- repeat 4-mer counts
- spacer count
- mean spacer length
- spacer length variability when available

Do not use organism or taxonomy as first-model features because they can cause shortcut learning.

## Validation

Validation should hold out whole genomes, species, or genera. Random array-level splitting is likely to overestimate performance because closely related genomes can share repeat families.

## PAM/PFS Plan

After subtype prediction exists, add a curated PAM/PFS rule table keyed by subtype. PAM support should be reported as evidence categories such as compatible, weak, absent, or not evaluated rather than as definitive resistance.
