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

The collector command combines multiple completed CCTyper output directories:

```powershell
python -m crispr_phage_predictor.ml.collect_cctyper_training data/training/cctyper_manifest.csv --output data/training/repeats_cas_types.csv
```

The manifest columns are:

```text
cctyper_output_dir
genome_id
organism
taxonomy
assembly_level
```

Before running CCTyper locally, check the environment:

```powershell
python -m crispr_phage_predictor.ml.check_cctyper_environment
```

Required runtime pieces:

- `cctyper`
- `prodigal`
- `hmmsearch` from HMMER
- `minced`
- `blastn`
- CCTyper database path via `CCTYPER_DB` or `--db`

## Local CCTyper Setup Status

Current Windows-native status:

- Miniconda was installed successfully.
- The existing Python CCTyper launcher was repaired enough to show help by pinning compatible `setuptools<81` and `drawSvg<2`.
- A pilot CCTyper run reached the runtime stage but failed because the CCTyper database was missing.
- Native Windows conda could not solve the full CCTyper environment because the bioconda dependency stack is Linux-oriented.

Current missing pieces for local CCTyper annotation:

- Prodigal
- HMMER / `hmmsearch`
- MinCED
- CCTyper database

Preferred reproducible path for annotation is Linux/WSL or another Linux environment:

```bash
conda env create -f envs/cctyper-linux.yml
conda activate crispr-cctyper
cctyper genome.fa output_dir --no_plot --simplelog -t 4
```

CCTyper's official quick start uses:

```bash
conda create -n cctyper -c conda-forge -c bioconda -c russel88 cctyper
```

This route is preferred because CCTyper's maintainers state conda installs the software, dependencies, and database together.

Current WSL note:

- Running `wsl --install` from the Codex shell did not complete.
- Windows reported that WSL is not installed and should be installed with `wsl.exe --install`.
- This likely needs an elevated/admin PowerShell session and may require a restart.
- User ran `wsl.exe --install` in Administrator PowerShell.
- WSL 2.7.3 and `VirtualMachinePlatform` installed successfully.
- Windows reported that the change will not take effect until reboot.
- Next step is to reboot Windows, finish Ubuntu first-run setup, then run the Linux conda/CCTyper setup.

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

Implemented feature extraction currently creates:

- `repeat_length`
- `repeat_gc_percent`
- `spacer_count`
- `mean_spacer_length`
- normalized DNA 2-mer, 3-mer, and 4-mer frequencies by default

The feature table keeps `genome_id`, `contig_id`, `cas_type`, and `cas_subtype` for audit/labels, but `feature_columns()` excludes them from model inputs.

## Validation

Validation should hold out whole genomes, species, or genera. Random array-level splitting is likely to overestimate performance because closely related genomes can share repeat families.

## Baseline Classifier

The first implemented classifier is `RepeatCasSubtypeClassifier` in `crispr_phage_predictor/ml/classifier.py`.

Current behavior:

- builds numeric repeat/array features from the local training table
- trains a balanced random forest classifier
- predicts `cas_subtype`
- returns a confidence value and full subtype probability map
- requires at least two subtype labels for training

This is a baseline model, not a final scientific classifier. It should be compared against nearest-repeat similarity and evaluated with genome/species/genus holdout splits before interpretation.

## Training/Evaluation Command

The first CLI evaluator is:

```powershell
python -m crispr_phage_predictor.ml.train_classifier data/training/repeats_cas_types.csv
```

It currently:

- loads the local repeat/Cas training CSV
- filters to high-confidence labels by default
- creates a train/test split
- evaluates nearest-repeat similarity and random forest methods
- prints accuracy, a classification report, and a confusion matrix for each method

Use `--include-medium-confidence` to evaluate all validated rows.

This is an initial sanity-check evaluator. Publication-quality evaluation still needs additional methods, repeated splits, and genome/species/genus holdout splitting.

## PAM/PFS Plan

After subtype prediction exists, add a curated PAM/PFS rule table keyed by subtype. PAM support should be reported as evidence categories such as compatible, weak, absent, or not evaluated rather than as definitive resistance.
