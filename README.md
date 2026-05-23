# SABR CRISPR Targeting Evidence Mapper

Early-stage bioinformatics tool for mapping candidate CRISPR spacer targeting evidence against bacteriophages and inferring likely CRISPR-Cas type/subtype from repeat and array features.

## Current Status
This repository is a scaffold. The first goal is to build a reproducible GUI pipeline that can:

- accept multiple bacterial FASTA files
- accept multiple phage FASTA files
- parse and validate sequences
- show upload diagnostics before running heavier analysis
- detect candidate exact-repeat CRISPR arrays
- extract candidate spacers
- match extracted spacers exactly against phage genomes
- show a bacteria-by-phage spacer-targeting evidence heatmap as the main result
- prepare Cas typing, PAM analysis, and scoring modules

See `context.md` for the scientific scope and project decisions.

## Planned Workflow
1. Detect CRISPR arrays in bacterial genomes.
2. Extract repeats and spacers.
3. Cross all bacterial spacers against all phage genomes.
4. Predict CRISPR-Cas type/subtype using repeat, array, and cas-gene features.
5. Select candidate PAM rules based on predicted system type.
6. Produce bacteria-by-phage CRISPR targeting evidence scores.

## Install
Create a virtual environment, then install dependencies:

```bash
pip install -r requirements.txt
```

Machine-learning dependencies are kept separate so the GUI remains easier to deploy:

```bash
pip install -r requirements-ml.txt
```

Optional external-backend dependencies are also separate:

```bash
pip install -r requirements-external.txt
```

## Train Baseline Cas-Subtype Classifier

After creating a repeat/Cas training CSV at `data/training/repeats_cas_types.csv`, run:

```bash
python -m crispr_phage_predictor.ml.train_classifier data/training/repeats_cas_types.csv
```

This compares the nearest-repeat baseline and the first random-forest repeat/array classifier, then prints accuracy, classification reports, and confusion matrices. This is an early sanity check; final evaluation should use genome/species/genus holdout splits.

## CCTyper Annotation Environment

For generating training labels from CCTyper, use a Linux/WSL conda environment:

```bash
conda env create -f envs/cctyper-linux.yml
conda activate crispr-cctyper
```

Then run CCTyper on genomes and collect the output folders using the training-data commands documented in `data/training/README.md`.

## Run GUI

```bash
streamlit run app.py
```

On Windows after dependencies are installed, you can also run:

```bat
run_app.bat
```

Then open:

```text
http://127.0.0.1:8501
```

After uploading files, first check the upload diagnostics table. It should show at least one parsed record for each bacterial and phage file. Then press **Run CRISPR-phage analysis** in the sidebar.

Each completed GUI analysis writes a timestamped local output folder under
`outputs/runs/`. The saved files include run metadata, FASTA record summaries,
candidate arrays, extracted spacers, spacer hits, the evidence matrix, and the
heatmap table.

The GUI supports optional external backends when they are available:

- CRISPR detection: auto recommended, internal exact-repeat MVP, or MinCED-compatible
- Spacer-phage matching: auto recommended, internal exact match, or BLASTN

Auto recommended currently keeps CRISPR detection on the internal detector for
responsive whole-genome uploads. MinCED-compatible detection can be selected
manually for benchmarking when the Diced Python package or a `minced` command is
available. For spacer-phage matching, auto recommended uses BLASTN when `blastn`
and `makeblastdb` are available and otherwise falls back to the internal exact
matcher. The internal methods remain available as reproducible baselines.

When BLASTN is active, the sidebar exposes a minimum identity threshold and a
full-spacer-alignment requirement. BLAST-derived hit tables include identity,
alignment length, spacer coverage, e-value, and bitscore.

## Current CRISPR Detection Method
The first detector is a transparent MVP baseline. It identifies exact direct repeats of CRISPR-like length separated by plausible spacer lengths. This is useful for early development and benchmarking, but it is not yet a replacement for established CRISPR callers.

## Current Spacer Matching Method
The first matcher searches every extracted spacer against every uploaded phage genome on both strands. It reports exact spacer-protospacer hits and a bacteria-by-phage evidence matrix. Approximate matching, PAM analysis, and Cas-type-aware scoring are planned next.

## Important Scientific Note
The tool reports evidence-based candidate CRISPR targeting, not definitive biological resistance. Spacer matches, repeat-derived Cas subtype calls, PAM/PFS support, and seed summaries are useful evidence layers, but experimental validation remains necessary before claiming resistance.
