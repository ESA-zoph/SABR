# CRISPR-Phage Resistance Predictor

Early-stage bioinformatics tool for estimating hypothetical CRISPR-mediated bacterial resistance against bacteriophages.

## Current Status
This repository is a scaffold. The first goal is to build a reproducible GUI pipeline that can:

- accept multiple bacterial FASTA files
- accept multiple phage FASTA files
- parse and validate sequences
- show upload diagnostics before running heavier analysis
- detect candidate exact-repeat CRISPR arrays
- extract candidate spacers
- match extracted spacers exactly against phage genomes
- show a bacteria-by-phage heatmap as the main result
- prepare Cas typing, PAM analysis, and scoring modules

See `context.md` for the scientific scope and project decisions.

## Planned Workflow
1. Detect CRISPR arrays in bacterial genomes.
2. Extract repeats and spacers.
3. Cross all bacterial spacers against all phage genomes.
4. Predict CRISPR-Cas type/subtype using repeat, array, and cas-gene features.
5. Select candidate PAM rules based on predicted system type.
6. Produce bacteria-by-phage resistance likelihood scores.

## Install
Create a virtual environment, then install dependencies:

```bash
pip install -r requirements.txt
```

Machine-learning dependencies are kept separate so the GUI remains easier to deploy:

```bash
pip install -r requirements-ml.txt
```

## Train Baseline Cas-Subtype Classifier

After creating a repeat/Cas training CSV at `data/training/repeats_cas_types.csv`, run:

```bash
python -m crispr_phage_predictor.ml.train_classifier data/training/repeats_cas_types.csv
```

This trains the first baseline repeat/array classifier and prints accuracy, a classification report, and a confusion matrix. This is an early sanity check; final evaluation should use genome/species/genus holdout splits.

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

## Current CRISPR Detection Method
The first detector is a transparent MVP baseline. It identifies exact direct repeats of CRISPR-like length separated by plausible spacer lengths. This is useful for early development and benchmarking, but it is not yet a replacement for established CRISPR callers.

## Current Spacer Matching Method
The first matcher searches every extracted spacer against every uploaded phage genome on both strands. It reports exact spacer-protospacer hits and a bacteria-by-phage evidence matrix. Approximate matching, PAM analysis, and Cas-type-aware scoring are planned next.

## Important Scientific Note
The tool should report evidence-based hypothetical CRISPR targeting or resistance likelihood, not definitive biological resistance. Experimental validation remains necessary.
