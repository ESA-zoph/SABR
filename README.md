# SABR CRISPR Targeting Evidence Mapper

SABR, Spacer Alignment-Based Recognition, is a bioinformatics tool for mapping candidate CRISPR spacer-targeting evidence against bacteriophages and inferring likely CRISPR-Cas type/subtype from repeat and array features.

SABR is an evidence mapper, not a direct phage-resistance caller. Spacer matches, repeat-derived subtype calls, PAM/PFS checks, and scores should be interpreted as computational evidence layers that require biological context and experimental validation before resistance claims.

## Current Status
This repository is an active public-beta research tool. The current GUI pipeline can:

- accept multiple bacterial FASTA files
- accept multiple phage FASTA files
- parse and validate sequences
- show upload diagnostics before running heavier analysis
- load a small real built-in PA14/JBD18/Lambda demo panel for first-time users
- detect candidate exact-repeat CRISPR arrays
- extract candidate spacers
- match extracted spacers exactly against phage genomes
- show a bacteria-by-phage spacer-targeting evidence heatmap as the main result
- label each bacteria-phage pair as no, weak, moderate, or strong candidate targeting evidence
- inspect heatmap cells with spacer/protospacer alignments, PAM/PFS context, seed summaries, and score components
- predict likely CRISPR-Cas subtype from FASTA-derived repeat/array features when the trained model artifact is available
- evaluate subtype-aware PAM/PFS support for the curated rule subset currently encoded in SABR
- save reproducible local run outputs and download a Markdown analysis report

See `context.md` for the scientific scope and project decisions.

## Workflow
1. Detect CRISPR arrays in bacterial genomes.
2. Extract repeats and spacers.
3. Cross all bacterial spacers against all phage genomes.
4. Predict CRISPR-Cas type/subtype using repeat and array features.
5. Select candidate PAM rules based on predicted system type.
6. Produce bacteria-by-phage CRISPR targeting evidence scores.

## Positioning

SABR is intended to complement, not replace, established CRISPR annotation and
phage-host prediction tools. Dedicated tools such as CRISPR array callers,
CRISPR-Cas locus typers, spacer-target search servers, and large-scale
phage-host predictors each address important parts of the workflow. SABR's
public-beta niche is an integrated panel-analysis interface: users can upload
bacterial and phage FASTA files together, obtain a bacteria-by-phage targeting
evidence matrix, inspect spacer/protospacer support, and interpret hits with
repeat-derived subtype and PAM/PFS evidence in one reproducible run.

The current publishable contribution is therefore the integrated evidence
mapping workflow and validation of its repeat-derived subtype component, not a
claim of definitive phage resistance prediction.

## Contact

For bugs, problems using SABR, or reports about the public demo, contact:

Esber Saba, Ph.D.  
The Phage Lab, Faculty of Medicine, American University of Beirut  
`es60@aub.edu.lb`

## Install For Local Use
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

## Model Artifact

The deployed subtype predictor expects a trained model artifact at:

```text
models/cas_subtype_extratrees.joblib
```

This artifact is intentionally ignored by Git because it is large. For local beta use, keep the artifact in `models/`. If the artifact is absent, SABR falls back to a nearest-repeat training-table baseline when a local training table is available, but public deployment should ship or download a frozen model artifact with documented metadata.

Current documented selected model:

- architecture: flat ExtraTrees
- training table: `data/training/repeats_cas_types_crisprcasdb_sql_candidate.csv`
- training rows after minimum-class filtering: 23,478
- artifact SHA-256: `ae0e8a5d56d6a2a4eb8206fb916fd9dee51f7fb9276528346b0cb8279b76cd32`
- internal genome-held-out accuracy: 0.9455
- independent enriched strict CCTyper validation: 23/25 expected-subtype rows correct, accuracy 0.9200
- unresolved validation gap: Type III-B still lacks an evaluable expected-subtype array-level validation row

See `MODEL_CARD.md` for artifact metadata, intended use, validation summary,
and public-release distribution guidance.

## Docker

Build and run the local public-beta container:

```bash
docker build -t sabr:local .
docker run --rm -p 7860:7860 \
  -e SABR_MODEL_URL="https://example.org/cas_subtype_extratrees.joblib" \
  -e SABR_MODEL_SHA256="ae0e8a5d56d6a2a4eb8206fb916fd9dee51f7fb9276528346b0cb8279b76cd32" \
  sabr:local
```

Open `http://127.0.0.1:7860`.

The Docker image intentionally excludes `models/*.joblib`. Set
`SABR_MODEL_URL` and `SABR_MODEL_SHA256`, or mount
`models/cas_subtype_extratrees.joblib` into `/app/models/` when running locally.
If no artifact is provided, SABR starts with its documented missing-artifact
fallback.

## Deployment Strategy

Recommended article-linked public-beta strategy:

1. Archive the frozen `v0.1-beta` source release and runtime model artifact on
   Zenodo so the article can cite a DOI.
2. Host a public browser demo on Hugging Face Spaces using the Dockerfile in
   this repository.
3. Configure the Space with `SABR_MODEL_URL` pointing to the archived model
   artifact and `SABR_MODEL_SHA256` set to the hash above.
4. Keep large CCTyper runs, raw validation packages, generated manuscript
   artifacts, local outputs, and private data outside the deployable app
   package.
5. Use targeting-evidence terminology throughout. SABR is not a definitive
   phage-resistance caller.

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

After uploading files, first check the upload diagnostics table. It should show at least one parsed record for each bacterial and phage file. Then press **Run SABR analysis** in the sidebar.

Each completed GUI analysis writes a timestamped local output folder under
`outputs/runs/`. The saved files include run metadata, FASTA record summaries,
candidate arrays, extracted spacers, spacer hits, the evidence matrix, and the
heatmap table.

The GUI supports optional external backends when they are available:

- CRISPR detection: auto recommended, internal exact-repeat MVP, or MinCED-compatible
- Spacer-phage matching: auto recommended, internal exact match, or BLASTN

The Streamlit selector defaults to the internal exact-repeat detector and
internal exact spacer matcher so the built-in real demo is portable and
predictable. Auto recommended uses MinCED-compatible detection when the Diced
Python package or a `minced` command is available, and otherwise falls back to
the internal exact-repeat detector. For spacer-phage matching, auto recommended
uses BLASTN when `blastn` and `makeblastdb` are available and otherwise falls
back to the internal exact matcher. The internal methods remain available as
reproducible baselines.

When BLASTN is active, the sidebar exposes a minimum identity threshold and a
full-spacer-alignment requirement. BLAST-derived hit tables include identity,
alignment length, spacer coverage, e-value, and bitscore.

## Current CRISPR Detection Method
The first detector is a transparent MVP baseline. It identifies exact direct repeats of CRISPR-like length separated by plausible spacer lengths. This is useful for early development and benchmarking, but it is not yet a replacement for established CRISPR callers.

## Current Spacer Matching Method
The internal matcher searches every extracted spacer against every uploaded phage genome on both strands. It reports exact spacer-protospacer hits and a bacteria-by-phage evidence matrix. BLASTN matching can be selected when BLAST+ is available. PAM/PFS-aware scoring is implemented for the current curated subtype rule subset and should be treated as partial biological evidence, not functional proof.

## Important Scientific Note
The tool reports evidence-based candidate CRISPR targeting, not definitive biological resistance. Spacer matches, repeat-derived Cas subtype calls, PAM/PFS support, and seed summaries are useful evidence layers, but experimental validation remains necessary before claiming resistance.
