# SABR Public-Beta Release Checklist

## Ready

- Streamlit GUI accepts bacterial and phage FASTA uploads.
- FASTA parsing and upload diagnostics are implemented.
- Internal exact-repeat CRISPR detection is implemented.
- Spacer extraction and exact spacer-phage matching are implemented.
- Optional BLASTN matching path exists when BLAST+ is installed.
- Repeat/array-derived Cas subtype prediction is implemented.
- PAM/PFS-aware scoring is implemented for the current curated rule subset.
- Output tables and run metadata are exported under `outputs/runs/`.
- Manuscript draft documents cautious targeting-evidence framing.
- Selected ExtraTrees model has documented internal and CCTyper-supported validation.
- Docker deployment files are present for local and Hugging Face Spaces-style runs.
- Docker startup can download the runtime model artifact and validate its SHA-256 hash.
- A lightweight real PA14/JBD18/Lambda demo panel is available under `data/examples/real_demo/`.
- 2026-06-17 release-readiness checks passed:
  - compileall for `app.py`, `crispr_phage_predictor`, `scripts`, and `docs`
  - local model artifact presence check
  - phage-host database audit with `issues 0`
  - full pytest suite with cache disabled: `111 passed`
  - Docker image build for `sabr:local`
  - Docker container HTTP check at `http://127.0.0.1:7860`
  - built-in real demo pipeline check: PA14/JBD18 strong evidence and PA14/Lambda no evidence

## Must Do Before Public Beta

- Smoke-test `streamlit run app.py` in a fresh virtual environment.
- Choose the public GitHub repository owner/name and update `CITATION.cff` and
  `.zenodo.json`.
- Choose the Hugging Face Space owner/name and update `CITATION.cff`,
  `.zenodo.json`, and `deployment/huggingface/README.md` if needed.
- Confirm public contact/support address is correct in README, deployment docs,
  Space README, `CITATION.cff`, and `.zenodo.json`.
- Archive the frozen source release plus runtime model artifact on Zenodo.
- After Zenodo publication, replace all example model URLs with the final
  DOI-backed model file URL.
- Configure Hugging Face Space variables:
  - `SABR_MODEL_URL`
  - `SABR_MODEL_SHA256`
- Add an example output folder or screenshots with no large raw data.
- Review every user-facing label for targeting-evidence terminology.
- Confirm the app behavior when `models/cas_subtype_extratrees.joblib` is missing.
- Remove or ignore bulky local research artifacts from the release package.

## Should Do Soon After Beta

- Add CI for linting and tests.
- Add a small command-line batch mode for users who do not want Streamlit.
- Improve PAM/PFS subtype coverage.
- Benchmark internal CRISPR detection against established callers.
- Expand independent III-B and III-D validation.

## Not Ready To Claim

- Phenotype-level phage resistance prediction.
- Complete CRISPR-Cas system functionality.
- Exhaustive PAM/PFS support.
- Strong Type III-B external subtype validation.
- Large-genome or large-batch web-hosted scalability without resource testing.
