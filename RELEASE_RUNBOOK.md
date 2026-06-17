# SABR v0.1-beta Article Release Runbook

This runbook defines the recommended publication release path:

1. GitHub is the development repository.
2. Zenodo is the citable frozen software/model archive.
3. Hugging Face Spaces is the public browser demo.

## 1. Finalize Public Metadata

Before tagging the release, replace placeholders in:

- `CITATION.cff`
- `.zenodo.json`
- `README.md`
- `DEPLOYMENT.md`
- `deployment/huggingface/README.md`

Required decisions:

- GitHub repository URL.
- Hugging Face Space URL.
- Zenodo DOI: `10.5281/zenodo.20737961`.
- Author names, affiliations, support contact email, and ORCID IDs where available.
- Open-source license identifier: `MIT`.

The selected public-beta software license is MIT. The current copyright and
support contact are assigned to Esber Saba, Ph.D., The Phage Lab, Faculty of
Medicine, American University of Beirut, `es60@aub.edu.lb`. Confirm this text
before tagging the release.

## 2. Confirm Runtime Model Artifact

Current selected model artifact:

```text
models/cas_subtype_extratrees.joblib
```

Expected SHA-256:

```text
ae0e8a5d56d6a2a4eb8206fb916fd9dee51f7fb9276528346b0cb8279b76cd32
```

The model card is `MODEL_CARD.md`. Keep this metadata with the archived release.

## 3. Create Clean Release

The deployable release should include:

- `app.py`
- `crispr_phage_predictor/`
- `assets/`
- `data/examples/real_demo/`
- `requirements.txt`
- `requirements-ml.txt`
- `requirements-external.txt`
- `Dockerfile`
- `.dockerignore`
- `README.md`
- `DEPLOYMENT.md`
- `RELEASE_CHECKLIST.md`
- `MODEL_CARD.md`
- `CITATION.cff`
- `.zenodo.json`
- tests needed for public maintenance

Do not include:

- `.venv/`
- `outputs/`
- raw CCTyper runs
- imported validation zip files
- large local example genomes
- generated manuscript assets unless intentionally archived as article material
- private or unpublished genomes

## 4. Verify Before Tagging

Run:

```bash
python -m py_compile app.py
python -m pytest --basetemp=.pytest-tmp -p no:cacheprovider
python scripts/ensure_model_artifact.py
streamlit run app.py
```

Then manually smoke-test:

- `data/examples/real_demo/bacteria/PA14_CRISPR_region.fasta`
- `data/examples/real_demo/phages/JBD18_positive_control.fasta`
- `data/examples/real_demo/phages/Lambda_negative_control.fasta`

Confirm:

- upload diagnostics parse both files;
- `Run SABR analysis` shows progress under the button;
- heatmap renders;
- hit details render without resetting the app;
- outputs are written under `outputs/runs/`;
- the public-demo privacy warning is visible.

## 5. Archive On Zenodo

Use Zenodo as the article citation target:

Recommended sequence:

1. Create the public GitHub repository and push the cleaned release branch.
2. Connect the repository in Zenodo's GitHub integration.
3. Create a GitHub release/tag such as `v0.1-beta`.
4. Confirm Zenodo ingests the GitHub release archive.
5. Upload `models/cas_subtype_extratrees.joblib` to the same Zenodo record, or
   create a linked software/data record for the runtime artifact.
6. Confirm the DOI, version metadata, and artifact SHA-256.
7. Update `README.md`, `MODEL_CARD.md`, `CITATION.cff`, `.zenodo.json`, and the
   article with the DOI after Zenodo reserves or publishes it.

Current Zenodo DOI:

```text
10.5281/zenodo.20737961
```

## 6. Deploy Hugging Face Space

Create a Docker Space using the contents of this repository.

Set Space variables:

```text
SABR_MODEL_URL=https://zenodo.org/records/20737961/files/cas_subtype_extratrees.joblib?download=1
SABR_MODEL_SHA256=ae0e8a5d56d6a2a4eb8206fb916fd9dee51f7fb9276528346b0cb8279b76cd32
```

Use `deployment/huggingface/README.md` as the Space README template.

After deployment, run the same built-in real demo smoke test in the browser.

The Docker Space README must include `sdk: docker` and `app_port: 7860` in its
YAML front matter. Runtime variables are configured in the Space settings, not
committed into the repository.

## 7. Cite In The Article

Use:

- Zenodo DOI: `10.5281/zenodo.20737961` for the frozen software/model release.
- Hugging Face Space URL: live public demo.
- GitHub URL: development repository and issue tracker.

Phrase the tool as CRISPR spacer-targeting evidence mapping with repeat-derived
subtype and PAM/PFS support, not as a definitive phage resistance predictor.
