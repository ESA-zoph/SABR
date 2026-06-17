# SABR Deployment Notes

SABR should be released first as a local Streamlit public beta. The scientific
positioning is important: SABR reports CRISPR spacer-targeting evidence and
repeat-derived Cas subtype predictions, not definitive phage resistance.

## Recommended Release Shape

Use a small deployable repository containing:

- `app.py`
- `crispr_phage_predictor/`
- `assets/`
- `models/cas_subtype_extratrees.joblib` or documented download instructions
- `requirements.txt`
- `requirements-ml.txt`
- `requirements-external.txt`
- `README.md`
- `DEPLOYMENT.md`
- example FASTA files and expected output snapshots
- tests needed for public maintenance

Keep bulky research/provenance material out of the deployable app package:

- raw CCTyper output folders
- large validation zip files
- generated manuscript figures unless needed for documentation
- training source dumps
- local `outputs/runs/`

## Local Public-Beta Run

```bash
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
pip install -r requirements-ml.txt
streamlit run app.py
```

Open:

```text
http://127.0.0.1:8501
```

Optional external matching/detection tools can be installed separately:

```bash
pip install -r requirements-external.txt
```

BLAST+ and CCTyper/MinCED-compatible tooling are system-level dependencies, not
ordinary Python-only dependencies.

## Docker Run

The repository includes a Dockerfile for local Streamlit deployment and Hugging
Face Spaces-style hosting. The image does not copy `models/*.joblib` by default.
Provide the runtime artifact at launch through environment variables:

```bash
docker build -t sabr:local .
docker run --rm -p 7860:7860 \
  -e SABR_MODEL_URL="https://example.org/cas_subtype_extratrees.joblib" \
  -e SABR_MODEL_SHA256="ae0e8a5d56d6a2a4eb8206fb916fd9dee51f7fb9276528346b0cb8279b76cd32" \
  sabr:local
```

Open:

```text
http://127.0.0.1:7860
```

If `SABR_MODEL_URL` is not set, the container still starts and the app displays
the missing-artifact fallback warning.

To use a local model artifact instead of downloading it, mount it into the
container:

```bash
docker run --rm -p 7860:7860 \
  -v "$PWD/models/cas_subtype_extratrees.joblib:/app/models/cas_subtype_extratrees.joblib:ro" \
  sabr:local
```

## Hugging Face Spaces

Use the Docker SDK and port `7860`. A Space README template is provided at
`deployment/huggingface/README.md`; copy its YAML front matter into the Space
repository README when creating the hosted demo.

Configure Space variables:

- `SABR_MODEL_URL`
- `SABR_MODEL_SHA256`

For public demos, include a visible warning that users should not upload
sensitive unpublished genomes. Recommend local Docker for private lab analysis.

## Article-Linked Release Plan

For publication, use Zenodo as the primary citable archive and Hugging Face
Spaces as the interactive demo:

1. Create a clean `v0.1-beta` source release from the deployable repository.
2. Archive the source release on Zenodo to obtain a DOI.
3. Upload `models/cas_subtype_extratrees.joblib` with the Zenodo release or as a
   linked Zenodo dataset/software artifact.
4. Set the Hugging Face Space variables:
   - `SABR_MODEL_URL`: stable HTTPS URL for `cas_subtype_extratrees.joblib`
   - `SABR_MODEL_SHA256`:
     `ae0e8a5d56d6a2a4eb8206fb916fd9dee51f7fb9276528346b0cb8279b76cd32`
5. In the article, cite the Zenodo DOI for the exact software/model release and
   provide the Hugging Face Space URL as a live demo.

The development GitHub repository can remain the working project and issue
tracker, but the DOI-backed Zenodo release should be the citation target.

## Model Artifact

The preferred runtime model path is:

```text
models/cas_subtype_extratrees.joblib
```

The current artifact is large, so do not commit it to ordinary Git history.
For article-linked distribution, prefer Zenodo DOI-backed archival storage for
the source release and the runtime model artifact. GitHub Releases,
institutional storage, or a Hugging Face model repository can be used as mirrors,
but the paper should cite the DOI-backed frozen release.

Every public model artifact should include metadata:

- training table name and hash
- training row count after filtering
- feature schema version
- model class and hyperparameters
- validation metrics
- known limitations, especially Type III-B validation

## Public-Beta Caveats

User-facing documentation must state:

- SABR is not a phenotype-level resistance predictor.
- Internal CRISPR detection is a transparent exact-repeat baseline.
- BLASTN and MinCED/CCTyper support depend on optional local tools.
- PAM/PFS rules are incomplete.
- Type III-B subtype prediction remains insufficiently validated externally.
- Scores are bounded evidence scores, not probabilities of resistance.

## Support Contact

For public-demo problems, bug reports, or user support:

Esber Saba, Ph.D.  
The Phage Lab, Faculty of Medicine, American University of Beirut  
`es60@aub.edu.lb`

## Minimum Pre-Release Checks

Run:

```bash
python -m py_compile app.py
python -m pytest --basetemp=.pytest-tmp -p no:cacheprovider
streamlit run app.py
```

Then manually verify:

- app loads without missing import errors
- model artifact metadata is displayed or missing-artifact fallback is clear
- bacterial and phage FASTA uploads parse correctly
- CRISPR arrays, spacers, hits, subtype predictions, and score matrices export
- output language stays within targeting-evidence framing

## Later Deployment Extensions

The current Docker image is intended for Streamlit plus the Python-only SABR
runtime. Docker remains the preferred route for bundling BLAST+,
MinCED-compatible detection, CCTyper, and exact tool versions later. A hosted
web app is possible, but genome uploads and model artifact size require resource
limits and careful data-handling language.
