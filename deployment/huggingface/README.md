---
title: SABR
colorFrom: blue
colorTo: green
sdk: docker
app_port: 7860
---

# SABR

SABR is a Streamlit public-beta app for CRISPR spacer-targeting evidence and
repeat-derived CRISPR-Cas subtype prediction.

Do not upload sensitive unpublished genomes to a public hosted demo. For private
or large lab analyses, run the same Docker image locally.

## Contact

For public-demo problems, bug reports, or user support, contact:

Esber Saba, Ph.D.  
The Phage Lab, Faculty of Medicine, American University of Beirut  
`es60@aub.edu.lb`

## Required Space Variables

Set these variables in the Hugging Face Space settings before launch:

- `SABR_MODEL_URL`:
  `https://zenodo.org/records/20737961/files/cas_subtype_extratrees.joblib?download=1`
- `SABR_MODEL_SHA256`:
  `ae0e8a5d56d6a2a4eb8206fb916fd9dee51f7fb9276528346b0cb8279b76cd32`

If `SABR_MODEL_URL` is unset, the app still starts but uses its documented
missing-artifact fallback instead of the selected ExtraTrees model.

## Citation

For article-linked use, cite the Zenodo DOI for the frozen SABR source/model
release: `10.5281/zenodo.20737961`. This Space is the live demo and may change
after the archived release.
