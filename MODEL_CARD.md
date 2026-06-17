# SABR Runtime Cas-Subtype Model Card

## Artifact

- File name: `cas_subtype_extratrees.joblib`
- Expected runtime path: `models/cas_subtype_extratrees.joblib`
- Size: `410202588` bytes
- SHA-256: `ae0e8a5d56d6a2a4eb8206fb916fd9dee51f7fb9276528346b0cb8279b76cd32`

## Intended Use

This artifact supports SABR's repeat/array-derived CRISPR-Cas subtype prediction
inside the CRISPR spacer-targeting evidence workflow. It is an evidence layer for
choosing and reporting candidate subtype-aware PAM/PFS support. It is not a
standalone phage-resistance predictor.

## Training Metadata

- Model method: ExtraTrees
- Training table: `data/training/repeats_cas_types_crisprcasdb_sql_candidate.csv`
- Training rows after minimum-class filtering: `23478`
- Minimum class count: `20`
- Random state: `42`
- Number of estimators: `400`
- Classes: `I-A`, `I-B`, `I-C`, `I-D`, `I-E`, `I-F`, `I-G`, `II-A`,
  `II-B`, `II-C`, `III-A`, `III-B`, `III-C`, `III-D`, `V-A`, `V-K`,
  `VI-B1`

## Validation Summary

- Internal genome-held-out accuracy documented for the selected model: `0.9455`.
- Independent enriched strict CCTyper validation: `23/25` expected-subtype rows
  correct, accuracy `0.9200`.
- Known external validation gap: Type `III-B` remains insufficiently tested at
  array level.

## Distribution Plan

For article-linked public releases, archive this model artifact with the frozen
SABR software release on Zenodo and use the Zenodo-hosted artifact URL as
`SABR_MODEL_URL` for the Hugging Face Spaces demo. Configure:

```text
SABR_MODEL_URL=https://zenodo.org/records/20737961/files/cas_subtype_extratrees.joblib?download=1
SABR_MODEL_SHA256=ae0e8a5d56d6a2a4eb8206fb916fd9dee51f7fb9276528346b0cb8279b76cd32
```

Archived release DOI: `10.5281/zenodo.20737961`.

The Docker startup script validates the SHA-256 before placing the artifact at
the runtime model path.

## Limitations

- The model predicts likely subtype from repeat/array features; it does not
  verify full CRISPR-Cas locus functionality.
- CRISPRCasdb SQL-derived labels are computational candidates rather than final
  gold-standard labels.
- Weak/confusable regions include Type III-heavy boundaries, especially
  `III-B`, `III-D`, and neighboring Type I-B/I-C/I-A cases.
- PAM/PFS support in SABR is based on the curated rule subset currently encoded
  in the app and should be interpreted as partial biological evidence.
