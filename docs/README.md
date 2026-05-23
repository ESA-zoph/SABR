# SABR Documentation Outputs

This folder contains generated project documentation for the SABR model and
dataset development track.

Current generated report:

```text
docs/SABR_model_development_report.docx
```

Figures used by the report are stored in:

```text
docs/figures/
```

Regenerate the Word report and figures after dataset or model changes:

```bash
python docs/generate_model_figures.py
python docs/generate_project_report.py
```

The report is intentionally a living document. It separates manually curated
gold seed rows from computational candidate rows, documents the current model
plan, and records why random row-level accuracy is only a smoke test rather
than publication-grade validation.
