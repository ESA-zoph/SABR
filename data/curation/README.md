# Proof-of-Concept Curation

This folder tracks scientifically supported bacteria-phage pairs for SABR
benchmarking. The goal is to separate:

- pairs with experimentally supported CRISPR-mediated resistance
- pairs with spacer/protospacer evidence but known biological caveats
- phage-host systems suitable for future acquisition of derived resistant
  strains or BIMs
- records that should be excluded or verified before benchmarking

Do not treat every accession in these manifests as a final benchmark label.
Rows marked `challenge_caution` or `needs_verification` are included because
they are useful for testing whether SABR avoids overclaiming spacer matches as
confirmed resistance.

Recommended next steps:

1. Download or verify missing FASTA records listed in
   `accession_download_manifest.tsv`.
2. Deduplicate local phage files by accession and sequence hash.
3. Add source-specific notes for any experimentally derived resistant mutants
   whose complete genomes are not deposited.
4. Add PAM/seed checks before using spacer hits as resistance evidence.

## Benchmark Label Schema

`benchmark_labels.tsv` is the stricter benchmark table for SABR score
calibration and validation. It separates:

- observed phenotype label
- CRISPR-mediated resistance label
- CRISPR evidence level
- PAM evidence level
- anti-CRISPR status
- host-range status
- expected SABR behavior
- curation confidence

Use `proof_of_concept_pairs.tsv` for broader notes and early candidate tracking.
Use `benchmark_labels.tsv` only when the row can be assigned explicit labels and
caveats. Validate the table with:

```bash
python -m pytest tests/test_benchmark.py
```
