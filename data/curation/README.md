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

## Phage-Host Interaction Schema

`phage_host_interactions.tsv` is the starting table for sensitivity/resistance
prediction. It should capture experimental bacteria-phage interaction outcomes
before they are converted into SABR training rows. Use one row per measured
strain-phage-assay result, not one row per paper.

Important fields:

- `raw_eop`: the value reported by the source, including inequalities such as
  `<1e-6`.
- `eop_value`: normalized numeric value when extractable.
- `eop_class`: normalized class using the first SABR thresholds:
  `high >= 0.5`, `medium >= 0.1`, `low >= 1e-3`, `trace < 1e-3`,
  and `none` for zero or below-detection values.
- `susceptibility_label`: phenotype label derived from the assay when possible:
  `susceptible`, `reduced_susceptibility`, `resistant`, `nonhost`, `mixed`, or
  `unknown`.
- `anti_crispr_status` and `anti_crispr_genes`: phage counter-defense evidence
  that can explain why a spacer/PAM hit does not produce resistance.
- `crispr_interference_evidence`: whether the phenotype is experimentally tied
  to CRISPR interference, inferred from spacer/PAM/seed evidence, or not
  evaluated.

Validate and normalize this schema with:

```bash
python -m pytest tests/test_interactions.py
```

### Current imported EOP sources

- `Magadan2012_CRISPR3`: hand-curated DGCC7710/phage 2972 BIM EOP rows
  from Table 3.
- `Cady2012_PA14_CRISPR`: hand-curated PA14/JBD plaque-assay rows where
  numeric EOP was not reported.
- `MirzaeiNilsson2015_EOP`: programmatically imported EOP rows from S1 Table
  using `scripts/import_mirzaei_nilsson_2015_eop.py`. The downloaded source
  file is `mirzaei_nilsson_2015_s1.docx`, and the derived import table is
  `mirzaei_nilsson_2015_eop.tsv`.
- `Sofy2021_KPP5`, `Latka2021_K23Klebsiella`, `Singh2025_KpCocktail`,
  `Abdelhadi2021_STM2`, `Wintachai2022_vABWU2101`, and
  `Wintachai2024_vECPW8`, `Kim2022_KP1_KP12`, `Ciacci2018_vB_Kpn_F48`,
  and `Peng2025_vB_Kp_XP4`: imported through
  `scripts/import_embedded_literature_eop.py` into
`embedded_literature_eop.tsv`.

## Accession Linkage

`accession_linkage.tsv` tracks whether interaction-table bacteria and phages
can be connected to genome accessions for feature extraction. Keep this separate
from the phenotype table because many EOP studies report clinical/panel isolate
names without deposited genome assemblies.

Run:

```bash
python scripts/report_accession_linkage_coverage.py
```

Current linkage interpretation:

- `phage_genome_linked`: enough information exists to extract phage genome
  features, anti-CRISPR candidates, tail-fiber/baseplate annotations, and
  phage k-mer features.
- `bacterium_genome_linked`: enough information exists to extract host CRISPR,
  PAM/spacer matching, and defense-system features.
- `pair_genome_ready`: both sides are linked; these rows can support full
  SABR susceptibility/resistance model features.

The first linkage pass intentionally resolves phage accessions before host
clinical-panel accessions, because one phage genome often unlocks many measured
interaction rows. The largest remaining bottleneck is exact host-strain
genome resolution for ECOR/ESBL/SARA/SARB and clinical Klebsiella panels.

Hybrid mode is now represented with `linkage_status = reference_proxy`.
Reference proxies are downloaded public genomes used only for provisional
species-level host features. They do not replace exact isolate genomes. Coverage
reports therefore distinguish:

- `pair_genome_ready`: exact/strain-alias host genome plus exact phage genome.
- `pair_hybrid_ready`: exact phage genome plus exact/strain-alias/reference
  proxy host genome.

Download linked public genomes with:

```bash
python scripts/download_accession_linkage_genomes.py
```

For exact host-strain resolution, generate a review table with:

```bash
python scripts/resolve_host_accession_candidates.py
```

This writes `host_accession_candidates.tsv`. Candidate rows are intentionally
not promoted to `accession_linkage.tsv` automatically because NCBI nucleotide
searches often return plasmids, partial genes, or WGS master records rather
than exact chromosome/genome assemblies.

For exact host assemblies, use the Assembly-based promoter:

```bash
python scripts/promote_host_assembly_links.py --dry-run
python scripts/promote_host_assembly_links.py
```

This writes `host_assembly_candidates.tsv`, verifies exact strain identity from
the NCBI assembly report, downloads accepted genome FASTA files, and appends
accepted rows to `accession_linkage.tsv`. The current promoted set adds 71
exact host assemblies and raises strict exact-pair genome-ready coverage to 74
rows.

Run the database QA audit before using the interaction database for modeling or
manuscript claims:

```bash
python scripts/audit_phage_host_database.py
```

The audit validates schemas, EOP class consistency, accession-linkage paths,
downloaded FASTA metadata, coverage/feature-table alignment, and promoted host
assembly evidence.
