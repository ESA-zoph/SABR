# CRISPRCasdb 34 Raw Source

Local raw files moved here from the temporary `crisprdb/` folder:

- `dr_34.zip`: direct-repeat FASTA exports.
- `spacer_34.zip`: spacer FASTA exports.
- `ccpp_db.zip`: PostgreSQL dump containing CRISPRCasdb-style tables.
- `dbschem.pdf`: database schema reference.

Observed archive contents:

- `dr_34.zip`
  - `direct_repeat_id.fsa`
  - `direct_repeat_seqName.fsa`
  - `direct_repeat_taxon.fsa`
- `spacer_34.zip`
  - `spacer_id.fsa`
  - `spacer_seqName.fsa`
  - `spacer_taxon.fsa`
- `ccpp_db.zip`
  - `home/pa.charbit/20220414_ccpp_recette_chromo_complete.sql`

Quick inventory from local inspection:

- Direct repeats in `direct_repeat_id.fsa`: 28,712 FASTA records.
- Spacers in `spacer_id.fsa`: 353,377 FASTA records.
- SQL dump includes tables such as `crisprlocus`, `crisprlocus_region`,
  `region`, `sequence`, `taxon`, `clustercas`, and `clustercas_gene`.

Expected SABR value:

- Good source for expanding repeat-sequence coverage and nearest-repeat
  baselines.
- Potential source for computational candidate Cas type/subtype labels if the
  SQL dump can connect repeats or loci to nearby `clustercas` records.
- Useful for audit/provenance and comparison against the current
  CRISPRCasdb-derived Vink import path.

Caveats:

- Treat as raw computational source data, not gold-standard training labels.
- Do not use taxonomy, accession, or database-only fields as runtime model
  features.
- Before training, create a derived table that records source file, accession or
  sequence ID, repeat sequence, label source, label confidence, and filtering
  criteria.
- Large raw files are intentionally ignored by git.
