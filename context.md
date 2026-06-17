# Project Context

## Working Name
SABR

Expanded name:

- Spacer Alignment-Based Recognition

Branding direction:

- SABR should be the app-facing name.
- The interface should retain cautious scientific language such as CRISPR-phage targeting evidence, repeat-derived Cas typing, PAM/PFS support, and candidate interaction.
- The header should use AUB-inspired colors with the SABR name on one side and The Phage Lab/AUB logo on the other.
- Current logo asset: `assets/aub-logo.png`.

## Goal
Build an easy-to-deploy bioinformatics tool for the scientific community that reports CRISPR spacer-targeting evidence against bacteriophages and predicts likely CRISPR-Cas type/subtype from repeat and array features.

The tool should accept multiple bacterial FASTA files and multiple phage FASTA files, identify CRISPR arrays, extract spacers and repeats, compare spacers against phage genomes, infer possible PAM/PFS compatibility, and produce an evidence-based bacteria-by-phage CRISPR targeting matrix. Resistance should remain a later interpretation layer only when phenotype, literature, or benchmark evidence supports it.

## Current Snapshot: 2026-05-23

Current scientific framing:

- SABR is now a CRISPR spacer-targeting evidence and repeat-based Cas-typing tool, not a direct resistance caller.
- Preferred main score name: `crispr_targeting_score`.
- Historical `hypothetical_resistance_score` aliases remain only for backward compatibility with older artifacts and tests.
- Taxonomy/species/genus calling should not be inferred from CRISPR evidence; it would require a separate genome identity module later.

Current best Cas subtype model:

- Runtime artifact: `models/cas_subtype_extratrees.joblib`.
- Method: ExtraTrees on repeat/array features.
- Training table: `data/training/repeats_cas_types_crisprcasdb_sql_candidate.csv`.
- Selected SABR ExtraTrees model, trained on the CRISPRCasdb-derived
  annotation table, genome-holdout performance: accuracy 0.9455.
- Transfer to compatible rows in the current SABR table: accuracy 0.9811.
- Reciprocal current-table-to-CRISPRCasdb transfer: accuracy 0.9559.
- Weak/confusable regions remain Type III-heavy, especially `III-B`, `III-D`, and adjacent Type I-B/I-C/I-A boundaries.
- Caveat: CRISPRCasdb SQL labels are computational candidates, not final gold-standard truth.

Latest independent CCTyper validation:

- VirtualBox/CCTyper package imported from `cctyper_full_balanced_validation_package.zip`.
- Extracted under `data/validation/cctyper_full_balanced_validation_package/`.
- Reproducible evaluator: `docs/evaluate_cctyper_strict_validation.py`.
- Outputs: `data/validation/cctyper_strict_validation/`.
- Strict CCTyper screen retained 20 genomes: 6 I-B, 5 I-C, 5 III-A, 1 III-B, and 3 III-D.
- Expected-subtype array-level evaluation retained 25 rows: 8 I-B, 6 I-C, 6 III-A, and 5 III-D.
- Selected SABR ExtraTrees evaluation on these rows: 23/25 correct, accuracy 0.9200.
- Errors: two I-C rows predicted as III-B.
- Important caveat: the strict III-B genome did not contribute an expected-subtype III-B array row; Type III-B remains insufficiently tested at array level.

Current deployment/package direction:

- Release target should be `v0.1-beta` as a local Streamlit app first, not a hosted public web service yet.
- Scientific positioning must remain: CRISPR spacer-targeting evidence and repeat-derived subtype prediction, not phenotype-level phage resistance prediction.
- Added public-facing deployment docs:
  - `DEPLOYMENT.md`
  - `RELEASE_CHECKLIST.md`
- Updated `README.md` with current model-artifact expectations, validation summary, and public-beta strategy.
- Added tiny synthetic smoke-test FASTAs under `data/examples/basic_demo/`.
- Smoke-test result for those examples: 1 candidate array, 2 extracted spacers, 1 spacer-protospacer hit.
- Updated app sidebar to show whether `models/cas_subtype_extratrees.joblib` is present.
- Updated pipeline stage text so Cas subtype prediction and PAM/PFS evaluation are no longer described as only planned.
- Updated `.gitignore` to avoid accidentally packaging bulky CCTyper runs, imported validation zip files, extracted raw validation packages, local outputs, and temp/cache folders.
- Focused verification passed after cleanup:
  - `py_compile` for edited Python files.
  - 11 focused tests covering pipeline output, Cas prediction, external dataset evaluation, and CCTyper import.
- Full release-readiness verification now passed:
  - `python -m compileall -q app.py crispr_phage_predictor docs/evaluate_cctyper_strict_validation.py docs/generate_manuscript_docx.py`
  - full pytest suite: 86 passed.
  - Streamlit headless startup on port 8502 returned HTTP 200 OK.
- Runtime model I/O improvement:
  - `crispr_phage_predictor/ml/model_artifact.py` now caches model artifact loading and SHA-256 hashing with bounded `lru_cache`.
  - `save_artifact()` clears those caches to avoid stale artifacts after retraining/export.
- Main packaging blocker: the selected runtime model artifact is about 410 MB; public distribution needs Git LFS, a release asset, Zenodo/institutional storage, or Docker image packaging rather than ordinary Git.
- Next practical release steps:
  1. Run the full test suite in a clean workspace.
  2. Create a fresh virtual environment and verify install from `requirements.txt` plus `requirements-ml.txt`.
  3. Run `streamlit run app.py` and manually test the synthetic example FASTAs.
  4. Decide model artifact distribution.
  5. Keep bulky research/provenance files out of the deployable release package while documenting reproducibility paths.

Public deployment decision:

- Best route for users to open a browser and run analysis:
  - Public beta/demo: Hugging Face Spaces using Docker + Streamlit.
  - Serious lab/private-genome use: same app as a local Docker image.
- Rationale:
  - Hugging Face Spaces is designed for hosted ML/demo apps and supports Docker-based apps.
  - Docker gives control over Python version, Streamlit startup, model download/cache paths, and future system dependencies such as BLAST+/MinCED.
  - The 410 MB `models/cas_subtype_extratrees.joblib` artifact should not be committed to ordinary Git; distribute through a Hugging Face model repo, GitHub Release asset, Zenodo/institutional storage, or include it in a controlled Docker image layer.
  - Public hosted demo should include a visible warning not to upload sensitive unpublished genomes.
  - Local Docker should be the recommended path for real lab analysis involving private or large genomes.
- Next session should start with:
  1. Create Docker deployment files (`Dockerfile`, optional `.dockerignore`, startup command).
  2. Add model artifact download/cache strategy.
  3. Add Hugging Face Spaces README/app metadata if deploying there.
  4. Test Docker build/run locally if Docker is available.

Recent production behavior:

- PAM/PFS-unsupported rows are capped at 39.0 when PAM/PFS was evaluated and no hit supports the expected rule.
- Rescored benchmark artifact: `data/curation/benchmark_evaluation_20260520_210025_rescored.csv`.
- PA14/JBD18, PA14/JBD25, and PA14/JBD67 remain high and meet expectations.
- PA14/DMS3 and DGCC7710/Phi2972 are now treated as low-score challenge/control rows because current evidence lacks PAM/PFS support.

Recent diagnostic-only features:

- `pam_subtype_support.csv` challenges repeat-derived subtype calls using observed PAM/PFS flank rules but does not override the repeat model.
- `best_pam_compatibility_score` and `mean_pam_compatibility_score` report probabilistic PAM/PFS compatibility.
- `experimental_pam_weighted_score` was tested but should not replace production scoring because it over-credits ambiguous PA14/DMS3.

Current CCTyper decision:

- CCTyper is not required for SABR runtime and is not a blocker for current development.
- It is useful as an offline Linux/WSL curation and cross-check tool, especially for improving weak Type III labels.
- Native Windows installation was problematic; use WSL Ubuntu.
- Existing repo support includes `envs/cctyper-linux.yml`, `crispr_phage_predictor/ml/check_cctyper_environment.py`, and CCTyper import/collector scaffolding.
- Immediate next practical step is to verify WSL Ubuntu, install Miniforge/CCTyper in WSL, run CCTyper on a small PA14 example, then import `crisprs_near_cas.tab` into SABR training tables.

## Scientific Framing
The tool should not claim that spacer matches prove resistance. A spacer-protospacer match is evidence of possible CRISPR targeting or prior exposure, but biological resistance also depends on:

- CRISPR-Cas system type and subtype
- PAM compatibility
- spacer-protospacer identity
- seed-region integrity
- functional cas genes
- array expression and processing
- phage escape mutations
- anti-CRISPR genes
- genome assembly quality

Output language should use terms such as:

- CRISPR targeting evidence score
- repeat-derived Cas subtype prediction
- PAM/PFS compatibility evidence
- candidate CRISPR-phage immunity interaction
- evidence of CRISPR targeting

## Main User Workflow
1. User uploads one or more bacterial FASTA files.
2. User uploads one or more phage FASTA files.
3. Tool detects CRISPR arrays in bacterial genomes.
4. Tool extracts repeats and spacers from each bacterial genome.
5. Tool compares all bacterial spacers against all phage genomes.
6. Tool predicts or infers CRISPR-Cas type/subtype when possible.
7. Tool selects appropriate candidate PAM rules based on predicted system type.
8. Tool scores every bacterium-phage pair.
9. Tool exports summary and detailed evidence tables.

## Core Output
The main output should be a bacteria-by-phage matrix.

Each matrix cell should contain:

- CRISPR targeting evidence score
- number of matching spacers
- best spacer-protospacer identity
- PAM support level
- predicted CRISPR-Cas type/subtype
- classifier confidence
- evidence interpretation

The GUI should emphasize a heatmap as the first analysis output. Detailed arrays, spacers, and spacer-hit tables should remain available, but they should not dominate the first view.

## Version 0.1 Scope
Initial version should be transparent and reproducible rather than biologically complete.

Planned features:

- Streamlit GUI
- multiple bacterial FASTA uploads
- multiple phage FASTA uploads
- FASTA validation and sequence statistics
- exact-repeat MVP CRISPR array detection
- exact spacer-protospacer matching against uploaded phage genomes
- initial scoring model structure
- project layout suitable for later tests and publications

Version 0.1 should not yet make final resistance claims.

## Initial Algorithm Decision
The first CRISPR detection implementation should be an internal, transparent MVP detector rather than an external dependency. It searches for exact direct repeats of plausible CRISPR repeat length separated by plausible spacer lengths.

Rationale:

- easy to deploy
- no external binary required
- reproducible and inspectable
- useful as a baseline for publication benchmarking

Known limitations:

- does not yet handle approximate or degenerate repeats
- does not infer orientation
- may miss divergent arrays
- may produce false positives in repetitive sequence
- must be benchmarked against established CRISPR callers before scientific claims

Future versions should compare against tools such as CRISPRDetect, CRISPRCasFinder, and PILER-CR where licensing and deployment constraints permit.

## Initial Spacer-Phage Matching Decision
The first matching implementation uses exact spacer-protospacer search on both forward and reverse-complement strands.

Rationale:

- transparent baseline
- fast enough for early uploaded examples
- no external BLAST dependency
- simple to validate in tests

Known limitations:

- does not yet allow mismatches
- does not evaluate seed-region conservation
- does not evaluate PAM flanking sequence
- does not use BLAST, MMseqs2, or indexed search for large datasets

Future matching should add approximate matching and scalable search while keeping exact matches as the baseline evidence layer.

## Machine Learning Direction
A major project goal is to create a curated dataset of bacterial genomes with known CRISPR-Cas type/subtype labels, extract their CRISPR repeat and array features, and train classifiers that can predict system type for new arrays.

The classifier target should be:

- CRISPR-Cas type, such as Type I, II, III, V, VI
- subtype when enough labeled data exists, such as I-E, I-F, II-A, V-A
- associated PAM rule where known

The classifier should use:

- repeat sequence
- repeat length
- repeat k-mer features
- spacer count
- spacer length statistics
- array structure

Runtime feature constraint:

- The deployed SABR model should use only features derivable from plain uploaded FASTA files.
- Trusted annotation, GenBank, literature, CRISPRCasdb, CCTyper, or other tools may be used to create training labels, but those annotations should not be required as model inputs at runtime.
- Organism name, taxonomy, accession, and source database labels should be retained for provenance, auditing, and validation split design, but not used as first-pass model features.
- Nearby Cas gene profiles are biologically powerful but require annotation or gene prediction, so they should be treated as a later optional annotation-aware model, not as the default FASTA-only SABR runtime model.

## Model Strategy
For publication quality, compare several model classes:

1. Similarity baseline
   - nearest known repeat sequence
   - interpretable and important benchmark

2. Classical ML
   - random forest, gradient boosting, SVM, or logistic regression
   - features from repeats, arrays, and cas genes

3. Neural model
   - sequence CNN or transformer-style classifier
   - only useful after enough labeled data is curated

The project should not assume that a neural network is automatically better. Model choice should be justified by benchmark performance, calibration, interpretability, and deployment cost.

## Repeat/Cas-Type Dataset Direction
The next development track is to build a curated CRISPR repeat-to-Cas-type dataset for training a model that can infer likely CRISPR-Cas type/subtype from array sequence features when genome annotation is absent.

Preferred strategy:

1. Build a curated panel of well-characterized bacterial and archaeal genomes or loci with trusted CRISPR-Cas type/subtype annotation from papers, RefSeq/GenBank annotation, CRISPRCasdb/CRISPRCasFinder-style records, or other documented expert sources.
2. Use one consistent repeat/array extraction method, currently MinCED-compatible detection through Diced where practical, to extract consensus repeats, spacer counts, spacer length statistics, and array structure from those genomes.
3. Link each extracted array/repeat to a Cas type/subtype label only when the annotation supports a clear nearby or otherwise documented Cas system association.
4. Use repeat-based predictions as a model output and comparison baseline, not as ground truth labels.
5. Exclude ambiguous, hybrid, orphan, low-confidence, incomplete, or multi-system records from the first training set unless the array-to-Cas association is unambiguous.
6. Treat CCTyper/CRISPRCasTyper as optional future scale-up and cross-check tools, not as a blocking dependency for the next development phase.
7. Later, build a reproducible RefSeq/GenBank complete-genome annotation pipeline using fixed tool versions.

The first local schema is tracked in `crispr_phage_predictor/ml/dataset.py`, with documentation in `data/training/README.md` and technical notes in `handout.md`.

Implemented dataset scaffolding now includes:

- repeat/Cas training table schema
- schema validation
- high-confidence label filtering
- curated-label plus MinCED/Diced/internal detector collector for building repeat/Cas training rows from trusted labeled genomes or loci
- empty curated Cas-type manifest template at `data/training/curated_cas_type_manifest.tsv`
- optional CCTyper `crisprs_near_cas.tab`-style importer into the local schema
- optional CCTyper output-directory collector driven by a manifest CSV
- optional CCTyper environment checker for database and external binaries
- tests for schema validation and CCTyper conversion
- tests for curated-label training collection
- repeat feature extraction for model-ready numeric tables
- baseline random-forest repeat/array classifier for Cas subtype prediction
- nearest-repeat similarity classifier as an interpretable baseline
- command-line train/test evaluator comparing nearest-repeat and random-forest methods

Current external-backend integration status:

- The Streamlit app now exposes selectable CRISPR detection and spacer-matching backends.
- Completed GUI analyses now automatically save timestamped output folders under `outputs/runs/` for follow-up evaluation.
- Saved run artifacts include metadata, bacterial/phage record summaries, candidate arrays, extracted spacers, spacer hits, evidence matrix, and heatmap table.
- CRISPR detection options are currently internal exact-repeat MVP and MinCED.
- Spacer-phage matching options are currently internal exact match and BLASTN.
- BLAST+ is installed on Windows and visible to the app (`blastn` and `makeblastdb` are available on PATH).
- MinCED is not currently installed on Windows PATH.
- The GUI now uses an "Auto recommended" backend mode by default:
  - use MinCED-compatible detection when Diced or `minced` is available
  - fall back to the internal exact-repeat detector when no MinCED-compatible backend is available
  - current Windows-friendly MinCED-compatible backend is the Python package Diced (`diced>=0.1.3`), installed locally in the project venv and Docker image
  - if Diced is absent, the app can still use a `minced` command on PATH
  - use BLASTN for spacer-phage matching when available, otherwise fall back to the internal exact matcher
  - keep manual backend selection as an advanced/debugging option
- The internal detector and exact matcher should remain available as transparent reproducible baselines for benchmarking and troubleshooting rather than being deleted.
- BLASTN runs now expose configurable minimum identity and a full-spacer-alignment toggle.
- BLASTN hit outputs now include identity, alignment length, spacer length, coverage, e-value, and bitscore.
- Evidence matrix outputs now include best identity percent and best coverage percent.
- Streamlit analysis runs now show progress bars for CRISPR detection and spacer-phage matching.
- Internal CRISPR detection progress updates during repeat-length scan steps so long single-record genomes no longer look idle.

Latest full-run output review:

- Latest reviewed run folder: `outputs/runs/20260519_124610`.
- Run configuration:
  - CRISPR detection: MinCED-compatible via Diced (`detection_method=minced`, `detection_backend_detail=diced`)
  - Spacer-phage matching: BLASTN
  - BLAST minimum identity: 90%
  - Require full spacer alignment: true
- Runtime:
  - CRISPR detection: 4.858 seconds
  - BLAST matching: 161.254 seconds
  - Total: 166.141 seconds
- Input scale:
  - 19 bacterial FASTA files
  - 20 bacterial sequences
  - 38 phage FASTA files/sequences
- Output:
  - 36 candidate CRISPR arrays
  - 983 extracted spacers
  - 20 BLAST spacer-protospacer hits
- Positive bacteria-phage evidence rows:
  - `Paeruginosa_PA14.fasta` vs `PhageJBD18.fasta`: 5 unique matching spacers, best identity 100%, coverage 100%
  - `Paeruginosa_PA14.fasta` vs `PhageJBD67.fasta`: 4 unique matching spacers, best identity 100%, coverage 100%
  - `Sislandicus_REY15A.fasta` vs `SIRV2.fasta`: 4 unique matching spacers, best identity 97.44%, coverage 100%
  - `Sthermo_DGCC7710.fasta` vs `Phi2972.fasta`: 2 unique matching spacers, best identity 100%, coverage 100%
  - `Ssolfataricus_P2.fasta` vs `SIRV2.fasta`: 1 unique matching spacer, best identity 90%, coverage 100%
  - `Paeruginosa_PA14.fasta` vs `DMS3.fasta` / `PhageDMS3.fasta`: duplicate-looking DMS3 evidence, one exact spacer each
  - `Paeruginosa_PA14.fasta` vs `PhageJBD25.fasta`: 1 exact spacer

Scientific challenge of latest run:

- The strongest PA14 signals are biologically plausible. Published work reports that `JBD18`, `JBD25`, and `JBD67` are blocked by the native *Pseudomonas aeruginosa* PA14 CRISPR-Cas system, with targets in PA14 CRISPR arrays.
- The tool recovers PA14 spacer evidence for `JBD18`, `JBD25`, `JBD67`, and DMS3-like phages, which is an encouraging sanity check.
- The number of PA14 BLAST hits for `JBD18` and `JBD67` may be higher than the set of experimentally validated functional protospacers, so PAM/seed filtering is needed before interpreting all hits as functional interference.
- The `Sthermo_DGCC7710.fasta` vs `Phi2972.fasta` signal is scientifically important but cautionary: published work indicates spacer matching can occur without interference when PAM context is absent. This makes PAM analysis a high-priority next feature.
- `Sislandicus_REY15A.fasta` and `Ssolfataricus_P2.fasta` vs `SIRV2.fasta` are plausible archaeal CRISPR/SIRV signals but should be framed as candidate historical exposure or targeting evidence unless strain-specific literature or PAM/Cas context confirms resistance.
- `PhageDMS3.fasta` and `DMS3.fasta` appear duplicated by length/accession context and should be deduplicated by accession or sequence hash before benchmark summaries.
- `KpV41.fasta` is only 1,479 bp and looks suspiciously small for a complete phage genome; verify whether it is a fragment, wrong file, or intentionally included record before using it in benchmarking.

Immediate scientific next priorities:

1. Add phage/bacteria input deduplication by sequence hash and accession where available.
2. Add PAM/protospacer-flank extraction around BLAST hits.
3. Add initial PAM compatibility rules for common CRISPR-Cas types where system type is known or inferable.
4. Add seed-region mismatch summaries for BLAST hits.
5. Add a benchmark/validation table comparing tool evidence against known literature labels for PA14/JBD phages, DGCC7710/phi2972, and Sulfolobus/SIRV systems.

Priority correction for the ML track:

- The immediate model-building priority is not phage-host infection outcome data.
- The immediate priority is a well-labeled CRISPR repeat-to-Cas-type/subtype training panel.
- The training panel should contain bacterial/archaeal genomes or CRISPR-Cas loci with clear Cas type/subtype calls, associated CRISPR repeats, spacer statistics, and source provenance.
- Phage-bacteria infection/resistance datasets remain useful later for validating resistance scoring, but they are not the primary input for training the repeat/type classifier.
- The repeat/type classifier should predict likely CRISPR-Cas type/subtype from repeat and array features when genome annotation is absent.
- First labels should come from high-confidence Cas-gene-supported or expert-curated calls, not from phage matching and not from MinCED alone.
- MinCED/Diced should be used to extract repeat and array features, while labels should come from trusted annotation or literature.
- Preferred label sources remain:
  - well-characterized literature examples with explicit CRISPR-Cas type/subtype calls.
  - RefSeq/GenBank complete genomes with reliable Cas gene and CRISPR-Cas subtype annotation.
  - CRISPRCasdb/CRISPRCasFinder-derived records with associated Cas genes and type/subtype annotation.
  - CCTyper/CRISPRCasTyper outputs only as optional future scale-up or label cross-checks where installation and runtime are practical.
- Exclude ambiguous, orphan, hybrid, low-confidence, incomplete, or multi-system records from the first training set unless the array-to-Cas association is unambiguous.

Curated proof-of-concept data track:

- A new curation folder has been added at `data/curation/`.
- `data/curation/proof_of_concept_pairs.tsv` tracks literature-supported bacteria-phage pairs and challenge controls.
- `data/curation/literature_sources.tsv` tracks primary or source URLs used for labels.
- `data/curation/accession_download_manifest.tsv` tracks public accessions to verify/download.
- Current curated pair tiers:
  - Tier A validated positives: PA14 vs JBD18, JBD25, JBD67.
  - Tier B challenge/caution: PA14 vs DMS3, DGCC7710 wild type vs phage 2972, Sulfolobus/SIRV candidate systems.
  - Tier C future acquisition/derived-strain models: DGCC7710 vs phages 858/DT1/2972 BIMs and SMQ-301 non-host context.
  - Tier D QC exclusion: KpV41 until the 1,479 bp local file is verified.
- Current accession manifest status:
  - 10 records already present locally.
  - 1 likely duplicate present locally: DMS3/PhageDMS3.
  - 1 missing high-priority phage: Streptococcus phage DT1 (`AF085222`, also check `NC_002072`).
  - 1 missing host/accession lookup: *Streptococcus thermophilus* SMQ-301.
- These curation files are not final benchmark truth. They are meant to support a proof-of-concept panel with explicit labels, caveats, and citations.

Online data gathered for proof-of-concept expansion:

- Downloaded accession-backed FASTA files into `data/curation/downloads/`.
- New downloaded records:
  - *Streptococcus thermophilus* SMQ-301, `CP011217.1`, 1,861,792 bp, `data/curation/downloads/bacteria/Streptococcus_thermophilus_SMQ-301.fasta`
  - *Streptococcus thermophilus* LMD-9, `NC_008532.1`, 1,856,368 bp, `data/curation/downloads/bacteria/Streptococcus_thermophilus_LMD-9.fasta`
  - Streptococcus phage DT1, `AF085222.2`, 34,815 bp, `data/curation/downloads/phages/Streptococcus_phage_DT1.fasta`
- Download inventory is tracked in `data/curation/downloaded_records.tsv`.
- Added new curation manifest rows:
  - `smq301_2972_downloaded`
  - `smq301_dt1_downloaded`
  - `lmd9_dt1_downloaded`
- Ran a curation-only pilot using downloaded SMQ-301/LMD-9 plus DT1/2972/858:
  - bacteria: SMQ-301 and LMD-9
  - phages: DT1, Phi2972, Phi858
  - detection: MinCED-compatible/Diced
  - matching: BLASTN, 90% identity, full-spacer alignment required
  - arrays: 6
  - spacers: 61
  - hits: 7
  - positive evidence rows:
    - LMD-9 vs DT1: 1 unique spacer, best identity 90%, coverage 100%
    - SMQ-301 vs Phi2972: 2 unique spacers, best identity 96.67%, coverage 100%
    - SMQ-301 vs Phi858: 2 unique spacers, best identity 100%, coverage 100%
    - SMQ-301 vs DT1: 2 unique spacers, best identity 90%, coverage 100%
- Pilot outputs saved under `data/curation/pilot_outputs/`:
  - `strep_curation_evidence_matrix.tsv`
  - `strep_curation_arrays.tsv`
  - `strep_curation_spacer_hits.tsv`

Current CCTyper decision:

- CCTyper installation was difficult on native Windows because the full runtime depends on Linux-oriented bioinformatics tools and databases.
- CCTyper should not be a blocking dependency for SABR runtime development or for the first ML dataset.
- Existing CCTyper import/collection scaffolding can remain in the repository for future Linux/WSL use, but the immediate workflow should not require it.
- The near-term training data workflow is:
  - select well-characterized genomes or loci with trusted Cas type/subtype labels;
  - extract CRISPR repeats and array features consistently with MinCED-compatible detection through Diced where practical;
  - manually or semi-automatically link arrays to the trusted labels only when the association is clear;
  - train and validate repeat/array-feature models using genome, species, or genus holdout splits.
- The local command for the curated workflow is:
  - `python -m crispr_phage_predictor.ml.collect_curated_minced_training data/training/curated_cas_type_manifest.tsv --output data/training/repeats_cas_types.csv`
- Initial curated seed entries have been added to `data/training/curated_cas_type_manifest.tsv`:
  - *Pseudomonas aeruginosa* PA14, `NC_008463.1`, curated Type I-F label, genome-scope label, source keys `Cady2012_PA14_CRISPR;BondyDenomy2013_Acr`
  - *Escherichia coli* K-12 MG1655, `NC_000913.3`, curated Type I-E label, genome-scope label, source key `Ecoli_K12_TypeIE_Review`
- The seed manifest was expanded with additional curated examples:
  - *Streptococcus pyogenes* SF370, `NC_002737.2`, curated Type II-A label from `Spyogenes_SF370_TypeIIA`, restricted to the canonical six-spacer array coordinates `860819-861250`
  - *Francisella tularensis* subsp. *novicida* U112, `NC_008601.1`, curated Type V-A/Cpf1-Cas12a label from `Francisella_U112_Cpf1_TypeVA`
  - *Staphylococcus epidermidis* RP62A, `NC_002976.3`, curated Type III-A label from `Sepidermidis_RP62A_TypeIIIA`
- Additional downloaded but not yet included training candidate:
  - *Leptotrichia shahii* JCM16776, `AP019827`, Type VI-A candidate source `Leptotrichia_shahii_TypeVIA`; excluded from the first seed table because literature indicates multiple CRISPR-Cas systems, so locus-level labels are required before inclusion.
- Running the curated collector on the current seed manifest produced `data/training/repeats_cas_types_seed.csv` with eight rows:
  - two PA14 Type I-F arrays
  - two E. coli K-12 Type I-E arrays
  - one S. pyogenes SF370 Type II-A array
  - two F. novicida U112 Type V-A arrays
  - one S. epidermidis RP62A Type III-A array
- Next seed-expansion targets should add more Type II-A/II-C, Type III, Type V-A, and carefully curated Type VI-A examples. Multi-locus systems such as *Streptococcus thermophilus* and *Leptotrichia shahii* should use contig or array-coordinate label scope rather than broad genome labels.
- A scaled computational candidate importer was added for Vink et al. 2021 / CRISPRCasdb-derived supplementary data:
  - command: `python -m crispr_phage_predictor.ml.import_vink2021_repeats data/training/external_sources/vink_2021_additional_file_2.csv --output data/training/repeats_cas_types_vink2021_candidate.csv --max-per-subtype 200`
  - source table: Vink et al. 2021 Additional file 2, downloaded locally under `data/training/external_sources/`
  - output: `data/training/repeats_cas_types_vink2021_candidate.csv`
  - current output size: 2,161 computational candidate rows
  - filter: repeat subtype must agree with `subtypesinproximity`; rows are collapsed from spacer-level records into accession + repeat + subtype records; at least two unique spacers per accession/repeat/subtype are required; each subtype is capped at 200 rows for balance.
  - label confidence: `computational_proximity`, not manually curated gold truth.
  - current type distribution: Type I 1,132 rows; Type II 422 rows; Type III 530 rows; Type V 32 rows; Type VI 45 rows.
  - current subtype distribution: I-A 154, I-B 200, I-C 200, I-D 113, I-E 200, I-F 200, I-U 58, I-V 7, II-A 200, II-B 22, II-C 200, III-A 200, III-B 200, III-C 16, III-D 114, V-A 28, V-B 4, VI-A 2, VI-B1 41, VI-C 2.
  - random row-level sanity check on the candidate table with `--include-medium-confidence`: nearest-repeat accuracy 0.8281; random-forest accuracy 0.8503.
  - This random split is only a smoke test because repeat-level and genome-level leakage are likely. Publication-quality evaluation still needs genome/species/genus holdout splits and separation between manually curated gold rows and computational candidate rows.
- A larger 1,000-row-per-subtype candidate table was also generated:
  - command: `python -m crispr_phage_predictor.ml.import_vink2021_repeats data/training/external_sources/vink_2021_additional_file_2.csv --output data/training/repeats_cas_types_vink2021_candidate_1k.csv --max-per-subtype 1000`
  - output size: 4,582 computational candidate rows
  - type distribution: Type I 3,166 rows; Type II 776 rows; Type III 563 rows; Type V 32 rows; Type VI 45 rows.
  - subtype distribution: I-A 154, I-B 549, I-C 639, I-D 113, I-E 1,000, I-F 646, I-U 58, I-V 7, II-A 340, II-B 22, II-C 414, III-A 217, III-B 216, III-C 16, III-D 114, V-A 28, V-B 4, VI-A 2, VI-B1 41, VI-C 2.
  - random row-level smoke test with `--include-medium-confidence`: nearest-repeat accuracy 0.9188; random-forest accuracy 0.9180.
  - This larger source is useful for model development, but Type V and Type VI remain underrepresented and need additional targeted curation or additional data sources.
- An uncapped version of the same conservative Vink/CRISPRCasdb-derived import was generated:
  - command: `python -m crispr_phage_predictor.ml.import_vink2021_repeats data/training/external_sources/vink_2021_additional_file_2.csv --output data/training/repeats_cas_types_vink2021_candidate_full.csv --max-per-subtype 0`
  - output size: 4,877 computational candidate rows
  - type distribution: Type I 3,461 rows; Type II 776 rows; Type III 563 rows; Type V 32 rows; Type VI 45 rows.
  - The cap was not the main bottleneck; the conservative proximity-agreement and minimum-spacer filters are. The next meaningful scale-up needs additional sources, not just a higher cap.
- A GenBank annotation-signature importer scaffold was added:
  - path: `crispr_phage_predictor/ml/collect_genbank_annotated_training.py`
  - purpose: query NCBI nuccore by subtype-specific Cas signature terms, download GenBank records, require subtype signature genes in annotations, reject records matching multiple subtype profiles, and create candidate training rows.
  - initial searches show GenBank can be useful for targeted acquisition, especially signature-gene queries such as `cas13a`, `cas12a`, `cas9 csn2`, `cas10 csm`, and `cas10 cmr`.
  - direct large-scale GenBank scanning is currently inefficient because some NCBI hits are assembly project records or very large genomes without annotated CRISPR repeat features; the importer needs batching, record-length filters, and preferably joining GenBank Cas-gene labels to existing repeat databases before it becomes the main scale-up route.
  - Tests cover subtype signature matching, but no GenBank-derived table should be treated as production-quality yet.
- Project documentation/reporting:
  - A living Word report was generated at `docs/SABR_model_development_report.docx`.
  - The report generator is `docs/generate_project_report.py`.
  - Report figures are stored under `docs/figures/`: Cas type distribution, Cas subtype distribution, and the model/validation plan.
  - The report documents the dataset layers, current feature set, planned model ladder, smoke-test results, and validation plan.
  - Regenerate after dataset/model changes with `python docs/generate_project_report.py`.
  - Model comparison metrics are kept in `docs/model_comparison_current.csv`.
  - Required presentation figures to maintain: model accuracy comparison, best-model confusion matrix, per-class F1, ROC curves for probability-capable models, and feature importance for the best tree ensemble.
  - Current best predictor to emphasize unless superseded: ExtraTrees on `data/training/repeats_cas_types_augmented_vink_genbank_targeted.csv` with genus holdout and `--min-class-count 20`.
  - Added model-figure generator: `docs/generate_model_figures.py`.
  - Current model/presentation figures:
    - `docs/figures/model_accuracy_comparison.png`
    - `docs/figures/best_model_confusion_matrix.png`
    - `docs/figures/best_model_per_class_f1.png`
    - `docs/figures/best_model_roc_curve.png`
    - `docs/figures/best_model_feature_importance.png`
    - `docs/figures/best_model_error_by_subtype.png`
    - `docs/figures/best_model_confidence_correct_vs_wrong.png`
    - `docs/figures/best_model_top_errors.png`
  - The best-model held-out prediction export is stored at `docs/best_model_predictions.csv`.
    - Current export size: 1,215 genus-holdout predictions.
    - Correct predictions: 1,112.
    - Wrong predictions: 103.
    - This corresponds to the current best ExtraTrees genus-holdout accuracy of about 0.9152.
  - The Word report now embeds the model comparison, best-predictor, ROC, feature-importance, and error-analysis figures.
- Model evaluation was expanded beyond nearest-repeat and random forest:
  - `crispr_phage_predictor.ml.train_classifier` now supports `nearest_repeat`, `logistic_regression`, `linear_svm`, `gradient_boosting`, and `random_forest`.
  - CLI method selection is available with `--methods`.
  - `--split-strategy row_random` remains a smoke test.
  - `--split-strategy group_holdout --group-column genome_id` performs accession/genome holdout and is the current stronger validation mode.
  - Current model comparison results are documented in `docs/model_comparison_current.csv` and included in `docs/SABR_model_development_report.docx`.
  - On `data/training/repeats_cas_types_vink2021_candidate_full.csv`, random row split accuracies were: nearest-repeat 0.9172, logistic regression 0.8664, linear SVM 0.8475, random forest 0.9148.
  - On the same table with genome/accession holdout, accuracies were: nearest-repeat 0.8942, logistic regression 0.8376, linear SVM 0.8162, random forest 0.8868.
  - Interpretation: nearest-repeat remains a very strong baseline; random forest is close but not clearly better yet. This supports keeping the model ladder honest and prioritizing harder species/genus holdouts before adding neural models.
- FASTA-only feature expansion:
  - The feature extractor now includes additional variables that can be computed from plain FASTA-derived CRISPR arrays: repeat AT percent, GC skew, repeat count, array length estimate, spacer/repeat length ratio, min/max/median/std spacer length fallbacks, repeat terminal base composition, terminal k-mers, repeat self reverse-complement identity, longest terminal inverted stem, and a simple hairpin-like score.
  - These features do not require organism name, taxonomy, accession, GenBank annotation, or Cas-gene annotation.
  - On genome/accession holdout of `repeats_cas_types_vink2021_candidate_full.csv`, random forest improved slightly from 0.8868 to 0.8950 after adding the expanded FASTA-only features, while nearest-repeat remained 0.8942.
  - The current result suggests expanded FASTA-only structure features help modestly, but the nearest-repeat baseline remains extremely competitive.
- Immediate next validation priority:
  - Add species/genus-style holdout validation using metadata only for splitting, not as model input features.
  - Compare row-random, genome/accession holdout, and genus/taxonomy holdout results in the Word report.
  - Use the comparison to decide whether the FASTA-only model generalizes beyond close relatives or mainly memorizes repeat families.
  - If genus holdout performance drops substantially, prioritize more diverse curated data and additional FASTA-only local-context features before adding complex neural models.
- Genus-holdout validation status:
  - The evaluator can now derive `genus` and `species` groups from existing `organism` metadata when those columns are requested for `--group-column`.
  - Metadata are used only for split design, not as model input features.
  - The evaluator also supports `--min-class-count` to remove rare classes before validation; this is needed because rare subtypes such as VI-A can otherwise be held out entirely.
  - Current command: `python -m crispr_phage_predictor.ml.train_classifier data/training/repeats_cas_types_vink2021_candidate_full.csv --include-medium-confidence --test-size 0.25 --split-strategy group_holdout --group-column genus --methods nearest_repeat,random_forest --min-class-count 20`
  - Current genus-holdout result after dropping subtypes with fewer than 20 rows: nearest-repeat accuracy 0.8869; random-forest accuracy 0.9066.
  - Interpretation: random forest now outperforms nearest-repeat under genus holdout, suggesting the expanded FASTA-only features add useful non-nearest-neighbor signal. Type III subtypes remain a major confusion area.
- Stronger controlled model test:
  - Added `extra_trees` to the evaluator as a stronger tree-ensemble model that still works on interpretable FASTA-only tabular features.
  - Current genus-holdout comparison with `--min-class-count 20`: nearest-repeat 0.8869, random forest 0.9066, ExtraTrees 0.9123.
  - ExtraTrees is the best current model under genus holdout, but the gain is modest. This supports trying stronger tabular ensembles, while still prioritizing data diversity and Type III/Type V/Type VI coverage before neural models.
- Hybrid nearest-neighbor feature test:
  - Added `hybrid_extra_trees`, which augments the expanded FASTA-only feature table with nearest-neighbor confidence features computed from the training repeat database.
  - Added neighbor features: best repeat identity, second-best identity, best-vs-second margin, count of neighbors above 90% identity, top subtype vote fraction among close neighbors, and subtype entropy among close neighbors.
  - These are still FASTA-compatible at runtime because they compare the uploaded repeat to the trained reference-repeat database rather than using organism metadata or annotation.
  - Current genus-holdout comparison with `--min-class-count 20`: ExtraTrees 0.9123; hybrid ExtraTrees 0.9131.
  - Interpretation: the hybrid gives only a very small improvement. Keep it as an interpretable confidence-enhanced candidate, but the next large improvement is more likely to come from better Type III/Type V/Type VI data coverage and local FASTA context features than from more model complexity.
- Targeted GenBank+MinCED scale-up:
  - The GenBank annotation collector was tightened to filter NCBI summaries before download/scanning: keep complete genomic records, exclude WGS/project/MAG/scaffold records, and limit to reasonable bacterial sequence lengths.
  - A targeted run for Type V-A and Type VI-A used GenBank annotation signatures as labels and MinCED/Diced extraction as the FASTA-derived feature source.
  - Clean targeted outputs:
    - `data/training/repeats_cas_types_genbank_targeted_va.csv`: 2 Type V-A rows from complete GenBank records.
    - `data/training/repeats_cas_types_genbank_targeted_via.csv`: 5 Type VI-A rows from complete GenBank records.
  - Combined output:
    - `data/training/repeats_cas_types_augmented_vink_genbank_targeted.csv`
    - total rows: 4,884
    - type distribution: Type I 3,461; Type II 776; Type III 563; Type V 34; Type VI 50.
    - subtype increases: V-A increased to 30 rows; VI-A increased to 7 rows.
  - Type III-A/III-B targeted GenBank runs produced no accepted rows because candidate complete genomes often matched multiple subtype profiles, for example I-E plus III-A or III-A plus III-B. These should not be labeled at genome scope; they need locus-level array-to-Cas association before use.
  - Genus-holdout comparison on the augmented table with `--min-class-count 20`: nearest-repeat 0.8922, ExtraTrees 0.9152, hybrid ExtraTrees 0.9119.
  - Interpretation: targeted GenBank+MinCED acquisition slightly improved the best model and helped Type V/VI coverage, but much more targeted acquisition is still needed. Type III remains the main difficult class because many candidate genomes contain multiple CRISPR-Cas systems.
- WSL/CCTyper can be revisited later as a scaling and cross-validation route after the curated proof-of-concept dataset is working.

The first model should use repeat and array features only:

- repeat sequence
- repeat length
- repeat GC percent
- repeat k-mer counts
- spacer count
- mean spacer length
- spacer length variability where available
- min, max, median, and standard deviation of spacer length where available
- repeat terminal motifs and terminal k-mer composition
- repeat self-complementarity and simple hairpin-like scores
- array length, repeat count, and spacer/repeat length ratio where available

Taxonomy and organism metadata should be retained for auditing and split design, but should not be used as first-model features because they may cause shortcut learning.

Validation should hold out whole genomes, species, or genera rather than randomly splitting individual arrays.

## Dataset Plan
Training data should be stored in CSV or TSV form with fields such as:

```csv
genome_id,organism,contig_id,array_start,array_end,repeat_sequence,repeat_length,spacer_count,mean_spacer_length,cas_type,cas_subtype,pam,source
```

Potential data sources to evaluate later:

- CRISPRCasdb
- CRISPRDetect outputs
- CRISPRCasFinder outputs
- curated literature datasets
- RefSeq genomes with known CRISPR-Cas annotations

Each dataset source must be documented for reproducibility.

## Current Local Example Data
The local `data/examples` directory now contains two outcome-labeled groups:

- `data/examples/resistant`
- `data/examples/susceptible`

Current parsing check:

- resistant bacteria: 4 FASTA files, 4 records, 13,904,092 bp total
- resistant phages: 8 FASTA files, 8 records, 280,423 bp total
- susceptible bacteria: 16 FASTA files, 17 records, 75,216,641 bp total
- susceptible phages: 30 FASTA files, 30 records, 2,009,566 bp total

All current files parse as FASTA.

The susceptible folder spelling has been standardized.

This dataset can be used as an early validation set for:

- upload and parsing behavior
- candidate CRISPR array detection
- spacer extraction
- spacer-phage matching
- heatmap behavior
- preliminary separation of resistant versus susceptible examples

Before using it in publication-quality benchmarking, each resistance/susceptibility label should be traceable to a source such as a paper, database entry, or experimental note.

## Deployment Goals
The tool should be easy to deploy and high performing.

Preferred path:

- local Streamlit app for early versions
- command-line pipeline behind the GUI
- Docker support later
- optional external tools such as BLAST+ or MMseqs2 behind clear configuration
- pure-Python fallback where practical

## Publication Goals
To support publication in a bioinformatics journal, the project should eventually include:

- reproducible pipeline
- documented algorithm
- benchmark datasets
- validation experiments
- comparison to existing CRISPR array and CRISPR-Cas typing tools
- calibration analysis for resistance scores
- tests
- versioned example outputs
- citation-ready methods text

## Open Questions
- Which bacterial species or genera should be prioritized first?
- Which phage groups should be prioritized first?
- Should the first working version use BLAST+ or a pure-Python approximate matcher?
- Which CRISPR-Cas types/subtypes are most important for the first classifier?
- Which curated source should be used as the first labeled repeat/type dataset?
- How should unknown or ambiguous Cas types be represented?

## Latest Implementation Status: PAM and Automatic Subtype Bridge

Recent code updates added the first PAM/PFS evaluation layer.

Implemented:

- Spacer hits now carry strand-aware protospacer context:
  - `protospacer_5p_flank`
  - `protospacer_3p_flank`
  - `genomic_upstream_flank`
  - `genomic_downstream_flank`
- Flank extraction is implemented for both internal exact matching and BLASTN matching.
- A new `crispr_phage_predictor/pam.py` module evaluates simple PAM/PFS rules such as:
  - `5prime:CCN`
  - `5prime:AWG`
  - `3prime:NGG`
  - `5prime:TTTN`
- PAM evaluation supports IUPAC ambiguity codes.
- Spacer-hit exports now include:
  - predicted Cas subtype
  - Cas subtype confidence
  - subtype prediction method
  - PAM rule
  - PAM rule source
  - PAM sequence
  - PAM match
  - PAM support level
- Evidence matrix outputs now include:
  - `pam_supported_hits`
  - `pam_evaluated_hits`
  - `pam_support_level`
- A new `crispr_phage_predictor/cas_prediction.py` module provides the first automatic bridge from CRISPR repeat sequence to candidate PAM rule.
- Runtime subtype prediction currently uses the existing nearest-repeat classifier against `data/training/repeats_cas_types_augmented_vink_genbank_targeted.csv`.
- The Streamlit app now defaults to `Auto from predicted subtype` for PAM/PFS evaluation.
- Manual PAM entry remains available only as an expert override, not as the intended default user workflow.

Current conservative curated subtype-to-PAM catalog:

- `I-E -> 5prime:AWG`
- `I-F -> genomic_3prime:GG`
- `I-A -> 5prime:CCN`
- `II-A -> 3prime:NGG`
- `V-A -> 5prime:TTTN`

Important caveats:

- SABR still should not claim confirmed resistance.
- Automatic PAM/PFS support is only applied when subtype confidence is high enough and the predicted subtype has a curated rule.
- Many subtype/PAM combinations remain intentionally unsupported until they are curated from reliable sources.
- Array orientation remains a biological uncertainty for some detections, so PAM interpretations should remain cautious.
- Database-derived PAM strings in the large candidate training table are noisy and should not be blindly applied as runtime rules.

Validation after this implementation:

- Full local test suite passes: 42 tests.
- The Streamlit app was reachable locally at `http://localhost:8501`.

## Latest Implementation Status: Input Deduplication

Recent code updates added sequence-hash deduplication for uploaded FASTA records.

Implemented:

- `FastaRecord` now exposes a normalized SHA-256 `sequence_hash`.
- `crispr_phage_predictor/io.py` now includes:
  - `sequence_hash`
  - `deduplicate_records`
  - `summarize_duplicate_records`
- Streamlit upload diagnostics still show all uploaded records, but analysis runs now use unique records only.
- Duplicate bacterial and phage records are shown in expandable duplicate tables.
- Run summaries now show the number of duplicates excluded.
- Saved run metadata now includes:
  - `bacterial_duplicate_record_count`
  - `phage_duplicate_record_count`
- Saved bacterial/phage record summaries now include `sequence_hash`.

Rationale:

- Duplicate sequences such as `DMS3.fasta` and `PhageDMS3.fasta` should not inflate bacteria-phage matrix cells or spacer-hit counts.
- Deduplication is based on normalized sequence content rather than filename, because filenames may vary for the same underlying accession/sequence.

Validation after this implementation:

- Full local test suite passes: 45 tests.

## Latest Implementation Status: Seed-Region Mismatch Summaries

Recent code updates added a first seed-region mismatch summary layer.

Implemented:

- `SpacerHit` now includes:
  - `seed_region`
  - `seed_length`
  - `seed_mismatches`
  - `seed_mismatch_positions`
- A new `summarize_seed_mismatches` helper compares the spacer and protospacer sequence near the PAM/PFS-proximal end.
- Seed-proximal side is inferred from the evaluated PAM/PFS rule:
  - `5prime:*` rules evaluate the 5-prime end of the spacer/protospacer alignment.
  - `3prime:*` rules evaluate the 3-prime end.
- Seed fields are left unevaluated when no PAM/PFS side is known.
- The evidence matrix now includes:
  - `seed_evaluated_hits`
  - `best_seed_mismatches`
- Streamlit now exposes a `Seed-region length` setting, defaulting to 8 nt.
- Saved run metadata now records the seed length.

Caveats:

- This is a transparent first approximation, not a subtype-specific seed model.
- Different CRISPR-Cas systems can have different seed definitions and mismatch tolerance.
- BLAST hit parsing currently uses ungapped subject sequence for comparison; full gapped alignment reconstruction remains a future improvement for complex indel-containing alignments.

Validation after this implementation:

- Full local test suite passes: 47 tests.

## Latest Implementation Status: Evidence Scoring

Recent code updates replaced the scoring placeholder with a transparent evidence-scoring layer.

Implemented:

- `crispr_phage_predictor/scoring.py` now implements `score_resistance_likelihood`.
- The score is reported as `hypothetical_resistance_score` on a 0-100 scale.
- The score combines:
  - unique matching spacer count
  - best spacer-protospacer identity
  - best spacer alignment coverage
  - PAM/PFS support
  - seed-region mismatch burden
  - Cas subtype prediction confidence, when available
- The evidence matrix now includes:
  - `hypothetical_resistance_score`
  - `current_evidence_level`
  - `evidence_summary`
  - cautious `interpretation`
- Evidence-level labels are phrased as:
  - `strong candidate CRISPR targeting evidence`
  - `moderate candidate CRISPR targeting evidence`
  - `weak candidate CRISPR targeting evidence`
  - `no spacer-match evidence`

Important caveats:

- This is not a calibrated biological resistance probability.
- The score should be interpreted as a transparent ranking of CRISPR targeting evidence.
- Confirmed resistance still requires experimental phenotype or strong literature support.
- Future calibration should compare scores against curated bacteria-phage challenge labels.

Validation after this implementation:

- Full local test suite passes: 50 tests.

## Latest Implementation Status: Accession Metadata and Conflict Checks

Recent code updates added accession extraction and accession conflict diagnostics.

Implemented:

- `FastaRecord` now exposes an `accession` property.
- FASTA summaries now include accession values when they can be extracted from NCBI-style headers.
- Supported accession examples include:
  - `JX495041.1`
  - `AF085222.2`
  - `NC_008463.1`
  - pipe-delimited forms such as `gi|...|ref|NC_008463.1|`
- Saved bacterial/phage record summaries now include both:
  - `accession`
  - `sequence_hash`
- Streamlit diagnostics now report accession conflicts where the same accession appears with more than one distinct sequence hash.

Rationale:

- Sequence hash deduplication handles exact duplicate sequence content.
- Accession metadata improves provenance tracking and benchmark curation.
- Accession conflicts should not be automatically merged, because the same accession paired with different sequence content signals a data-quality issue requiring review.

Validation after this implementation:

- Full local test suite passes: 52 tests.

## Latest Documentation: Development Scheme

A development roadmap has been added at `docs/development_scheme.md`.

The project is now organized into five tracks:

1. Core analysis pipeline
   - FASTA parsing
   - sequence hash/accession diagnostics
   - CRISPR array detection
   - spacer extraction
   - spacer-phage matching
   - protospacer flank extraction
   - subtype prediction
   - PAM/PFS evaluation
   - seed mismatch summary
   - hypothetical targeting score
   - exports and heatmap
2. Data expansion
   - grow repeat/Cas subtype rows
   - grow curated bacteria-phage benchmark pairs
   - prioritize Type III, Type V, and Type VI gaps
3. Model development
   - nearest-repeat baseline
   - tabular FASTA-only models
   - genome/species/genus holdout validation
   - confidence calibration
4. Benchmarking and calibration
   - formal bacteria-phage labels
   - positives, negatives, host-range controls, PAM-failure controls, and anti-CRISPR cases
   - score calibration once labels are large enough
5. Productization
   - Streamlit GUI
   - CLI workflow
   - saved artifacts
   - documentation
   - packaging

## Latest Implementation Status: Benchmark Label Schema

Recent code updates added a formal benchmark-label schema for curated bacteria-phage validation pairs.

Implemented:

- New module: `crispr_phage_predictor/benchmark.py`
- New curation table: `data/curation/benchmark_labels.tsv`
- New tests: `tests/test_benchmark.py`
- Updated curation documentation: `data/curation/README.md`

Benchmark labels now separate:

- observed phenotype label
- CRISPR-mediated resistance label
- CRISPR evidence level
- PAM evidence level
- anti-CRISPR status
- host-range status
- expected SABR behavior
- curation confidence
- source keys and notes

Current benchmark-label table size:

- 8 rows total
- 3 validated rows
- 2 challenge rows
- 2 rows needing literature review
- 1 exclusion/QC row

Important details:

- The benchmark table is stricter than `proof_of_concept_pairs.tsv`.
- Rows without source keys are not allowed as valid curated labels.
- Archaeal SIRV rows are currently marked `needs_literature_review` until stronger source attribution is added.
- `KpV41` remains `exclude_until_verified` because the local file is suspiciously short.

Validation after this implementation:

- Benchmark TSV validates through `load_benchmark_label_table`.
- Full local test suite passes: 56 tests.

## Latest Implementation Status: Benchmark Run Evaluation

Recent code updates added tooling to join saved SABR run outputs to curated benchmark labels.

Implemented:

- `crispr_phage_predictor.benchmark.evaluate_benchmark_run`
- `crispr_phage_predictor.benchmark.summarize_benchmark_evaluation`
- CLI entrypoint:
  - `python -m crispr_phage_predictor.evaluate_benchmark <run_dir> --output <joined_eval.csv>`
- The evaluator joins `evidence_matrix.csv` to `data/curation/benchmark_labels.tsv` using local bacteria/phage file names.
- It includes accession values from saved `bacterial_records.csv` and `phage_records.csv` when available.
- It supports both old and new run schemas:
  - new runs use `hypothetical_resistance_score`;
  - old runs without that column use a legacy unique-spacer proxy for expectation checks.
- It reports `score_expectation_result` values such as:
  - `meets_expectation`
  - `below_expected`
  - `above_expected`
  - `not_evaluated`
  - `excluded`

First joined evaluation artifact:

- `data/curation/benchmark_evaluation_20260519_124610.csv`

Evaluation of historical run `outputs/runs/20260519_124610`:

- benchmark rows evaluated: 8
- matched rows: 7
- expectation results:
  - `meets_expectation`: 6
  - `below_expected`: 1
  - `above_expected`: 1
  - `not_evaluated`: 1

Important interpretation:

- This historical run predates automatic PAM, seed, subtype, deduplication, and numeric scoring.
- The initial joined evaluation is useful for checking benchmark wiring, but it should not be treated as final model/scoring performance.
- The next meaningful benchmark evaluation should be run after generating a fresh SABR output folder with the current pipeline.

Validation after this implementation:

- Full local test suite passes: 58 tests.

## Latest Fresh Benchmark Run

A fresh full example-data run was generated with the current pipeline after adding:

- sequence-hash deduplication
- accession metadata
- automatic subtype-to-PAM/PFS rule selection
- PAM/PFS evaluation
- seed mismatch summaries
- hypothetical resistance scoring

Corrected implementation detail:

- The initial curated Type I-F rule was too conservative/misoriented for the PA14 DNA-targeting benchmark.
- The current catalog now uses:
  - `I-F -> genomic_3prime:GG`
  - `I-A -> 5prime:CCN`
- Unsupported PAM/PFS evidence is now penalized more strongly in scoring, and seed mismatch bonuses are applied only when at least one PAM/PFS-supported hit exists.

Fresh run:

- Run folder: `outputs/runs/20260520_202705`
- Detection: MinCED-compatible via Diced
- Matching: BLASTN
- BLAST minimum identity: 90%
- Full-spacer alignment required: true
- PAM mode: automatic from predicted subtype
- Seed length: 8
- Runtime:
  - CRISPR detection: 1.739 seconds
  - BLAST matching: 71.805 seconds
  - total: 74.019 seconds
- Input after deduplication:
  - bacteria: 21 raw records, 19 unique records
  - phages: 38 raw records, 37 unique records
  - duplicate bacterial records excluded: 2
  - duplicate phage records excluded: 1
- Output:
  - 34 candidate CRISPR arrays
  - 965 extracted spacers
  - 19 spacer-phage hits

Benchmark evaluation:

- Joined evaluation artifact: `data/curation/benchmark_evaluation_20260520_202705.csv`
- Full test suite after these changes: 59 tests passing.

Key benchmark rows:

- `pa14_jbd18`: score 100.0, PAM compatible, meets high-score expectation.
- `pa14_jbd67`: score 100.0, PAM compatible, meets high-score expectation.
- `pa14_jbd25`: score 28.0, PAM not supported, below high-score expectation.
  - Needs targeted review because literature treats this as a validated PA14 CRISPR resistance positive.
  - Possible causes include array orientation, protospacer flank orientation, exact hit choice, or curated PAM rule details.
- `dgcc7710_2972_wt`: score 38.53, PAM not supported, meets low-score expectation.
  - This is the desired behavior for the PAM-failure control.
- `pa14_dms3`: score 28.0, below moderate-score expectation.
  - This may be acceptable depending on how DMS3 is ultimately labeled; it remains a challenge/control case.
- `rey15a_sirv2`: score 99.36, above moderate-score expectation.
  - This row is still `needs_literature_review`, so the high score should trigger deeper archaeal SIRV curation rather than be treated as validated resistance.
- `ssolfataricus_p2_sirv2`: score 25.5, below moderate-score expectation.
  - Also requires source-backed archaeal curation.

Immediate scientific/debugging priority from the fresh run:

1. Investigate the PA14/JBD25 hit and PAM context.
2. Confirm Type I-F PAM orientation/motif against the primary PA14/JBD literature and structural/review sources.
3. Review whether SABR should evaluate both canonical and shifted/slipped PAM windows for Type I-F.
4. Add benchmark notes for PA14/JBD25 once the exact validated protospacer/PAM context is reconciled.
5. Revisit archaeal Type I-A/SIRV labels and PAM rules before treating high SIRV scores as validated positives.

## Latest Implementation Status: Flexible BLAST Coverage Threshold

Recent code updates added configurable BLAST spacer coverage filtering.

Implemented:

- BLAST matching now accepts `min_coverage` in addition to `min_identity`.
- `find_spacer_hits_for_records` exposes `blast_min_coverage`.
- Streamlit now exposes:
  - `BLASTN minimum identity`
  - `BLASTN minimum spacer coverage`
  - `Require full spacer alignment`
- Default GUI behavior is now flexible high-coverage mode:
  - minimum identity: 90%
  - minimum spacer coverage: 95%
  - full-spacer alignment not required by default
- Strict full-length mode remains available with the checkbox.
- Saved run metadata now records `blast_min_coverage`.

Validation:

- Added parser tests for minimum coverage filtering.
- Full local test suite passes: 61 tests.

Flexible benchmark run:

- Run folder: `outputs/runs/20260520_204026`
- BLAST settings:
  - minimum identity: 90%
  - minimum spacer coverage: 95%
  - full-spacer alignment required: false
- Output:
  - 34 arrays
  - 965 spacers
  - 22 spacer-phage hits
- Joined benchmark artifact:
  - `data/curation/benchmark_evaluation_20260520_204026.csv`

Important finding:

- Flexible coverage increased total hits from 19 to 22, but it did not fix the PA14/JBD25 benchmark miss.
- PA14/JBD25 already has a full-length exact 32/32 hit at 100% identity and 100% coverage.
- The current low score is caused by unsupported PAM context:
  - detected JBD25 hit: `Paeruginosa_PA14.fasta|NC_008463.1|array_1|spacer_1`
  - JBD25 coordinates: 30465-30496 on strand `-`
  - current PAM rule before the genomic-coordinate fix: `3prime:GG`
  - extracted PAM sequence: `GC`
  - PAM support: `not_supported`
- Therefore, the next PA14/JBD25 issue is not full-length stringency; it is PAM orientation/window/literature-target reconciliation.

Next specific PA14/JBD25 debug tasks:

1. Extract both upstream/downstream and shifted PAM windows around the JBD25 protospacer.
2. Compare those windows against the PA14/JBD25 validated target in Cady et al. 2012.
3. Determine whether Type I-F PAM checking needs shifted/slipped windows or array-orientation correction.
4. Add a regression test once the correct validated PA14/JBD25 PAM context is resolved.

## Runtime Cas Subtype Prediction Clarification

Current runtime behavior:

- SABR now predicts array Cas subtype at runtime with the saved ExtraTrees artifact when `models/cas_subtype_extratrees.joblib` is present.
- If the artifact is missing or cannot be loaded, SABR falls back to the nearest-repeat classifier.
- It then maps the predicted subtype to a small curated subtype-to-PAM/PFS catalog.
- This means the current preferred runtime flow is:
  - CRISPR repeat and array features -> ExtraTrees subtype prediction -> curated PAM/PFS rule -> flank check.

Important correction:

- The best validated model in the project reports is not nearest-repeat.
- Current best documented model is ExtraTrees on `data/training/repeats_cas_types_augmented_vink_genbank_targeted.csv` with genus holdout and `--min-class-count 20`.
- Current best documented accuracy is about 0.9152 under genus holdout.
- Nearest-repeat remains an important interpretable baseline, but it is not the current best model.

Implemented model-artifact work:

- New module: `crispr_phage_predictor/ml/model_artifact.py`
- New export CLI: `python -m crispr_phage_predictor.ml.export_model --output models/cas_subtype_extratrees.joblib`
- Artifact written locally to `models/cas_subtype_extratrees.joblib`.
- Artifact export summary:
  - training rows: 4,848
  - classes: `I-A`, `I-B`, `I-C`, `I-D`, `I-E`, `I-F`, `I-U`, `II-A`, `II-B`, `II-C`, `III-A`, `III-B`, `III-D`, `V-A`, `VI-B1`
- `.gitignore` excludes `models/*.joblib` because the local artifact is large.
- `cas_prediction.py` now reports `prediction_method` as `extra_trees` when artifact inference is used and `nearest_repeat` when fallback inference is used.

Artifact-backed runtime sanity check:

- PA14 arrays predict `I-F` with confidence `1.000`, method `extra_trees`, PAM rule `genomic_3prime:GG`.
- *S. thermophilus* arrays include `II-C`, `III-A`, `I-E`, and `II-A` calls.

Validation after this implementation:

- Full local test suite passes: 63 tests after adding explicit genome-coordinate PAM rule support.

Latest artifact-backed full run:

- Previous artifact-backed run folder: `outputs/runs/20260520_205322`
- Detection: MinCED-compatible via Diced
- Matching: BLASTN
- BLAST minimum identity: 90%
- BLAST minimum spacer coverage: 95%
- Full-spacer alignment required: false
- PAM mode: automatic from predicted subtype
- Seed length: 8
- Runtime:
  - CRISPR detection: 1.613 seconds
  - BLAST matching: 65.103 seconds
  - total: 67.070 seconds
- Input after deduplication:
  - bacteria: 19 unique records, with 2 duplicate bacterial records excluded
  - phages: 37 unique records, with 1 duplicate phage record excluded
- Output:
  - 34 candidate CRISPR arrays
  - 965 extracted spacers
  - 22 spacer-phage hits
- Joined benchmark artifact:
  - `data/curation/benchmark_evaluation_20260520_205322.csv`

Important benchmark status from artifact-backed run:

- Validated PA14/JBD18 and PA14/JBD67 remain strong.
- PA14/JBD25 was still below expectation in run `20260520_205322` because the full-length exact hit failed the previous strand-oriented `I-F -> 3prime:GG` PAM check.
- Flexible BLAST thresholds did not solve JBD25 because the problem was PAM window/orientation/literature-target reconciliation, not spacer length or identity.
- Debugging the JBD25 window showed the raw phage-genome downstream flank starts with `GG`, while the strand-oriented 3-prime flank starts with `GC`.
- The Type I-F catalog has therefore been updated to `genomic_3prime:GG`; this needs a fresh run and regression test to confirm PA14/JBD25 behavior without breaking other benchmark controls.

Fresh run after genome-coordinate Type I-F PAM fix:

- Run folder: `outputs/runs/20260520_210025`
- Joined benchmark artifact:
  - `data/curation/benchmark_evaluation_20260520_210025.csv`
- Benchmark expectation counts:
  - `meets_expectation`: 3
  - `below_expected`: 2
  - `above_expected`: 2
  - `not_evaluated`: 1
- Validated PA14 rows:
  - `pa14_jbd18`: score 100.0, PAM compatible, meets high-score expectation.
  - `pa14_jbd25`: score 83.0, PAM compatible, meets high-score expectation.
  - `pa14_jbd67`: score 100.0, PAM compatible, meets high-score expectation.
- Remaining control/calibration issue:
  - `dgcc7710_2972_wt`: score 45.0, PAM not supported, above low-score expectation because flexible matching finds 3 spacers. This row should drive the next scoring-calibration pass, not the PAM-orientation fix.

## Current Data Point Counts

Current local counts as of the latest review:

- Local example FASTA files: 58
- Local example FASTA records: 59
  - 21 bacterial records
  - 38 phage records
- Curated proof-of-concept bacteria-phage pairs: 15
- Literature/source tracking rows: 14
- Accession manifest rows: 14
- Downloaded accession-backed records: 7
- Latest full SABR run matrix cells: 703
- Latest full SABR run detected CRISPR arrays: 34
- Latest full SABR run extracted spacers: 965
- Latest full SABR run spacer-phage hits: 22
- Current main ML repeat/Cas training table rows: 4,884
- High-confidence curated seed ML rows: 8
- Current saved runtime ExtraTrees artifact training rows: 4,848
- Current benchmark-label rows: 8

Interpretation:

- The Cas subtype ML table is useful for development, but it is mostly computational candidate data rather than final gold-standard truth.
- The curated bacteria-phage benchmark panel is still too small for publication-quality resistance scoring validation.
- The main data scale-up priority is to increase both:
  - repeat/Cas subtype training rows, especially Type III, Type V, and Type VI systems;
  - curated bacteria-phage benchmark pairs, including positives, negatives, anti-CRISPR examples, and PAM-failure controls.

## Near-Term Execution Plan

Immediate engineering priorities:

1. Calibrate the hypothetical resistance score around PAM-unsupported multi-spacer controls, starting with `dgcc7710_2972_wt`.
2. Expand curated bacteria-phage benchmark pairs toward an initial 100-300 pair panel.
3. Expand data acquisition for repeat/Cas subtype training, targeting 20,000-50,000 rows over time.
4. Add model-artifact metadata reporting to exported SABR runs so each run records the model version/hash used for subtype prediction.
5. Add stricter accession-aware review workflows for conflicts and version mismatches.
6. Revisit archaeal SIRV labels and Type I-A/SIRV PAM rules before treating high SIRV scores as validated positives.

Data scale-up strategy:

- Prioritize CRISPRCasdb/CRISPRCasFinder-style records with array-to-Cas associations.
- Use locus-level labels rather than broad whole-genome labels when multiple CRISPR-Cas systems exist.
- Target underrepresented systems using subtype-specific signatures such as `cas12a`, `cas12b`, `cas13a`, `cas13b`, `cas10 csm`, and `cas10 cmr`.
- Use CCTyper/CRISPRCasTyper later in WSL/Linux as a scale-up and cross-check route, but do not make it a blocking SABR runtime dependency.
- Keep all source, accession, confidence, and caveat metadata attached to every curated row.

## Current Terminology Decision: Targeting Evidence First

SABR should now be framed as a CRISPR spacer-targeting evidence and repeat-based
Cas-typing tool, not as a direct resistance caller.

Current preferred wording:

- CRISPR targeting evidence
- spacer-protospacer targeting evidence
- candidate CRISPR-phage interaction
- repeat-derived Cas type/subtype prediction
- PAM/PFS compatibility evidence
- CRISPR targeting evidence score

Avoid presenting the primary output as a resistance score. Biological resistance
can be caused or modified by many factors beyond spacer matching, including Cas
gene functionality, expression, phage escape, anti-CRISPR genes, host range, and
experimental conditions. Resistance should remain a later interpretation layer
only when phenotype, literature, or benchmark evidence supports it.

Implementation note:

- New code should prefer `crispr_targeting_score`,
  `score_crispr_targeting_evidence`, and
  `build_crispr_targeting_evidence_matrix`.
- Backward-compatible aliases and the historical
  `hypothetical_resistance_score` column remain for older benchmark artifacts
  and tests during the transition.

## Latest Implementation Status: PAM-Unsupported Score Calibration

Recent code updates calibrated the CRISPR targeting evidence score for
PAM/PFS-unsupported controls.

Implemented:

- When PAM/PFS was evaluated and no spacer hit supports the expected rule, the
  targeting score is capped at 39.0.
- This prevents multiple exact spacer matches from being elevated to moderate
  evidence when the required PAM/PFS context is absent.
- PAM/PFS-supported hits are not capped by this rule.
- Added a regression test for a three-spacer, high-identity,
  PAM-unsupported case.

Rationale:

- The `dgcc7710_2972_wt` benchmark control demonstrates that spacer matches can
  occur without interference when PAM context is absent.
- The score should therefore keep PAM-unsupported multi-spacer cases in the weak
  targeting-evidence range unless additional evidence supports targeting.

Validation after this implementation:

- Full local test suite passes: 65 tests.

Rescored benchmark check:

- Recomputed the evidence matrix for historical run `outputs/runs/20260520_210025`
  without rerunning BLAST, using the current scoring code.
- Rescored output folder:
  `outputs/runs/20260520_210025_rescored`
- Joined benchmark artifact:
  `data/curation/benchmark_evaluation_20260520_210025_rescored.csv`
- Key results:
  - `pa14_jbd18`: score 100.0, meets high-score expectation.
  - `pa14_jbd25`: score 83.0, meets high-score expectation.
  - `pa14_jbd67`: score 100.0, meets high-score expectation.
  - `dgcc7710_2972_wt`: score 39.0, PAM not supported, now meets
    low-score expectation.
- Remaining benchmark interpretation issues:
  - `pa14_dms3` remains below moderate expectation and should stay a
    challenge/control row until anti-CRISPR and duplicate-file context are
    clarified.
  - `rey15a_sirv2` remains above moderate expectation, but the row still needs
    stronger archaeal/SIRV literature curation before being treated as a
    validated benchmark.
  - `ssolfataricus_p2_sirv2` remains below moderate expectation and also needs
    source-backed curation.

## Latest Implementation Status: Model Artifact Metadata in Run Outputs

Recent code updates added runtime Cas subtype model metadata to saved SABR run
metadata.

Implemented:

- `crispr_phage_predictor/ml/model_artifact.py` now exposes
  `model_artifact_metadata`.
- The metadata includes:
  - artifact path
  - whether the artifact exists
  - artifact SHA256 hash
  - model method
  - training table
  - training row count
  - min-class count
  - random state
  - estimator count
  - class labels
  - load error, if artifact loading fails
- `save_analysis_run` now writes this under `cas_model_artifact` in
  `run_metadata.json`.
- The Streamlit app passes current artifact metadata into saved runs.
- Missing artifacts are handled gracefully with `artifact_exists=false`.

Current local artifact metadata check:

- Artifact path: `models/cas_subtype_extratrees.joblib`
- Artifact exists: true
- SHA256 prefix: `0021905ca2eb2c4d`
- Method: `extra_trees`
- Training rows: 4,848
- First classes: `I-A`, `I-B`, `I-C`, `I-D`, `I-E`

Validation after this implementation:

- Full local test suite passes: 67 tests.

## Latest Curation Update: PA14/DMS3 Benchmark Row

Recent curation updates changed the `pa14_dms3` benchmark row from a moderate
score expectation to a low-score challenge expectation.

Rationale:

- The current local DMS3 evidence is a single spacer hit with PAM/PFS not
  supported.
- The local files combine duplicate-looking DMS3 records/accessions
  (`PhageDMS3.fasta` and `DMS3.fasta`).
- DMS3 appears in PA14 CRISPR and anti-CRISPR/lysogen literature, so it should
  not be treated as a simple validated CRISPR targeting positive without
  splitting the exact accession/construct and anti-CRISPR context.

Implemented:

- `data/curation/benchmark_labels.tsv`
  - `pa14_dms3` now uses `expected_sabr_behavior=low_score_expected`.
  - `pam_evidence_level` is set to `unknown`.
  - Notes now explicitly require splitting DMS3/PhageDMS3 and curating exact
    protospacer, PAM, and anti-CRISPR context before stronger labeling.
- `data/curation/proof_of_concept_pairs.tsv`
  - DMS3 notes now state that the current SABR hit lacks PAM support and that
    duplicate files should be split before stronger benchmark use.

Updated rescored benchmark status:

- `pa14_dms3`: score 28.0, PAM not supported, now meets low-score challenge
  expectation.
- `dgcc7710_2972_wt`: score 39.0, PAM not supported, meets low-score
  expectation.
- Remaining curation problems are now concentrated in the archaeal SIRV rows:
  `rey15a_sirv2` remains above expectation and `ssolfataricus_p2_sirv2` remains
  below expectation.

Validation:

- `tests/test_benchmark.py` passes.
- Updated joined benchmark artifact:
  `data/curation/benchmark_evaluation_20260520_210025_rescored.csv`.

## Latest Data Expansion Experiment: Type V/VI GenBank Candidates

Recent data-expansion work targeted underrepresented Cas subtype training rows,
especially Type V and Type VI.

Generated artifacts:

- Download/cache directory:
  `data/training/genbank_sources_expansion/`
- Candidate expansion table:
  `data/training/repeats_cas_types_genbank_expansion_iii_v_vi.csv`
- Expanded combined training table:
  `data/training/repeats_cas_types_augmented_vink_genbank_targeted_expanded.csv`
- Count summary:
  `docs/cas_subtype_expansion_summary.csv`
- Updated model comparison table:
  `docs/model_comparison_current.csv`

Collection details:

- Used the existing GenBank signature-candidate collector with Diced/MinCED
  array extraction.
- Target profiles: `III-A`, `III-B`, `V-A`, `VI-A`.
- The completed pass produced 69 candidate rows.
- After collapsing duplicate `(genome_id, repeat_sequence, cas_subtype)` keys
  and removing rows already present in the current augmented table, 47 unique
  new rows were added.

Subtype count changes:

- `V-A`: 30 -> 37 rows, delta +7.
- `VI-A`: 7 -> 47 rows, delta +40.
- `III-A`: unchanged at 217 rows.
- `III-B`: unchanged at 216 rows.
- `III-D`: unchanged at 114 rows.
- `VI-B1`: unchanged at 41 rows.
- `VI-C`: unchanged at 2 rows.

Genus-holdout evaluation on the expanded table:

- Rows used after `--min-class-count 20`: 4,902.
- Train rows: 3,673.
- Test rows: 1,229.
- Nearest-repeat accuracy: 0.8950.
- ExtraTrees accuracy: 0.9032.
- Hybrid ExtraTrees accuracy: 0.9040.

Interpretation:

- This expansion improves V-A and especially VI-A row coverage, but it does not
  improve the current overall best model.
- The current best documented runtime candidate remains ExtraTrees on
  `data/training/repeats_cas_types_augmented_vink_genbank_targeted.csv`, with
  genus-holdout accuracy about 0.9152.
- Do not replace `models/cas_subtype_extratrees.joblib` with a model trained on
  the expanded table yet.
- The Type III weakness remains unresolved; the GenBank signature pass did not
  add new unique Type III rows. Type III likely needs more targeted
  locus-level curation rather than broad whole-record signature queries.

## Latest Exploratory Feature: PAM/PFS Subtype Support Diagnostic

Recent code updates added a diagnostic-only analysis that tests whether observed
protospacer flanks support curated PAM/PFS rules for known subtypes.

Purpose:

- Challenge the repeat-based Cas subtype prediction using independent
  spacer-hit flank evidence.
- Report agreement, conflict, or lack of PAM/PFS support.
- Do not change the primary Cas subtype call.
- Do not change targeting scores.

Implemented:

- New pipeline function: `summarize_pam_subtype_support`.
- Saved run output: `pam_subtype_support.csv`.
- Streamlit expander: `Exploratory PAM/PFS subtype support`.
- Tests cover agreement and conflict cases.

Current rule catalog tested by this diagnostic:

- `I-E -> 5prime:AWG`
- `I-F -> genomic_3prime:GG`
- `I-A -> 5prime:CCN`
- `II-A -> 3prime:NGG`
- `V-A -> 5prime:TTTN`

Diagnostic result on `outputs/runs/20260520_210025_rescored`:

- `Paeruginosa_PA14` array 1:
  - repeat prediction `I-F`
  - top PAM-supported subtype `I-F`
  - agreement: `agrees_top_subtype`
- `Paeruginosa_PA14` array 2:
  - repeat prediction `I-F`
  - top PAM-supported subtype `II-A`
  - support counts: `I-F:3;II-A:4`
  - agreement: `conflicts_with_top_subtype`
  - interpretation caveat: this is likely affected by broad/generic `NGG`
    matching and should not override repeat-based I-F prediction.
- `Sislandicus_REY15A` array 3, `Ssolfataricus_P2` array 10, and
  `Sthermo_DGCC7710` array 12:
  - repeat-predicted subtype was not supported by current curated PAM/PFS rules.
- `Sthermo_DGCC7710` array 15:
  - no PAM subtype support detected.

Interpretation:

- The diagnostic is useful for highlighting repeat/PAM agreement and possible
  conflicts.
- Generic PAM motifs can create misleading subtype support, especially when a
  broad rule such as `II-A -> NGG` overlaps with narrower DNA-targeting motifs.
- This supports using PAM/PFS subtype support as an audit layer, not as a
  replacement for repeat-based Cas subtype prediction.

Validation:

- Full local test suite passes: 69 tests.

## Latest Modeling Experiment: Hierarchical ExtraTrees

Recent code updates added an experimental hierarchical model evaluator.

Model idea:

- First predict broad Cas type, such as Type I, Type II, Type III, Type V, or
  Type VI.
- Then predict subtype inside the predicted broad type.
- Evaluate final subtype accuracy against the same holdout splits used for the
  flat models.

Implemented:

- New evaluator method: `hierarchical_extra_trees`.
- The method trains:
  - one global ExtraTrees model for `cas_type`;
  - one per-type ExtraTrees model for `cas_subtype` when a type has multiple
    subtypes;
  - constant subtype fallback for broad types with one subtype in training;
  - global subtype fallback only if a predicted type has no subtype model.
- Added unit coverage in `tests/test_ml_train_classifier.py`.

Genus-holdout result on the current best dataset:

- Dataset: `data/training/repeats_cas_types_augmented_vink_genbank_targeted.csv`
- Rows used after `--min-class-count 20`: 4,848.
- Train rows: 3,633.
- Test rows: 1,215.
- Flat ExtraTrees accuracy: 0.9152.
- Hierarchical ExtraTrees accuracy: 0.9012.

Interpretation:

- The first hierarchical model underperformed the flat ExtraTrees model.
- It did not solve the Type III weakness:
  - `III-A` recall dropped from 0.43 in flat ExtraTrees to 0.30.
  - `III-B` recall dropped from 0.57 to 0.48.
  - `III-D` recall dropped from 0.40 to 0.27.
- Do not replace the runtime model with this hierarchy.
- Future hierarchical work should focus on uncertainty/abstention or
  type-level reporting, not assume a naive hard type-then-subtype split will
  improve accuracy.

Validation:

- `tests/test_ml_train_classifier.py` passes after adding the method.

## Latest Data Expansion Experiment: Type III-A Rescue Rows

Recent data-expansion work attempted to add more Type III training rows.

What was tried:

- Re-ran the GenBank signature-candidate collector with a larger Type III query:
  `III-A,III-B`, max 60 records per profile.
- Fixed the collector to batch NCBI ESummary requests after the larger query hit
  `HTTP 414: Request-URI Too Long`.
- The strict collector wrote zero rows because most downloaded Type III
  candidates matched both `III-A` and `III-B` signatures and were correctly
  rejected as ambiguous.
- A manual diagnostic found 5 unambiguous `III-A` Corynebacterium records with
  Diced/MinCED-detected arrays.

Generated artifacts:

- GenBank cache:
  `data/training/genbank_sources_typeiii_expansion/`
- Strict Type III expansion output:
  `data/training/repeats_cas_types_genbank_typeiii_expansion.csv`
  - 0 rows, because strict filtering rejected ambiguous Type III records.
- Type III rescue candidate rows:
  `data/training/repeats_cas_types_genbank_typeiii_rescue.csv`
  - 13 rows.
- Combined expanded table:
  `data/training/repeats_cas_types_augmented_vink_genbank_targeted_expanded_typeiii.csv`
  - 4,944 rows.
- Count summary:
  `docs/cas_subtype_typeiii_expansion_summary.csv`

Added rows relative to current best dataset:

- `III-A`: 217 -> 230 rows, delta +13.
- `V-A`: 30 -> 37 rows, delta +7.
- `VI-A`: 7 -> 47 rows, delta +40.
- `III-B`, `III-C`, `III-D`, `VI-B1`, and `VI-C`: unchanged.

Genus-holdout evaluation on expanded Type III table:

- Rows used after `--min-class-count 20`: 4,915.
- Train rows: 3,686.
- Test rows: 1,229.
- Nearest-repeat accuracy: 0.8950.
- ExtraTrees accuracy: 0.9032.
- Hybrid ExtraTrees accuracy: 0.9032.

Interpretation:

- The rescue rows modestly improved Type III-A recall in the expanded
  experiment compared with the prior V/VI-only expansion, but overall
  performance remains below the current best flat ExtraTrees model at 0.9152.
- Do not replace the runtime model.
- The main blocker is still high-quality, locus-level Type III-B/III-D curation.
- Broad GenBank signatures are not enough for Type III because many records
  contain annotations that satisfy both `csm` and `cmr` profiles or otherwise
  imply multiple/ambiguous systems.

Validation:

- `tests/test_ml_collect_genbank_annotated_training.py` now covers ESummary
  batching.

## Latest Exploratory Feature: Probabilistic PAM/PFS Compatibility

Recent code updates added a diagnostic PAM/PFS compatibility score.

Implemented:

- `PamEvaluation` now includes `compatibility_score`.
- `SpacerHit` now includes `pam_compatibility_score`.
- Detailed spacer-hit output now includes `pam_compatibility_score`.
- Evidence matrix output now includes:
  - `best_pam_compatibility_score`
  - `mean_pam_compatibility_score`
- The current targeting score logic is unchanged.

Current scoring heuristic:

- Exact motif compatibility scores 1.0.
- Incompatible motif positions score 0.0.
- IUPAC-compatible ambiguous positions count as compatible.
- The score is the mean per-position compatibility across the evaluated motif.
- Missing, invalid, or insufficient flank evaluations keep score as null.

Diagnostic check on `outputs/runs/20260520_210025_rescored`:

- `PA14/JBD18`: best 1.0, mean 0.6.
- `PA14/JBD25`: best 1.0, mean 1.0.
- `PA14/JBD67`: best 1.0, mean 0.5.
- `PA14/DMS3`: best 0.5, mean 0.5.
- `DGCC7710/Phi2972`: best 0.333333, mean 0.333333.
- `Ssolfataricus_P2/SIRV2`: best 0.333333, mean 0.333333.

Interpretation:

- The diagnostic separates perfect PAM-supported hits from partial or weak
  motif compatibility without changing current binary PAM support.
- This may be useful later for ranking or calibration, but it should remain
  diagnostic until validated against curated PAM/protospacer examples.

Validation:

- Full local test suite passes: 72 tests.

## Latest Scoring Experiment: Experimental PAM-Weighted Score

Recent code updates added an experimental comparator score:
`experimental_pam_weighted_score`.

Purpose:

- Test whether probabilistic PAM/PFS compatibility improves separation of
  validated positives from PAM-failure controls.
- Keep the production `crispr_targeting_score` unchanged.

Implemented:

- New diagnostic function:
  `score_experimental_pam_weighted_evidence`.
- Evidence matrix now includes `experimental_pam_weighted_score`.
- Benchmark evaluation now carries through:
  - `experimental_pam_weighted_score`
  - `best_pam_compatibility_score`
  - `mean_pam_compatibility_score`

Experiment artifact:

- Recomputed run:
  `outputs/runs/20260520_210025_pam_weighted_experiment`
- Joined benchmark:
  `data/curation/benchmark_evaluation_20260520_210025_pam_weighted_experiment.csv`

Current comparison:

- `pa14_jbd18`: current 100.0, experimental 91.0.
- `pa14_jbd25`: current 83.0, experimental 78.0.
- `pa14_jbd67`: current 100.0, experimental 90.0.
- `pa14_dms3`: current 28.0, experimental 53.0.
- `dgcc7710_2972_wt`: current 39.0, experimental 39.0.
- `rey15a_sirv2`: current 99.36, experimental 91.03.
- `ssolfataricus_p2_sirv2`: current 25.5, experimental 39.0.

Interpretation:

- The experimental PAM-weighted score keeps validated PA14 positives high and
  keeps the DGCC7710/Phi2972 PAM-failure control low.
- However, it raises ambiguous `pa14_dms3` to 53.0 because partial PAM
  compatibility of 0.5 receives too much credit.
- Therefore this exact formula should remain diagnostic only and should not
  replace the current production score.
- A future version would need stricter handling of partial PAM compatibility,
  likely requiring curated positive/negative PAM examples before calibration.

Validation:

- Full local test suite passes: 73 tests.

## Current Next Step: CCTyper in WSL

Why this matters:

- CCTyper/CRISPRCasTyper is not needed for SABR runtime.
- It can improve offline curation by providing independent locus-level
  CRISPR-Cas subtype calls, especially for Type III systems where the
  repeat-only model remains weak.
- It should be used as a training-label support and audit tool, not as a
  required user-facing dependency.

Recommended install route:

- Use WSL Ubuntu, not native Windows.
- Install Miniforge or another conda-compatible distribution inside Ubuntu.
- Create a Linux environment with the Bioconda/Russel88 CCTyper package:

```bash
mamba create -n cctyper -c conda-forge -c bioconda -c russel88 cctyper
conda activate cctyper
```

If `mamba` is unavailable, use:

```bash
conda create -n cctyper -c conda-forge -c bioconda -c russel88 cctyper
conda activate cctyper
```

Environment checks:

```bash
cctyper -h
which cctyper
which prodigal
which hmmsearch
which minced
which blastn
```

Small SABR test run from WSL:

```bash
cd /mnt/c/Users/es60/Desktop/Codex
mkdir -p outputs/cctyper_test
cctyper data/examples/resistant/bacteria/Paeruginosa_PA14.fasta outputs/cctyper_test/pa14 --no_plot --simplelog -t 4
```

If the PA14 example path is missing, find available FASTA files with:

```bash
find data -name "*.fasta" -o -name "*.fa" -o -name "*.fna" | head
```

The key CCTyper output for SABR import is usually `crisprs_near_cas.tab`.

Planned SABR import workflow after CCTyper runs:

```powershell
.\.venv\Scripts\python.exe -m crispr_phage_predictor.ml.collect_cctyper_training data\training\cctyper_manifest.csv --output data\training\repeats_cas_types_cctyper.csv
```

Practical next command to inspect the machine state:

```powershell
wsl -l -v
```

Latest WSL installation status:

- WSL Ubuntu is installed and starts correctly:
  - distro: `Ubuntu-24.04`
  - WSL version: `2.7.3.0`
  - kernel: `6.6.114.1-1`
- Conda already exists inside WSL:
  - `conda 26.3.2`
  - path: `/home/es60/miniconda3/bin/conda`
  - `bzip2` is also available inside conda.
- WSL has routing and DNS configured:
  - default route via `10.169.24.1`
  - DNS nameserver `10.255.255.254`
  - search domain `win2k.aub.edu.lb`
- Windows host internet works:
  - PowerShell `Invoke-WebRequest https://www.google.com` returns status 200.
  - `Test-NetConnection www.google.com -Port 443` succeeds.
  - `Test-NetConnection repo.anaconda.com -Port 443` succeeds.
- WSL outbound internet does not work:
  - `curl https://www.google.com` times out.
  - `curl https://repo.anaconda.com` times out.
  - `curl https://conda.anaconda.org/bioconda` times out.
- Tried WSL config route:
  - `.wslconfig` with `networkingMode=mirrored`, `dnsTunneling=true`,
    `autoProxy=true`, and `firewall=false`.
  - User reported this did not fix WSL outbound access.

Current interpretation:

- The CCTyper blocker is not CCTyper itself and not missing conda.
- The blocker is WSL outbound networking, likely due to Windows/enterprise/AUB
  firewall, VPN, endpoint security, or WSL/Hyper-V traffic policy.
- Do not spend more time trying normal `conda create` inside WSL until WSL can
  reach the internet.

Preferred alternatives if WSL networking remains blocked:

1. Run CCTyper on another Linux environment with working internet, such as a
   lab Linux workstation, HPC node, cloud VM, or collaborator machine, then copy
   the `crisprs_near_cas.tab` outputs back into this repo.
2. Use Google Colab or another hosted Linux notebook to install CCTyper with
   conda/mamba, run CCTyper on uploaded FASTA files, then download the output
   folder.
3. Use Docker only if Docker Desktop can pull images or build containers through
   Windows networking; if Docker relies on the same broken WSL network path, it
   may fail too.
4. As a last resort, use Windows to download Linux conda packages or a prepared
   Linux environment for offline transfer into WSL, but this is more fragile
   because compiled bioinformatics tools may have post-install and path issues.

## Latest Repository Organization and Backup

Repository cleanup completed:

- Raw CRISPRCasdb files from the temporary root `crisprdb/` folder were moved
  under `data/training/external_sources/crisprcasdb_34/`.
- Large raw files are ignored by git:
  - CRISPRCasdb ZIP/PDF payloads.
  - downloaded GenBank source folders.
  - downloaded FASTA curation inputs.
  - generated benchmark evaluation CSVs.
  - generated report CSV/DOCX/figure outputs.
- Small source/provenance files remain tracked:
  - code and tests.
  - curated metadata TSVs.
  - README files.
  - report-generation scripts.
  - `assets/aub-logo.png`.

Checkpoint:

- Full local test suite passed: 73 tests.
- Commit created: `5023ddb Organize SABR pipeline and training data sources`.
- GitHub remote configured:
  `https://github.com/ESA-zoph/SABR`
- Local `master` was pushed to `origin/master`.
- Working tree was clean immediately after push.

## Latest Implementation Status: CRISPRCasdb Direct-Repeat Inventory

New CRISPRCasdb release 34 source assessment:

- `dr_34.zip` contains direct-repeat FASTA exports.
- `spacer_34.zip` contains spacer FASTA exports.
- `ccpp_db.zip` contains a PostgreSQL dump with CRISPR locus, sequence, taxon,
  and Cas-cluster tables.
- Local quick counts:
  - `direct_repeat_id.fsa`: 28,712 records.
  - `spacer_id.fsa`: 353,377 records.

Implementation decision:

- Direct-repeat FASTA alone is useful, but it is not a defensible Cas
  type/subtype training label source by itself.
- Add an unlabeled inventory importer first, then later join records against
  SQL-derived Cas-cluster/locus metadata if the relationship is clear.

Implemented:

- `crispr_phage_predictor.ml.import_crisprcasdb_repeats`
  creates an unlabeled direct-repeat inventory from `dr_34.zip`.
- Output columns include source, release, FASTA member, record ID, accession
  count, repeat sequence, repeat length, sequence hash, IUPAC validity, and
  whether the repeat is currently usable by SABR repeat-feature extraction.
- Documentation updated in:
  - `data/training/README.md`
  - `data/training/external_sources/crisprcasdb_34/README.md`

Current command:

```powershell
.\.venv\Scripts\python.exe -m crispr_phage_predictor.ml.import_crisprcasdb_repeats data\training\external_sources\crisprcasdb_34\dr_34.zip --output data\training\crisprcasdb_34_direct_repeats_inventory.csv
```

Local inventory results:

- Default accession-oriented member `direct_repeat_seqName.fsa`:
  - output: `data/training/crisprcasdb_34_direct_repeats_inventory.csv`
  - rows: 7,683
  - currently usable by SABR repeat features: 7,652
  - not usable by current repeat features: 31
- ID-oriented member `direct_repeat_id.fsa`:
  - output: `data/training/crisprcasdb_34_direct_repeats_by_id_inventory.csv`
  - rows: 28,712
- Both output CSV files are derived data under `data/training/*.csv` and remain
  ignored by git.

## Latest Implementation Status: CRISPRCasdb SQL Candidate Labels

SQL dump structure:

- `crisprlocus.drconsensus` points to a direct-repeat `region.id`.
- `crisprlocus.sequence` and `clustercas.sequence` share the same sequence UUID.
- `clustercas.class` contains labels such as `CAS-TypeI-F`,
  `CAS-TypeIII-B`, and ambiguous `CAS`.
- `crisprlocus_region` links each locus to repeat/spacer `region` rows.
- `sequence.strain` links to `strain`, which provides GenBank/RefSeq accessions
  and assembly status.

Implemented:

- `crispr_phage_predictor.ml.import_crisprcasdb_sql`
  builds computational candidate repeat/Cas rows from the extracted PostgreSQL
  dump.
- Conservative filters:
  - CRISPR locus evidence level >= 4 by default.
  - valid direct-repeat consensus sequence.
  - unambiguous `CAS-Type...` cluster subtype.
  - nearest same-sequence Cas cluster within 20,000 bp by default.
  - spacer count and mean spacer length derived from linked region rows.
- Output is in the standard SABR repeat/Cas table schema but must be treated as
  computational candidates, not curated gold labels.

Current command:

```powershell
.\.venv\Scripts\python.exe -m crispr_phage_predictor.ml.import_crisprcasdb_sql data\training\external_sources\crisprcasdb_34\home\pa.charbit\20220414_ccpp_recette_chromo_complete.sql --output data\training\repeats_cas_types_crisprcasdb_sql_candidate.csv
```

Local output:

- `data/training/repeats_cas_types_crisprcasdb_sql_candidate.csv`
- Rows: 23,507.
- File size: about 6.7 MB.
- Top subtype counts:
  - `I-E`: 10,982
  - `I-F`: 2,172
  - `I-C`: 1,847
  - `III-A`: 1,637
  - `I-B`: 1,517
  - `II-C`: 1,451
  - `II-A`: 1,106
  - `III-B`: 833
  - `I-A`: 545
  - `III-D`: 422
- Type counts:
  - Type I: 17,689
  - Type III: 2,944
  - Type II: 2,612
  - Type V: 135
  - Type VI: 127

Interpretation:

- This gives a large, traceable computational candidate source that may help
  Type III coverage.
- It is still not manually curated truth and should be used separately from
  gold seed rows in reports and validation.

## Latest Analysis: CRISPRCasdb Candidate Audit Against Current Training

Audit script:

- `crispr_phage_predictor.ml.audit_crisprcasdb_candidates`
- Compares the current best training table against the CRISPRCasdb SQL
  candidate table.
- Reports:
  - row and unique repeat counts
  - repeat-sequence overlap
  - candidate repeat-to-subtype conflicts
  - novel non-conflicting candidate rows
  - a proposed balanced candidate subset capped by subtype

Command run:

```powershell
.\.venv\Scripts\python.exe -m crispr_phage_predictor.ml.audit_crisprcasdb_candidates data\training\repeats_cas_types_augmented_vink_genbank_targeted.csv data\training\repeats_cas_types_crisprcasdb_sql_candidate.csv --output-dir data\training\audits\crisprcasdb_sql_candidate
```

Audit summary:

- Current training rows: 4,884.
- CRISPRCasdb candidate rows: 23,507.
- Current unique repeat hashes: 2,880.
- Candidate unique repeat hashes: 5,355.
- Overlapping repeat hashes: 2,801.
- Candidate repeat hashes with subtype conflicts: 125.
- Candidate rows after conflict filtering: 22,378.
- Novel non-conflicting candidate rows: 4,563.
- Proposed balanced candidate rows with cap 500/subtype: 3,535.

Most useful balanced candidate additions by subtype:

- `I-E`: 500 capped from 1,528 novel non-conflicting rows.
- `II-C`: 489.
- `I-C`: 447.
- `I-B`: 359.
- `III-A`: 266.
- `III-B`: 253.
- `I-G`: 251; not present in current training.
- `I-F`: 235.
- `III-D`: 180.
- `II-A`: 179.
- `I-A`: 132.
- `I-D`: 95.
- `V-K`: 66; not present in current training.

Conflict examples:

- Some identical repeats map to multiple subtypes, including `I-E;I-F`,
  `I-F;III-B`, `I-B;III-B`, and `I-B;II-C;III-B`.
- These conflict repeat hashes should be excluded from first-pass augmented
  training.

Interpretation:

- CRISPRCasdb SQL candidates are highly overlapping with the current
  Vink/GenBank-derived table, so a naive merge would overcount known repeats.
- The useful next dataset should add only novel, non-conflicting candidates and
  cap overrepresented subtypes.
- Candidate labels should remain explicitly marked as computational nearby-Cas
  labels in reports and model comparisons.

## Latest Modeling Experiment: CRISPRCasdb-Augmented Balanced Dataset

Built filtered augmented table:

- Builder: `crispr_phage_predictor.ml.build_crisprcasdb_augmented_dataset`.
- Current table:
  `data/training/repeats_cas_types_augmented_vink_genbank_targeted.csv`
- Candidate table:
  `data/training/repeats_cas_types_crisprcasdb_sql_candidate.csv`
- Output:
  `data/training/repeats_cas_types_augmented_crisprcasdb_sql_balanced.csv`
- Candidate additions only:
  `data/training/repeats_cas_types_crisprcasdb_sql_balanced_additions.csv`
- Added rows: 3,535.
- Total rows: 8,419.
- Filtering:
  - remove candidate repeats already present in current training
  - remove candidate repeat hashes that map to multiple subtypes
  - cap added rows at 500 per subtype

Evaluation command pattern:

```powershell
.\.venv\Scripts\python.exe -m crispr_phage_predictor.ml.train_classifier <dataset.csv> --split-strategy group_holdout --group-column genus --methods extra_trees --min-class-count 20 --include-medium-confidence
```

Comparison:

- Current best table:
  - rows used after min-class filtering: 4,848
  - ExtraTrees genus-holdout accuracy: 0.9152
  - `III-A` f1: 0.58, recall: 0.43
  - `III-B` f1: 0.59, recall: 0.57
  - `III-D` f1: 0.40, recall: 0.40
- CRISPRCasdb-augmented balanced table:
  - rows used after min-class filtering: 8,378
  - ExtraTrees genus-holdout accuracy: 0.9050
  - `III-A` f1: 0.68, recall: 0.66
  - `III-B` f1: 0.55, recall: 0.48
  - `III-D` f1: 0.38, recall: 0.32

Interpretation:

- The balanced CRISPRCasdb additions improve `III-A` recall substantially.
- Overall accuracy drops by about 1.0 percentage point.
- `III-B`, `III-D`, and rare classes remain fragile.
- Do not replace the current best production model with this augmented model
  yet.
- Next modeling step should test targeted Type III-only augmentation or stricter
  distance/confidence filters rather than adding all balanced candidates.

Tracked comparison artifact:

- `docs/crisprcasdb_augmented_model_comparison.csv`

## Latest Modeling Experiment: Type III-Only CRISPRCasdb Augmentation

Built Type III-targeted augmented table:

- Builder option:
  `--include-subtypes III-A,III-B,III-C,III-D`
- Output:
  `data/training/repeats_cas_types_augmented_crisprcasdb_typeiii_balanced.csv`
- Candidate additions:
  `data/training/repeats_cas_types_crisprcasdb_typeiii_balanced_additions.csv`
- Added rows: 714.
- Total rows: 5,598.

Evaluation:

- rows used after min-class filtering: 5,578
- split: genus holdout
- method: ExtraTrees
- include medium-confidence labels: yes
- accuracy: 0.9004

Type III comparison:

- `III-A`:
  - current best f1/recall: 0.58 / 0.43
  - all-subtype CRISPRCasdb f1/recall: 0.68 / 0.66
  - Type III-only CRISPRCasdb f1/recall: 0.76 / 0.86
- `III-B`:
  - current best f1/recall: 0.59 / 0.57
  - all-subtype CRISPRCasdb f1/recall: 0.55 / 0.48
  - Type III-only CRISPRCasdb f1/recall: 0.49 / 0.57
- `III-D`:
  - current best f1/recall: 0.40 / 0.40
  - all-subtype CRISPRCasdb f1/recall: 0.38 / 0.32
  - Type III-only CRISPRCasdb f1/recall: 0.33 / 0.26

Interpretation:

- Type III-only augmentation strongly improves `III-A`.
- It does not improve `III-B` or `III-D`, and overall accuracy drops more than
  the all-subtype augmentation.
- This suggests CRISPRCasdb candidates are useful for `III-A` specifically, but
  not enough to resolve Type III broadly.
- Next experiment should try `III-A`-only augmentation or stricter nearby-Cas
  distance thresholds before considering any production model change.

## Latest Documentation and Visualization: Manuscript Draft and Feature-Space Projections

New manuscript draft:

- `docs/SABR_manuscript_draft.md`
- Scope:
  - SABR scientific framing as a CRISPR targeting evidence mapper, not a direct
    resistance caller.
  - current pipeline methods.
  - repeat-derived Cas subtype model development.
  - current best predictor results.
  - CRISPRCasdb direct-repeat and SQL candidate experiments.
  - probability/confidence interpretation.
  - limitations and future work.
- Tone is intentionally cautious and honest: computational candidates are not
  described as gold-standard labels, and spacer matches are not described as
  proof of resistance.

New dimensionality-reduction script:

- `docs/generate_dataset_reduction_figures.py`
- Uses current scikit-learn dependencies only.
- Generates:
  - PCA projections for all rows.
  - sampled t-SNE projections for local feature-space structure.
  - subtype-colored plots.
  - dataset-source-colored plots.
  - coordinate CSVs for downstream inspection.
- UMAP was not added yet to avoid introducing a new dependency; can be added
  later with `umap-learn`.

Generated assets:

- `docs/manuscript_assets/current_plus_crisprcasdb_balanced_pca_by_subtype.png`
- `docs/manuscript_assets/current_plus_crisprcasdb_balanced_pca_by_dataset_group.png`
- `docs/manuscript_assets/current_plus_crisprcasdb_balanced_tsne_by_subtype.png`
- `docs/manuscript_assets/current_plus_crisprcasdb_balanced_tsne_by_dataset_group.png`
- `docs/manuscript_assets/current_plus_crisprcasdb_typeiii_pca_by_subtype.png`
- `docs/manuscript_assets/current_plus_crisprcasdb_typeiii_pca_by_dataset_group.png`
- `docs/manuscript_assets/current_plus_crisprcasdb_typeiii_tsne_by_subtype.png`
- `docs/manuscript_assets/current_plus_crisprcasdb_typeiii_tsne_by_dataset_group.png`
- matching coordinate CSV files are also written in `docs/manuscript_assets/`.

Validation:

- Full test suite passes: 81 tests.

## Latest Modeling Experiment: CRISPRCasdb-Only Training Control

Question:

- Does a model trained solely on CRISPRCasdb SQL candidate rows perform better
  than the current SABR training table?

Rationale:

- This is an important internal consistency control for CRISPRCasdb.
- It does not prove the labels are biologically superior because the table is
  computationally derived from CRISPRCasdb nearest same-sequence Cas clusters.

Evaluation settings:

- Candidate table:
  `data/training/repeats_cas_types_crisprcasdb_sql_candidate.csv`
- Current table:
  `data/training/repeats_cas_types_augmented_vink_genbank_targeted.csv`
- split: `group_holdout`
- group column: `genome_id`
- min class count: 20
- include medium-confidence labels: yes

Current table, genome-holdout:

- nearest repeat:
  - rows used: 4,848
  - accuracy: 0.9084
  - `III-A` f1/recall: 0.66 / 0.69
  - `III-B` f1/recall: 0.60 / 0.60
  - `III-D` f1/recall: 0.27 / 0.32
- ExtraTrees:
  - rows used: 4,848
  - accuracy: 0.9266
  - `III-A` f1/recall: 0.77 / 0.75
  - `III-B` f1/recall: 0.63 / 0.52
  - `III-D` f1/recall: 0.15 / 0.11

CRISPRCasdb SQL candidate-only table, genome-holdout:

- nearest repeat:
  - rows used: 23,478
  - accuracy: 0.9389
  - `III-A` f1/recall: 0.86 / 0.86
  - `III-B` f1/recall: 0.52 / 0.50
  - `III-D` f1/recall: 0.54 / 0.49
- ExtraTrees:
  - rows used: 23,478
  - accuracy: 0.9455
  - `III-A` f1/recall: 0.89 / 0.88
  - `III-B` f1/recall: 0.58 / 0.53
  - `III-D` f1/recall: 0.54 / 0.43

Interpretation:

- The selected-model training annotation table performs better than the earlier
  development table under genome-holdout
  internal validation.
- The nearest-repeat baseline is also very high, suggesting strong repeat-family
  structure and possible database-specific consistency.
- This makes CRISPRCasdb valuable as a model-development source, but it should
  not automatically replace the current production training table without:
  - external validation against curated/literature labels
  - genus/species-level validation with organism metadata restored
  - conflict filtering and calibration
  - checking whether performance transfers to uploaded genomes outside the
    CRISPRCasdb distribution

## Latest Modeling Experiment: External Dataset Transfer Between Current and CRISPRCasdb Tables

Question:

- If the selected-model annotation training table looks better internally, does it transfer to the earlier
  SABR training table?
- Conversely, does the current model transfer to CRISPRCasdb candidate rows?

Implemented:

- `crispr_phage_predictor.ml.evaluate_external_dataset`
- Trains ExtraTrees on one repeat/Cas table and evaluates on another.
- Test rows with labels absent from the training classes are excluded and
  reported.

CRISPRCasdb train -> current table test:

- Training table:
  `data/training/repeats_cas_types_crisprcasdb_sql_candidate.csv`
- Train rows after min-class filtering: 23,478.
- Raw current-table test rows: 4,884.
- Evaluated current-table test rows: 4,806.
- Excluded rows: 78.
- Excluded subtypes: `I-U`, `I-V`, `V-B`, `VI-A`, `VI-C`.
- Accuracy: 0.9811.
- Type III:
  - `III-A` f1/recall: 0.91 / 0.94
  - `III-B` f1/recall: 0.93 / 0.96
  - `III-D` f1/recall: 0.89 / 0.87

Current table train -> CRISPRCasdb test:

- Training table:
  `data/training/repeats_cas_types_augmented_vink_genbank_targeted.csv`
- Train rows after min-class filtering: 4,848.
- Raw CRISPRCasdb test rows: 23,507.
- Evaluated CRISPRCasdb test rows: 22,981.
- Excluded rows: 526.
- Excluded subtypes: `I-G`, `III-C`, `V-B`, `V-B1`, `V-F4`, `V-K`,
  `VI-A`, `VI-C`.
- Accuracy: 0.9559.
- Type III:
  - `III-A` f1/recall: 0.93 / 0.90
  - `III-B` f1/recall: 0.69 / 0.57
  - `III-D` f1/recall: 0.70 / 0.59

Interpretation:

- The selected SABR ExtraTrees architecture trained on CRISPRCasdb-derived
  annotations transfers very well to rows in the earlier
  table whose subtypes exist in CRISPRCasdb training.
- The current table also transfers well to CRISPRCasdb rows, but is weaker for
  `III-B` and `III-D`.
- These numbers are encouraging but not a fully independent biological
  validation because the current table and CRISPRCasdb candidates share
  database ancestry.
- The result supports training a selected SABR runtime candidate from the
  CRISPRCasdb-derived annotation table
  artifact, but final replacement should wait for curated/literature/CCTyper
  external validation and probability calibration.

## Current Project Interpretation: What SABR Contributes

SABR appears to bring a real contribution, but the claim should be framed
carefully.

Strong contribution:

- SABR is an integrated CRISPR-phage targeting evidence mapper, not a direct
  resistance caller.
- It combines:
  - FASTA parsing and diagnostics
  - CRISPR array detection
  - spacer extraction
  - spacer-phage matching
  - PAM/PFS support diagnostics
  - seed-region mismatch summaries
  - repeat-derived Cas subtype prediction
  - bacteria-by-phage evidence matrices
  - Streamlit GUI and saved reproducible outputs
- It uses cautious scientific language:
  - CRISPR targeting evidence
  - candidate CRISPR-phage interaction
  - PAM/PFS compatibility evidence
  - repeat-derived Cas subtype prediction
- It avoids the common overclaim that spacer matches prove biological
  resistance.
- Its Cas subtype model uses FASTA-derived repeat/array features rather than
  runtime taxonomy, organism name, source database, or cas-gene shortcuts.
- The project now has transparent validation artifacts:
  - nearest-repeat baselines
  - ExtraTrees/random forest/hybrid/hierarchical comparisons
  - genome/genus holdout splits
  - controls using only the larger computational annotation-derived table
  - cross-dataset transfer tests
  - dimensionality-reduction figures
  - explicit negative results from naive data augmentation

Current best predictor interpretation:

- Safest current production candidate remains the flat ExtraTrees model trained
  on `repeats_cas_types_augmented_vink_genbank_targeted.csv`.
- Best documented genus-holdout accuracy for this model: 0.9152.
- The selected SABR ExtraTrees model trained on the larger computational
  annotation-derived table gives stronger internal genome-holdout accuracy
  at 0.9455 and transfers well to the current table, but it should remain an
  experimental candidate until externally validated against curated/literature
  or CCTyper-supported labels.

What SABR should not claim yet:

- confirmed phage resistance prediction
- clinical or phenotypic resistance prediction
- superiority over all existing CRISPR tools
- fully calibrated subtype probabilities
- final gold-standard Type III-B/III-D subtype resolution

Best manuscript framing:

- workflow integration
- cautious evidence reporting
- FASTA-only repeat-derived Cas subtype prediction
- reproducible model/data validation
- evidence that more computational training rows can improve one class while
  reducing overall performance
- transparent limitations and need for curated external validation

## Recommended Next Step

The next best technical step is probability calibration and external-validation
packaging, before changing the production runtime model.

Priority order:

1. Add calibration analysis for the current best ExtraTrees model:
   - reliability curve
   - confidence bins
   - expected calibration error
   - accuracy by confidence threshold
   - per-subtype confidence behavior
2. Export a selected SABR model artifact trained on the larger computational
   annotation-derived table separately from the
   production artifact:
   - keep `models/cas_subtype_extratrees.joblib` as production for now
   - write something like `models/cas_subtype_extratrees_crisprcasdb_experimental.joblib`
   - compare metadata and predictions without silently replacing runtime
3. Update the manuscript with:
   - clear contribution statement
   - model comparison table
   - PCA/t-SNE figures
   - calibration figure once available
4. Build a small independent validation panel:
   - manually curated/literature rows
   - CCTyper-supported labels if WSL/networking can be solved or run elsewhere
   - focus on Type III-B and Type III-D
5. Only after calibration and independent validation decide whether the
   selected SABR candidate trained on the larger annotation-derived table should
   replace the earlier production model.

## Latest Calibration and Experimental Model Artifact

Implemented:

- `crispr_phage_predictor.ml.calibrate_subtype_model`
- Outputs:
  - `summary.csv`
  - `predictions.csv`
  - `confidence_bins.csv`
  - `subtype_confidence.csv`
  - `accuracy_by_threshold.csv`
  - `reliability_curve.png`
  - `accuracy_by_confidence_threshold.png`
  - `subtype_confidence_accuracy.png`

Current production-candidate calibration:

- Dataset:
  `data/training/repeats_cas_types_augmented_vink_genbank_targeted.csv`
- Split: genus holdout.
- Train rows: 3,633.
- Test rows: 1,215.
- Accuracy: 0.915226.
- Mean predicted confidence: 0.689811.
- Expected calibration error: 0.225416.
- Interpretation: the model is substantially under-confident on this split;
  observed accuracy is higher than predicted confidence across most confidence
  bins.
- Output directory:
  `docs/calibration/current_best/`

Selected SABR model calibration on its annotation-derived internal holdout:

- Dataset:
  `data/training/repeats_cas_types_crisprcasdb_sql_candidate.csv`
- Split: genome_id holdout.
- Train rows: 17,607.
- Test rows: 5,871.
- Accuracy: 0.945495.
- Mean predicted confidence: 0.922668.
- Expected calibration error: 0.036156.
- Interpretation: the selected SABR ExtraTrees model is much better calibrated
  internally than the current production-candidate model, but this remains an
  internal computational-label validation.
- Output directory:
  `docs/calibration/crisprcasdb_only/`

Experimental model artifact:

- Command:

```powershell
.\.venv\Scripts\python.exe -m crispr_phage_predictor.ml.export_model --training-table data\training\repeats_cas_types_crisprcasdb_sql_candidate.csv --output models\cas_subtype_extratrees_crisprcasdb_experimental.joblib --min-class-count 20
```

- Output:
  `models/cas_subtype_extratrees_crisprcasdb_experimental.joblib`
- Training rows: 23,478.
- Classes:
  `I-A`, `I-B`, `I-C`, `I-D`, `I-E`, `I-F`, `I-G`, `II-A`, `II-B`,
  `II-C`, `III-A`, `III-B`, `III-C`, `III-D`, `V-A`, `V-K`, `VI-B1`.
- The artifact is local and ignored by git like other model binaries.
- Do not replace `models/cas_subtype_extratrees.joblib` yet.

Validation:

- Full local test suite passes: 82 tests.

## Latest Runtime Decision: SABR Uses CRISPRCasdb-Trained Model Locally

Runtime model decision:

- The local default runtime artifact `models/cas_subtype_extratrees.joblib` was
  regenerated from:
  `data/training/repeats_cas_types_crisprcasdb_sql_candidate.csv`
- Training rows after min-class filtering: 23,478.
- Classes:
  `I-A`, `I-B`, `I-C`, `I-D`, `I-E`, `I-F`, `I-G`, `II-A`, `II-B`,
  `II-C`, `III-A`, `III-B`, `III-C`, `III-D`, `V-A`, `V-K`, `VI-B1`.
- The model binary is about 410 MB and remains ignored by git under
  `models/*.joblib`.
- Source defaults were updated so future `export_model` runs also use the
  CRISPRCasdb SQL candidate table by default:
  - `crispr_phage_predictor/ml/model_artifact.py`
  - fallback nearest-repeat table in `crispr_phage_predictor/cas_prediction.py`

Caveat:

- This is now the local runtime model because it has the best internal and
  transfer performance so far.
- Manuscript language should describe it as the selected SABR ExtraTrees
  subtype model trained on computationally
  derived until independent curated or CCTyper-supported validation is added.

## Latest Manuscript Export

Generated manuscript Word document:

- `docs/SABR_manuscript.docx`
- Embeds PNG versions of the main figures for Word compatibility.
- Matching SVG files are stored under `docs/manuscript_assets/`.

New manuscript helper scripts:

- `docs/generate_manuscript_summary_figures.py`
- `docs/generate_manuscript_docx.py`

SVG manuscript figures generated:

- `benchmark_score_summary.svg`
- `model_comparison_summary.svg`
- `calibration_comparison.svg`
- `typeiii_performance_summary.svg`
- PCA subtype/source SVGs.
- sampled t-SNE subtype/source SVGs.

Validation:

- Full local test suite passes: 83 tests.

## Latest Error Projection: Where Wrong Calls Cluster

New script:

- `docs/generate_error_projection_figures.py`

Purpose:

- Project the selected-model held-out annotation-table test set into PCA and sampled t-SNE
  feature space.
- Overlay correct versus wrong calls.
- Plot wrong calls by true subtype.
- Summarize most common subtype confusions.

Generated files:

- `docs/manuscript_assets/crisprcasdb_model_error_projection_coordinates.csv`
- `docs/manuscript_assets/crisprcasdb_model_error_rows.csv`
- `docs/manuscript_assets/crisprcasdb_model_pca_correct_vs_wrong.png/.svg`
- `docs/manuscript_assets/crisprcasdb_model_tsne_correct_vs_wrong.png/.svg`
- `docs/manuscript_assets/crisprcasdb_model_tsne_errors_by_true_subtype.png/.svg`
- `docs/manuscript_assets/crisprcasdb_model_top_error_pairs.png/.svg`

Error summary:

- Wrong calls are concentrated in Type III and adjacent Type I regions.
- Top wrong-call pairs:
  - `III-B -> I-B`: 37
  - `III-A -> I-B`: 25
  - `III-D -> III-A`: 18
  - `III-B -> I-C`: 17
  - `III-B -> I-A`: 16
  - `I-B -> III-B`: 14
  - `III-D -> III-B`: 14
  - `III-A -> III-B`: 13
  - `III-B -> III-A`: 12
- Wrong calls by true subtype:
  - `III-B`: 94
  - `III-D`: 56
  - `III-A`: 49
  - `I-B`: 28
  - `I-A`: 22
  - `III-C`: 16
  - `I-D`: 16

Interpretation:

- The model is not failing randomly.
- Most residual errors are in biologically/plausibly overlapping repeat-feature
  regions involving Type III and some Type I subtypes.
- This supports targeted future validation and data expansion for Type III-B,
  Type III-D, and neighboring Type I-B/I-C/I-A boundaries.

## Latest Model Interpretability

New script:

- `docs/generate_model_interpretability.py`

Generated interpretability tables:

- `docs/interpretability/crisprcasdb_extratrees/builtin_feature_importance.csv`
- `docs/interpretability/crisprcasdb_extratrees/permutation_feature_importance.csv`
- `docs/interpretability/crisprcasdb_extratrees/feature_category_importance.csv`
- `docs/interpretability/crisprcasdb_extratrees/typeiii_error_feature_summary.csv`

Generated manuscript assets:

- `docs/manuscript_assets/crisprcasdb_builtin_feature_importance.png/.svg`
- `docs/manuscript_assets/crisprcasdb_permutation_feature_importance.png/.svg`
- `docs/manuscript_assets/crisprcasdb_feature_category_importance.png/.svg`
- `docs/manuscript_assets/crisprcasdb_typeiii_error_feature_summary.png/.svg`

Top built-in ExtraTrees features:

- `repeat_length`
- `spacer_repeat_length_ratio`
- `mean_spacer_length`
- `min_spacer_length`
- `max_spacer_length`
- `median_spacer_length`
- terminal repeat k-mers such as `terminal_kmer_4_CAAC`, `terminal_kmer_4_GTTG`, and `terminal_kmer_4_ACGT`
- `repeat_at_percent`, `repeat_gc_percent`, and `repeat_hairpin_score`

Category-level interpretation:

- Built-in importance is dominated by whole-repeat k-mers and terminal k-mers.
- Smaller but interpretable contributions come from terminal composition, array statistics, repeat composition, and repeat-structure features.
- Held-out permutation drops are small because many repeat-derived features are correlated; single-feature permutation ranks should be interpreted cautiously.
- The defensible conclusion is that the model uses repeat motifs, terminal motifs, array geometry, and simple repeat-structure proxies, not organism names, taxonomy, database source, accession, or cas-gene annotations as runtime inputs.

Manuscript update:

- `docs/SABR_manuscript_draft.md` now reflects the CRISPRCasdb runtime candidate and interpretability results.
- `docs/generate_manuscript_docx.py` now embeds interpretability text and figures.
- `docs/SABR_manuscript.docx` was regenerated with the new figures.

Validation:

- Full local test suite passes: 83 tests.

## Recommended Next Step

The next best scientific step is independent validation, not more internal
optimization using only the annotation-derived development source.

Priority order:

1. Build a CCTyper/literature-supported validation table focused on independent
   labels, especially `III-B`, `III-D`, and adjacent Type I-B/I-C/I-A cases.
2. Evaluate the selected SABR ExtraTrees subtype model on that independent
   table and report accuracy, calibration, and subtype-level errors.
3. If CCTyper remains blocked locally, run it on Linux/HPC/Colab and import the
   output with the existing SABR CCTyper collectors.
4. Expand the small resistance/sensitivity benchmark panel separately, keeping
   it framed as CRISPR targeting evidence validation rather than direct
   phenotypic resistance prediction.

## Latest Docker Runtime and UX Validation

Date:

- 2026-06-06

Docker status:

- Docker Desktop on Windows is working after `wsl --shutdown` resolved a
  temporary Docker Desktop Linux engine `500 Internal Server Error`.
- The Docker image now uses `python:3.11-slim-bookworm`.
- The Dockerfile installs:
  - Debian `ncbi-blast+`, exposing `blastn` and `makeblastdb`.
  - `requirements-external.txt`, currently `diced>=0.1.3`.
- Container-level verification showed:
  - `blastn`: available.
  - `makeblastdb`: available.
  - `diced`: available.
- Runtime model mounting works with:

```powershell
docker run --rm -p 7860:7860 `
  -v "${PWD}\models\cas_subtype_extratrees.joblib:/app/models/cas_subtype_extratrees.joblib:ro" `
  sabr:local
```

- The app reports the mounted selected ExtraTrees artifact as available
  at about `391.2 MB`.

Recommended Docker run command for persistent outputs:

```powershell
docker run --rm -p 7860:7860 `
  -v "${PWD}\models\cas_subtype_extratrees.joblib:/app/models/cas_subtype_extratrees.joblib:ro" `
  -v "${PWD}\outputs:/app/outputs" `
  sabr:local
```

Reason:

- Without the `/app/outputs` bind mount, completed Streamlit analyses are saved
  inside the short-lived Docker container rather than under the Windows
  workspace `outputs/runs/`.

User-facing runtime behavior:

- Sidebar now reports:
  - MinCED-compatible detection as `diced` when available.
  - BLAST+ as available when `blastn` and `makeblastdb` are installed.
  - Selected ExtraTrees artifact availability and size.
- `Auto recommended` now selects MinCED-compatible detection when Diced or a
  `minced` command is available, and falls back to the internal detector when
  no MinCED-compatible backend exists.
- BLASTN is selected automatically when available for spacer matching.
- The internal exact-repeat detector remains available manually as a transparent
  baseline.

Fragmented FASTA safeguard:

- A highly fragmented bacterial FASTA, `Ecol167.fasta`, was observed with
  `76,543` FASTA records.
- Internal CRISPR detection initially displayed a very large estimated runtime
  because progress was only updated after each record and the input contained
  tens of thousands of records.
- The Streamlit app now warns when a bacterial upload has at least `1,000`
  FASTA records and blocks internal CRISPR detection when any bacterial file
  has at least `5,000` records.
- User guidance is to use a better assembled genome, filter short contigs, or
  use an external detector before retrying very fragmented assemblies.

New helper script:

- `scripts/filter_fasta_records.py`
- Purpose:
  - filter a fragmented FASTA to longer records for SABR testing.
  - sort retained records by length.
  - optionally cap to the longest N records.
- Example:

```powershell
.\.venv\Scripts\python.exe scripts\filter_fasta_records.py Ecol167.fasta data\examples\Ecol167_filtered.fasta --min-length 10000 --max-records 500
```

Progress-bar improvement:

- Internal CRISPR detection now reports scan-level progress during repeat-length
  scanning.
- A long single-record genome no longer sits indefinitely at
  `Preparing CRISPR array detection...`; it advances through scan steps.
- The app progress text now uses `CRISPR detection: X/Y scan steps`.

Latest Docker smoke tests:

- Tiny `data/examples/basic_demo/` run succeeded.
- Real example run with Diced/BLASTN completed very quickly and produced
  plausible positive targeting-evidence results.
- Latest inspected Docker-internal run:
  `/app/outputs/runs/20260606_192711`.
- Run configuration:
  - CRISPR detection: MinCED-compatible (`detection_method=minced`,
    `detection_backend_detail=diced`).
  - Spacer matching: BLASTN.
  - BLAST minimum identity: `0.9`.
  - BLAST minimum coverage: `0.95`.
  - Full-query requirement: `false`.
  - PAM mode: auto from predicted subtype.
  - Cas predictions: 15.
- Runtime:
  - CRISPR detection: `0.233` seconds.
  - BLAST matching: `0.510` seconds.
  - Total: `0.829` seconds.
- Input scale:
  - 4 bacterial files/sequences.
  - 8 phage files/sequences.
  - Duplicate records excluded: 2 bacterial, 5 phage.
- Output:
  - 15 candidate arrays.
  - 715 extracted spacers.
  - 22 spacer hits.
  - 32 evidence-matrix rows.
- Top positive evidence rows:
  - `Paeruginosa_PA14.fasta` vs `PhageJBD18.fasta`:
    score `99.97`, 5 unique matching spacers, PAM compatible.
  - `Paeruginosa_PA14.fasta` vs `PhageJBD67.fasta`:
    score `99.97`, 4 unique matching spacers, PAM compatible.
  - `Sislandicus_REY15A.fasta` vs `SIRV2.fasta`:
    score `99.36`, 5 unique matching spacers, PAM compatible.
  - `Paeruginosa_PA14.fasta` vs `PhageJBD25.fasta`:
    score `82.95`, 1 unique matching spacer, PAM compatible.

Validation:

- Focused compile check passed for:
  - `app.py`
  - `crispr_phage_predictor/crispr.py`
  - `crispr_phage_predictor/pipeline.py`
  - `scripts/filter_fasta_records.py`
- Focused tests passed:
  - `tests/test_crispr_detection.py`
  - `tests/test_pipeline_outputs.py`
  - result: 8 passed.

Latest hit-detail UI addition:

- Added a first-pass hit-detail workflow for bacteria-phage heatmap cells.
- The custom heatmap itself is still static HTML, but the results page now
  provides a `Select heatmap cell` control listing nonzero bacteria-phage cells
  by score and spacer count plus an `Open hit details` button.
- Pressing `Open hit details` switches the Streamlit app to a dedicated
  hit-details view using session state, with a `Back to results` button.
- Follow-up navigation hardening:
  - selected heatmap cell is stored immediately when chosen
  - `Open hit details` also sets the URL/query state to `view=hit_details`
  - `Back to results` sets URL/query state to `view=results`
  - this makes the dedicated view less dependent on a single transient button
    event and easier to diagnose from the browser URL.
- Further fix after the details button appeared to reset the page:
  - latest run ID is stored after each analysis
  - `Open hit details` writes `view`, `bacterium`, `phage`, and `run` query
    parameters
  - if Streamlit reruns or reloads before in-memory results are available, the
    detail page reloads `evidence_matrix.csv` and `spacer_hits.csv` from the
    saved `outputs/runs/{run_id}` folder
  - this is why future Docker runs should include the persistent outputs bind
    mount.
- For the selected bacteria-phage pair, the dedicated page displays:
  - targeting score
  - unique matching spacers
  - total spacer hits
  - PAM/PFS summary
  - hit-level table with array ID, spacer ID, phage contig, coordinates, strand,
    identity, coverage, predicted subtype, PAM/PFS rule/sequence/support, and
    seed mismatches
  - per-hit spacer/protospacer alignment text
  - phage-reference context with the protospacer bracketed between upstream and
    downstream flanks
  - oriented 5-prime/3-prime protospacer flanks and genomic upstream/downstream
    flanks
- Circos-style genome visualization and optional GenBank annotation upload are
  intentionally deferred until the hit-detail workflow is stable.
- Validation:
  - `python -m py_compile app.py` passed.
  - Docker image rebuilt successfully.

## Latest Independent Validation Setup: VirtualBox CCTyper Run

Date:

- 2026-05-25

Linux/CCTyper environment:

- A VirtualBox Ubuntu environment is available for running
  CRISPRCasTyper/CCTyper independently of the Windows SABR runtime.
- CRISPRCasTyper version `1.8.0` was installed in a conda environment named
  `cctyper`.
- Initial startup failed because the launcher could not import
  `pkg_resources`; reinstalling/downgrading `setuptools` restored the launcher.
- `cctyper --help` now runs successfully.
- The recurring `pkg_resources` deprecation warning and a pandas
  `FutureWarning` during locus connection are dependency warnings and did not
  prevent successful output generation.

Input panel:

- The VM input folder is `~/cctyper/bacteria/`.
- It currently contains 19 bacterial FASTA files copied into one flat folder.
- Because the original resistant/susceptible subfolder organization was not
  retained in the VM folder, this run is suitable for independent
  CCTyper-supported Cas subtype validation, but not yet for direct
  resistant-versus-susceptible outcome comparisons.

Successful pilot run:

- `Sthermo_DGCC7710.fasta` was successfully processed with:

```bash
cctyper -t 4 bacteria/Sthermo_DGCC7710.fasta cctyper_runs/Sthermo_DGCC7710_test2 --simplelog
```

- The run completed ORF prediction, HMMER Cas profiling, MinCED CRISPR
  prediction, CRISPR-near-Cas BLAST analysis, repeat subtype prediction,
  locus connection, and plot generation.

Batch result:

- Main output directory: `cctyper_runs_batch1/`.
- Manifest initially created for all 19 FASTA inputs:
  `cctyper_manifest_batch1.csv`.
- Eleven batch-1 genomes generated `crisprs_near_cas.tab`, which are usable by
  the existing SABR CCTyper importer.
- The eight genomes without an importable batch-1
  `crisprs_near_cas.tab` were:
  - `aeruginosa_PAO1`
  - `aureus_NCTC_8325`
  - `cholerae_O1_biovar_El_Tor_N16`
  - `coli_K-12_MG1655_F+_deriv`
  - `faecalis_V583`
  - `monocytogenes_EGD-e`
  - `smegmatis_mc2_155`
  - `subtilis_subsp._subtilis_168`
- These eight were rerun under `cctyper_runs_batch2/`.
- Rerun interpretation:
  - `aeruginosa_PAO1`, `cholerae_O1_biovar_El_Tor_N16`,
    `faecalis_V583`, `smegmatis_mc2_155`, and
    `subtilis_subsp._subtilis_168` reported `No CRISPRs found.`
  - `coli_K-12_MG1655_F+_deriv` reported `No operons found.`
  - `aureus_NCTC_8325` and `monocytogenes_EGD-e` completed analysis but did
    not produce `crisprs_near_cas.tab` usable by the current importer.
  - Batch 2 therefore produced zero additional `crisprs_near_cas.tab` files.

New local SABR utility:

- Added `crispr_phage_predictor/ml/run_cctyper_batch.py`.
- Purpose:
  - scan a bacterial FASTA folder, including optional
    `resistant/`/`susceptible/` subfolders
  - execute CCTyper for each FASTA
  - write a `cctyper_manifest.csv` compatible with
    `crispr_phage_predictor.ml.collect_cctyper_training`
  - write a batch execution summary
- Added focused tests at `tests/test_ml_run_cctyper_batch.py`.
- Validation: 2 focused tests pass.

Immediate import step:

- In Ubuntu, create `cctyper_manifest_usable.csv` containing only the eleven
  batch-1 genomes whose output directories include
  `crisprs_near_cas.tab`.
- Copy `cctyper_runs_batch1/` and `cctyper_manifest_usable.csv` into
  `data/curation/` in the Windows SABR workspace.
- Then import the CCTyper-supported rows with:

```bash
python -m crispr_phage_predictor.ml.collect_cctyper_training data/curation/cctyper_manifest_usable.csv --output data/curation/cctyper_validation_repeats_cas_types.csv
```

Interpretation:

- This gives SABR its first independent CCTyper-supported evaluation material
  for repeat-derived Cas subtype comparisons.
- It is not yet an independent resistance/susceptibility validation panel,
  because the VM copy lost original phenotype grouping and CCTyper itself
  annotates CRISPR-Cas loci rather than phage-resistance phenotype.

## Latest CCTyper Import and External Subtype Evaluation

Imported files and outputs:

- Copied VM result directories into the SABR workspace root:
  - `cctyper_runs_batch1/`
  - `cctyper_runs_batch2/`
- Copied original all-input manifest:
  `cctyper_manifest_batch1.csv`.
- Generated the filtered, importer-safe manifest:
  `data/curation/cctyper_manifest_usable.csv`.
- The filtered manifest contains 11 genomes with available
  `crisprs_near_cas.tab` results.
- Imported validation table:
  `data/curation/cctyper_validation_repeats_cas_types.csv`.

Imported CCTyper-supported validation rows:

- Total array rows: 25.
- Subtype counts:
  - `I-E`: 11
  - `I-A`: 8
  - `I-F`: 3
  - `II-A`: 1
  - `II-C`: 1
  - `III-A`: 1

External evaluation against the selected SABR ExtraTrees subtype model:

- Training source:
  `data/training/repeats_cas_types_crisprcasdb_sql_candidate.csv`.
- Training rows after subtype count filtering: 23,478.
- External CCTyper test rows: 25.
- Evaluated rows: 25.
- Excluded rows/subtypes: 0.
- Accuracy: `1.0000` (25/25 rows).
- Each represented subtype was classified correctly:
  - `I-A`: 8/8
  - `I-E`: 11/11
  - `I-F`: 3/3
  - `II-A`: 1/1
  - `II-C`: 1/1
  - `III-A`: 1/1

Interpretation and limitation:

- This is an encouraging first independent CCTyper-supported evaluation of the
  selected SABR repeat-derived ExtraTrees subtype model.
- The validation panel is small and strongly concentrated in common Type I
  classes.
- It contains only one `III-A` row and no `III-B` or `III-D` rows, so it does
  not yet address the model's main weak/confusable Type III regions.
- These results support the model's current use as a subtype-prediction
  candidate, but do not establish general external accuracy or validate
  phage-resistance predictions.

Small importer correction:

- `crispr_phage_predictor/ml/collect_cctyper_training.py` now preserves blank
  optional manifest metadata as empty fields rather than importing them as
  the string `nan`.
- Added regression coverage in
  `tests/test_ml_collect_cctyper_training.py`.
- Focused CCTyper tests pass: 4 tests.

## Manuscript Framing and Advanced Figure Revision

Date:

- 2026-05-27

Required model terminology correction:

- From this point forward, describe the deployed/candidate subtype classifier as
  the **selected SABR ExtraTrees subtype model**.
- Do not describe it as a "CRISPRCasdb model" or imply that SABR is being
  compared against a CRISPR database model.
- `CRISPRCasdb` is a provenance label for computational training annotations
  used during SABR model development, not a competing predictive method.
- Avoid calling the model "optimal" until the enriched independent Type III /
  adjacent Type I validation panel is complete; "selected" is the defensible
  manuscript term at the current evidence level.

Manuscript revision:

- Rewrote `docs/SABR_manuscript_draft.md` as a fuller scientific first draft.
- Regenerated `docs/SABR_manuscript.docx` from that draft.
- The manuscript now covers:
  - SABR targeting-evidence framing and implemented workflow;
  - training-label provenance and model-selection logic;
  - architecture selection, candidate-annotation experiments, selected-model
    performance, calibration, error analysis, and CCTyper pilot validation;
  - explicit limitations and the ongoing enriched independent validation step.

Equations now included in the manuscript:

- Targeting-evidence score terms:
  - `S_spacer = min(35, 18 * log2(n_s + 1))`
  - `S_identity = min(25, 25 * I_max * C_max)`
  - piecewise `S_PAM` and `S_seed`
  - `S_subtype = min(10, 10 * q_max)` when confidence is available
  - `S_target = min(C_PAM, max(0, S_spacer + S_identity + S_PAM + S_seed + S_subtype))`
  - `C_PAM = 39` when PAM/PFS is evaluated but unsupported for every hit,
    otherwise `C_PAM = 100`
- Evaluation formulas:
  - subtype precision, recall, and F1
  - expected calibration error (ECE)

New manuscript figure generator:

- Added `docs/generate_manuscript_advanced_figures.py`.
- New figure assets in `docs/manuscript_assets/`:
  - `targeting_score_heatmap.png/.svg`
  - `architecture_selection_roc.png/.svg`
  - `selected_model_roc_focus_subtypes.png/.svg`
  - `selected_model_confusion_heatmap.png/.svg`
  - `selected_model_tsne_group_and_subtype_panels.png/.svg`
  - `architecture_selection_roc_metrics.csv`

Revised figure-generation naming:

- Existing manuscript-facing generators now label the classifier as the
  selected SABR model and label datasets as development/annotation sources
  rather than naming a database as if it were a separate model.
- New neutral selected-model interpretability and error assets were generated
  under `docs/manuscript_assets/selected_model_*` and
  `docs/interpretability/selected_sabr_extratrees/`.
- The Word manuscript now embeds a 13-figure sequence including targeting
  evidence heatmap, ROC curves, normalized confusion heatmap, t-SNE panels,
  calibration, feature interpretation, error pairs, and the CCTyper pilot.

Architecture-selection ROC results:

- Common genus-held-out development split, probability-capable models only.
- Macro AUC is averaged over test-evaluable one-vs-rest subtypes.
- Logistic regression:
  - accuracy `0.8337`
  - micro AUC `0.9515`
  - macro AUC `0.9144`
- Random forest:
  - accuracy `0.9128`
  - micro AUC `0.9916`
  - macro AUC `0.9719`
- ExtraTrees:
  - accuracy `0.9152`
  - micro AUC `0.9927`
  - macro AUC `0.9725`

Selected-model presentation:

- Genome-held-out selected SABR ROC figure highlights difficult boundary
  classes:
  - micro-average AUC `0.995`
  - `III-A` AUC `0.992`
  - `III-B` AUC `0.952`
  - `III-D` AUC `0.943`
  - `I-B` AUC `0.988`
  - `I-C` AUC `0.999`
- The normalized confusion heatmap and t-SNE subtype-focus panel reinforce
  that Type III-B/III-D and adjacent Type I classes remain the main unresolved
  independent-validation target.

Verification:

- Regenerated advanced figures and `docs/SABR_manuscript.docx`.
- Validated revised document/figure scripts with `python -m py_compile`.
- The earlier independent CCTyper pilot remains verified at `25/25` correct
  arrays and is explicitly presented as limited by the absence of `III-B` and
  `III-D` validation rows.

## Sensitivity/resistance prediction database expansion

New direction:

- Expand SABR from CRISPR targeting-evidence mapping toward bacteria-phage
  sensitivity/resistance and EOP prediction.
- Keep experimental phenotype/EOP observations separate from SABR evidence
  scores in `data/curation/phage_host_interactions.tsv`.
- Preserve anti-CRISPR and other defense evidence fields so model features can
  explain exceptions where spacer/PAM evidence does not imply resistance.

Current database foundation:

- Added `crispr_phage_predictor/interactions.py` with controlled schema,
  validation, EOP parsing, EOP class normalization, and susceptibility-label
  derivation.
- Added `tests/test_interactions.py` for schema and parser coverage.
- Added `scripts/normalize_phage_host_interactions.py`.
- Added `scripts/import_mirzaei_nilsson_2015_eop.py` for reproducible import
  of Mirzaei and Nilsson 2015 S1 DOCX EOP matrix.
- Current master table: `data/curation/phage_host_interactions.tsv`.
- Current imported source file: `data/curation/mirzaei_nilsson_2015_s1.docx`.
- Current derived source table: `data/curation/mirzaei_nilsson_2015_eop.tsv`.

Current interaction coverage:

- `562` total curated/measured interaction rows.
- `13` initial hand-curated rows:
  - PA14/JBD plaque-assay CRISPR-resistance rows from Cady et al. 2012.
  - DGCC7710/phage 2972 BIM EOP rows from Magadan et al. 2012 Table 3.
- `257` programmatically imported EOP rows from Mirzaei and Nilsson 2015 S1:
  - six phages `SU10`, `SU16`, `SU27`, `SU32`, `SU57`, `SU63`;
  - ECOR, ESBL E. coli, SARA, and SARB bacterial collections;
  - source cells with blanks were not treated as measured interactions.
- `292` programmatically generated rows from embedded open-literature tables:
  - `Sofy2021_KPP5`: 31 KPP-5 host-range/EOP-category rows.
  - `Latka2021_K23Klebsiella`: 24 K23 Klebsiella numeric EOP rows.
  - `Singh2025_KpCocktail`: 12 KPKp/KSKp/cocktail numeric EOP rows.
  - `Abdelhadi2021_STM2`: 108 Salmonella/non-Salmonella EOP-category rows.
  - `Wintachai2022_vABWU2101`: 22 Acinetobacter/control numeric EOP rows.
  - `Wintachai2024_vECPW8`: 20 MDR APEC numeric EOP rows.
  - `Kim2022_KP1_KP12`: 40 Klebsiella/non-Klebsiella percent-EOP rows,
    converted to ratios in the SABR schema.
  - `Ciacci2018_vB_Kpn_F48`: 14 K. pneumoniae numeric EOP rows.
  - `Peng2025_vB_Kp_XP4`: 21 K. pneumoniae percent-EOP host-range rows,
    converted to ratios in the SABR schema.

Next curation targets:

- Add source-ledger rows and import/curate additional EOP matrices from open
  supplementary tables in PLOS, Frontiers, BMC/Virology Journal, MDPI, and
  Microbiology Spectrum.
- Prioritize strain-level EOP matrices with explicit host strain names,
  phage names/accessions, reference host, and methods.
- Add anti-CRISPR annotation columns later by linking phage genomes through
  AcrFinder/PADLOC/DefenseFinder-style outputs where accessioned genomes are
  available.

Genome-linkage foundation:

- Added `crispr_phage_predictor/accession_linkage.py` with a controlled
  accession-linkage schema and coverage calculation.
- Added `tests/test_accession_linkage.py`.
- Added `data/curation/accession_linkage.tsv`.
- Added `scripts/report_accession_linkage_coverage.py`.
- Current linkage coverage for `562` interaction rows:
  - `145` rows have linked phage genome accessions.
  - `13` rows have linked bacterial genome accessions.
  - `13` rows are full pair-genome-ready.
- Hybrid mode enabled:
  - Added `reference_proxy` linkage status for public reference host genomes.
  - `562/562` rows now have either exact host genome linkage or a downloaded
    reference-proxy host genome.
  - `145/562` rows are hybrid-ready: exact phage genome plus exact/alias/proxy
    host genome.
  - `13/562` rows remain strict full pair-genome-ready.
- Added `scripts/download_accession_linkage_genomes.py`.
- Downloaded/linked public FASTA records under `data/curation/downloads/` and
  recorded them in `data/curation/downloaded_records.tsv`.
- New downloaded phage genomes include:
  - `JBD18`, `JBD25`, `JBD67`, `2972`, `vB_Kpn_F48`, `vB_Kp_XP4`, `KPP-5`,
    `KpS8`, `vB_KpnM_Seu621`, `vB_KpnP_Dlv622`, `vABWU2101`, and `vECPW8`.
- New downloaded host reference proxies include:
  - `E. coli K-12 MG1655` (`NC_000913.3`),
  - `Salmonella Typhimurium LT2` (`NC_003197.2`),
  - `K. pneumoniae HS11286` (`NC_016845.1`),
  - `P. aeruginosa PAO1` (`NC_002516.2`),
  - `S. aureus NCTC 8325` (`NC_007795.1`),
  - `A. baumannii ATCC 17978` (`NC_009085.1`),
  - `Citrobacter freundii` proxy (`NZ_CP007557.1`),
  - `Cronobacter sakazakii ATCC BAA-894` (`NC_009778.1`),
  - `Proteus mirabilis HI4320` (`NC_010554.1`).
- Full pair-ready rows currently come from:
  - `Cady2012_PA14_CRISPR`: PA14/JBD18/JBD25/JBD67.
  - `Magadan2012_CRISPR3`: DGCC7710/2972 BIM rows, using parent DGCC7710
    genome linkage for BIM derivatives.
- Newly resolved phage accessions include:
  - `vB_Kpn_F48`: `MG746602`.
  - `vB_Kp_XP4`: `PP663283`.
  - `KPP-5`: `MW600722.1`.
  - `KpS8`: `MT178275`.
  - `vB_KpnM_Seu621`: `MT939253`.
  - `vB_KpnP_Dlv622`: `MT939252`.
  - `vABWU2101`: `OK546191.1`.
  - `vECPW8`: `PQ362703`.
- Main remaining bottleneck: exact host-strain genome accessions for large
  phenotype panels, especially ECOR/ESBL/SARA/SARB, MDR Klebsiella clinical
  panels, MDR APEC `PW*` isolates, and A. baumannii `ABPW*` isolates.
- Remaining high-value unresolved phage groups:
  - `MirzaeiNilsson2015_EOP`: `SU10`, `SU16`, `SU27`, `SU32`, `SU57`, `SU63`.
  - `Abdelhadi2021_STM2`: `vB_STS_1`, `vB_STM_2`, `vB_STS_3`.
  - `Kim2022_KP1_KP12`: `KP1`, `KP12` still need direct accession
    confirmation.
  - `Singh2025_KpCocktail`: `KPKp`, `KSKp`, and cocktail rows need phage
    accession resolution or source-specific no-genome flagging.

Accession/linkage improvement pass:

- Resolved additional complete phage genome accessions:
  - `SU10`: `KM044272` (RefSeq alternate `NC_027395.1`).
  - `SU57`: `MT511058`.
  - `KSKp`: `PQ306550`.
- Updated `accession_linkage.tsv` and downloaded corresponding FASTA and
  GenBank records.
- Hybrid-ready rows increased from `145` to `213`.
- Current tier counts:
  - `tier1_exact_pair`: `13`.
  - `tier2_proxy_host_exact_phage`: `200`.
  - `tier4_host_only_or_proxy`: `349`.
- Added `dataset_tier` to `accession_linkage_coverage.tsv` and propagated it
  into generated feature tables as metadata. It is excluded from model features.
- Added `scripts/resolve_host_accession_candidates.py`.
- Generated `data/curation/host_accession_candidates.tsv`:
  - queried `80` unresolved host strain/source pairs;
  - produced `190` candidate rows;
  - `166` candidate rows require review;
  - many ECOR hits are plasmids, partial genes, or WGS master records, so no
    automatic exact host links were promoted.

First model-ready phage-host feature table:

- Added `crispr_phage_predictor/interaction_features.py`.
- Added `scripts/build_interaction_feature_table.py`.
- Added `tests/test_interaction_features.py`.
- Generated `data/training/phage_host_interaction_features.tsv`.
- Current table after phage-accession improvement:
  - `213` hybrid-ready interaction rows.
  - `155` columns.
  - strict exact/strain-alias host rows: `13`.
  - reference-proxy host rows: `200`.
  - binary susceptibility labels:
    - `resistant`: `77`.
    - `susceptible`: `68`.
  - EOP classes:
    - `high`: `43`.
    - `medium`: `15`.
    - `low`: `10`.
    - `trace`: `7`.
    - `none`: `67`.
    - `not_reported`: `3`.
- Feature families currently included:
  - phenotype labels and EOP class/value,
  - exact/proxy linkage provenance,
  - host/phage accessions and local FASTA paths,
  - host and phage genome length,
  - host and phage GC percentage,
  - host and phage 3-mer frequencies,
  - phage/host length ratio,
  - phage-host GC delta.
- This is an auditable baseline feature table only. It does not yet include
  CRISPR spacer/PAM matches, anti-CRISPR predictions, receptor features, or
  defense-system annotations.

First phage-host baseline classifier:

- Added `crispr_phage_predictor/phage_host_baseline.py`.
- Added `scripts/train_phage_host_baseline.py`.
- Added `tests/test_phage_host_baseline.py`.
- Generated outputs under `data/training/phage_host_baseline/`.
- Target: `binary_susceptibility` after excluding `eop_class = not_reported`.
- Metadata leakage controls:
  - excludes interaction IDs, source keys, organism names, strain names, phage
    names, accessions, local paths, labels, and EOP values from model features;
  - retains numeric genome/k-mer features and `uses_reference_proxy_host`.
- Evaluation summary:
  - row-random split:
    - majority baseline accuracy `0.5116`, macro F1 `0.3385`;
    - logistic regression accuracy `0.7674`, macro F1 `0.7643`;
    - random forest accuracy `0.7209`, macro F1 `0.7196`;
    - ExtraTrees accuracy `0.7209`, macro F1 `0.7171`.
  - phage-held-out split:
    - majority baseline accuracy `0.3774`, macro F1 `0.2740`;
    - logistic regression accuracy `0.6226`, macro F1 `0.3837`;
    - random forest/ExtraTrees accuracy `0.3774`, macro F1 `0.3437`.
  - source-held-out split:
    - majority baseline accuracy `0.4000`, macro F1 `0.2857`;
    - logistic/random forest/ExtraTrees accuracy `0.4833`, macro F1 `0.4762`.
- Interpretation:
  - row-random performance is optimistic and should not be used for manuscript
    claims;
  - grouped splits show weak generalization, especially with small, proxy-heavy
    data and no CRISPR/PAM/Acr/receptor/defense features yet;
  - next model step should add biologically meaningful features and resolve
    more exact phage/host accessions before optimizing algorithms.

First CRISPR targeting feature layer:

- Added `crispr_phage_predictor/interaction_targeting_features.py`.
- Added `scripts/add_targeting_features.py`.
- Added `tests/test_interaction_targeting_features.py`.
- Generated `data/training/phage_host_interaction_features_with_targeting.tsv`.
- The feature builder uses Diced/MinCED-compatible detection when available,
  exact spacer matching, and SABR's CRISPR targeting-evidence score.
- Internal exact-repeat fallback is limited to small FASTA records because the
  transparent internal scanner is too slow for multi-Mb reference genomes.
- Current augmented table after phage-accession improvement:
  - `213` rows.
  - `169` columns.
  - `13` rows with spacer hits.
  - max CRISPR targeting score `60.0`.
- Spacer-hit rows are exactly the strict CRISPR-positive systems:
  - PA14/JBD18/JBD25/JBD67 rows.
  - DGCC7710/2972 BIM rows.
- Re-running `scripts/train_phage_host_baseline.py` on the targeting-augmented
  table produced the same summary metrics as the composition-only table because
  targeting features are nonzero only in 13 rows and most hybrid rows use
  reference-proxy hosts with no source-specific spacer acquisition.
- Conclusion: targeting layer works technically and recovers known positives,
  but model improvement requires exact assayed host genomes, phage accessions
  for remaining panels, anti-CRISPR features, receptor/tail-fiber features, and
  broader CRISPR-positive/negative phenotype diversity.

Phage-side GenBank annotation feature layer:

- Added `crispr_phage_predictor/phage_annotation_features.py`.
- Added `scripts/download_phage_genbank_records.py`.
- Added `scripts/add_phage_annotation_features.py`.
- Added `tests/test_phage_annotation_features.py`.
- Downloaded `12` linked phage GenBank records to
  `data/curation/downloads/phages_genbank/`.
- Generated `data/training/phage_host_interaction_features_with_annotations.tsv`.
- Current annotation-augmented table after phage-accession improvement:
  - `213` rows.
  - `191` columns.
  - annotation features include CDS count, hypothetical CDS fraction, integrase,
    recombinase, excisionase, repressor, terminase, portal, capsid, tail,
    tail-fiber/receptor-binding, baseplate, holin, lysin/endolysin,
    depolymerase, DNA methyltransferase, restriction-evasion, anti-CRISPR
    keyword, temperate marker, structural marker, and lysis marker counts.
- Annotation results:
  - no anti-CRISPR keyword hits in the current linked GenBank product
    annotations;
  - temperate-marker keywords are present in several phages, including JBD
    phages, 2972, vB_KpnP_Dlv622, and vB_Kpn_F48;
  - tail/structural/lysis markers are broadly detected.
- Baseline retraining on the 213-row annotation-augmented table:
  - row-random logistic regression accuracy `0.7460`, macro F1 `0.7429`;
  - row-random random forest accuracy `0.7460`, macro F1 `0.7444`;
  - phage-held-out random forest accuracy `0.7000`, macro F1 `0.6875`;
  - phage-held-out ExtraTrees accuracy `0.7000`, macro F1 `0.7000`;
  - source-held-out logistic/random forest/ExtraTrees accuracy `0.4400`,
    macro F1 `0.4400`.
- Interpretation:
  - annotation features did not improve generalization;
  - this is expected because most annotation features are constant per phage
    while rows mix susceptible and resistant outcomes for the same phage;
  - next useful feature direction is interaction-specific host/phage matching:
    receptor/tail-fiber similarity, exact host genomes, defense-system calls,
    and curated Acr databases rather than raw keyword counts alone.

Exact host assembly promotion pass:

- Added `scripts/promote_host_assembly_links.py`.
- The script queries NCBI Assembly for unresolved host strain/source pairs,
  fetches assembly reports, and promotes only candidates whose report contains
  the exact expected strain name and `Genome representation: full`.
- The nucleotide candidate table remains useful for review, but Assembly is the
  correct source for whole-host FASTA files because many nucleotide hits are
  plasmids, partial genes, or WGS master placeholders.
- Dry run over the first 200 unresolved strain/source pairs:
  - `230` assembly candidate rows reviewed.
  - `71` exact assemblies passed report-level checks.
- Committed run:
  - appended `71` exact bacterial assembly rows to
    `data/curation/accession_linkage.tsv`.
  - downloaded their FASTA files under `data/curation/downloads/bacteria/`.
  - wrote candidate review evidence to
    `data/curation/host_assembly_candidates.tsv`.
- Current linkage coverage after promotion:
  - total interaction rows: `562`.
  - exact/strain-alias bacterial genome linked rows: `197`.
  - exact/proxy bacterial genome linked rows: `562`.
  - phage genome linked rows: `213`.
  - strict exact-pair genome-ready rows: `74` (up from `13`).
  - hybrid-ready rows: `213`.
  - dataset tiers:
    - `tier1_exact_pair`: `74`.
    - `tier2_proxy_host_exact_phage`: `139`.
    - `tier4_host_only_or_proxy`: `349`.
- Rebuilt downstream outputs:
  - `data/curation/accession_linkage_coverage.tsv`.
  - `data/training/phage_host_interaction_features.tsv`.
  - `data/training/phage_host_interaction_features_with_targeting.tsv`.
  - `data/training/phage_host_interaction_features_with_annotations.tsv`.
  - `data/training/phage_host_baseline_with_annotations/`.
- Current modelable table remains `213` rows because exact phage genome linkage
  is now the limiting factor, not host genome availability for the promoted
  ECOR rows.
- Baseline metrics on the combined targeting + annotation feature table are
  unchanged from the prior annotation run:
  - row-random logistic/random-forest/ExtraTrees accuracy about `0.746`.
  - phage-held-out random forest and ExtraTrees accuracy `0.700`.
  - source-held-out models accuracy `0.440`.
- Interpretation:
  - this pass mainly improves provenance and enables true exact-host feature
    extraction for many Mirzaei/Nilsson ECOR rows;
  - it does not add new modelable rows until more phages are resolved;
  - the next bottleneck is phage accession resolution for SU16, SU27, SU32,
    SU63 and other unresolved panels, followed by host defense-system and
    receptor/tail-fiber compatibility features.

Database QA audit:

- Added `scripts/audit_phage_host_database.py`.
- The audit checks:
  - interaction and accession-linkage schemas;
  - duplicate IDs and duplicate exact linkage targets;
  - EOP class/value consistency;
  - susceptibility-label consistency warnings;
  - accession-linkage local FASTA paths;
  - downloaded FASTA record counts and base-pair totals;
  - coverage table alignment with interaction IDs;
  - feature table alignment with hybrid-ready coverage rows;
  - promoted host assembly evidence.
- First stricter run found one stale metadata issue in
  `data/curation/downloaded_records.tsv`:
  - `NC_002737.2` was recorded as `1852441` bp;
  - the local FASTA contains `1852433` bp.
- Corrected the `NC_002737.2` downloaded-record total to match the FASTA.
- Current QA status:
  - `python scripts/audit_phage_host_database.py` reports `audit_status clean`
    and `issues 0`.
  - `python -m py_compile scripts/audit_phage_host_database.py scripts/promote_host_assembly_links.py`
    passes.
  - focused phage-host tests report `19 passed`.

Fuzzy CRISPR targeting and threshold experiment:

- Added fuzzy spacer/protospacer matching to
  `crispr_phage_predictor/interaction_targeting_features.py`.
- New targeting variables include:
  - `fuzzy_spacer_candidate_count`;
  - `fuzzy_unique_spacers`;
  - `fuzzy_targeted_phage_contig_count`;
  - `best_fuzzy_spacer_identity`;
  - `best_fuzzy_spacer_mismatches`;
  - `best_fuzzy_seed_edge_mismatches`;
  - `best_fuzzy_distal_mismatches`;
  - `fuzzy_seed_perfect_hit_count`;
  - `fuzzy_near_perfect_hit_count`;
  - `fuzzy_high_confidence_hit_count`;
  - `graded_crispr_interference_score`.
- The fuzzy scan allows imperfect spacer hits up to four mismatches and uses the
  better of the two spacer ends as a conservative possible PAM-proximal seed
  signal when subtype/PAM orientation is not available.
- Added `scripts/evaluate_targeting_thresholds.py` to evaluate CRISPR targeting
  thresholds directly against resistance/susceptibility labels.
- Added `--restrict-tier` to `scripts/train_phage_host_baseline.py` so exact
  host-phage rows can be evaluated separately from proxy-host rows.
- Regenerated:
  - `data/training/phage_host_interaction_features_with_targeting.tsv`
    (`213` rows, `180` columns).
  - `data/training/phage_host_interaction_features_with_annotations.tsv`
    (`213` rows, `202` columns).
  - `data/training/phage_host_baseline_with_fuzzy_targeting/`.
  - `data/training/phage_host_baseline_with_fuzzy_targeting_tier1_exact_pair/`.
  - `data/training/targeting_threshold_evaluation.tsv`.
  - `data/training/targeting_threshold_evaluation_tier1_exact_pair.tsv`.
- Current targeting coverage:
  - exact spacer-hit rows remain `13`.
  - fuzzy spacer-candidate rows are also `13`.
  - fuzzy high-confidence rows are `13`.
  - this means the improved alignment logic works, but current CRISPR-detected
    host arrays still only explain the PA14/JBD and DGCC7710/2972 resistance
    systems.
- Threshold interpretation:
  - On all modelable rows, CRISPR targeting thresholding reaches accuracy
    `0.542857`, macro F1 `0.428312`.
  - On tier1 exact-pair rows, CRISPR targeting thresholding reaches accuracy
    `0.704225`, macro F1 `0.639942`.
  - In the tier1 exact-pair subset it predicts `10` resistant rows and has `0`
    susceptible false positives, but misses `21` resistant rows.
  - Therefore CRISPR targeting is currently high-precision but low-recall for
    resistance in this database.
- Baseline retraining with fuzzy variables:
  - all modelable rows:
    - row-random random forest accuracy `0.714286`, macro F1 `0.713636`;
    - phage-held-out random forest/ExtraTrees accuracy `0.700000`, macro F1
      about `0.700000`;
    - source-held-out models accuracy `0.440000`, macro F1 `0.440000`.
  - tier1 exact-pair only:
    - row-random ExtraTrees accuracy `0.818182`, macro F1 `0.811966`;
    - grouped splits are unstable/collapsed because the exact-pair subset has
      only `71` modelable rows after excluding not-reported EOP rows and is
      source/phage biased.
- Scientific interpretation:
  - imperfect matching features are now represented, which addresses the concern
    that resistance can occur without a perfect spacer/protospacer match;
  - current data do not yet show additional fuzzy-only positives beyond exact
    CRISPR systems;
  - next real model improvement requires more CRISPR-positive assayed host
    genomes, calibrated subtype-specific PAM rules, and Acr/defense-system
    homology features.
- Current validation:
  - `python scripts/audit_phage_host_database.py` reports clean with `0` issues.
  - focused phage-host tests report `20 passed`.

Release-hardening pass for v0.1-beta:

- Added explicit README positioning: SABR complements existing CRISPR array
  callers, Cas-locus typers, spacer-target search tools, and phage-host
  predictors by integrating uploaded bacterial/phage FASTA panel analysis into
  a bacteria-by-phage CRISPR targeting-evidence matrix with repeat-derived
  subtype and PAM/PFS support.
- Updated `RELEASE_CHECKLIST.md` so Docker support, model SHA-256 checking, and
  tiny synthetic example FASTAs are listed as ready rather than future tasks.
- Updated `DEPLOYMENT.md` to use a clean-checkout-safe pytest temp directory:
  `python -m pytest --basetemp=.pytest-tmp -p no:cacheprovider`.
- Verification on 2026-06-09:
  - `python -m compileall -q app.py crispr_phage_predictor scripts docs` passed.
  - `python scripts/ensure_model_artifact.py` found the local runtime artifact.
  - focused release tests passed: `11 passed`.
  - full pytest suite passed with cache disabled: `106 passed`.
  - `python scripts/audit_phage_host_database.py` reported `audit_status clean`
    and `issues 0`.
  - foreground Streamlit startup reached the serving state on
    `http://localhost:8502`; hidden background launch exited in this shell
    environment, so the persistent local server was not left running.

UI update after release-hardening:

- The Streamlit CRISPR detection selector now defaults to
  `MinCED-compatible` instead of `Auto recommended`.
- Removed the sidebar warning phrase that described MinCED-compatible detection
  as benchmarking-only and potentially slow on whole-genome uploads.
- Renamed the primary button to `Run SABR analysis`.
- Analysis progress bars now render in the sidebar directly under the run button
  rather than lower on the main results page, so users see immediate feedback
  while detection and spacer-phage matching are running.
- Fixed result/detail navigation after a completed run:
  - the app now keeps rendering the latest heatmap/results on ordinary widget
    reruns instead of stopping at the pre-run upload message;
  - `Open hit details` uses session-state navigation instead of normal query
    parameter mutation, avoiding half-populated rerun states;
  - if session state is missing, the app reloads the latest saved run from
    `outputs/runs/` using `evidence_matrix.csv`, `heatmap.csv`, and
    `spacer_hits.csv`.
- Fixed brittle PAM/PFS reporting:
  - previous logic only checked the immediately adjacent motif-length window,
    so a `3prime:NGG` rule on `TTCGAGGGTA` reported adjacent `TTC` as
    unsupported even though nearby `AGG` was visible in the displayed flank;
  - PAM/PFS evaluation now prefers the adjacent motif, then scans the displayed
    flank for the nearest compatible motif and reports `compatible_nearby`;
  - spacer-hit outputs and the detail page now include
    `pam_offset_from_protospacer` so users can see whether support is adjacent
    (`0`) or offset within the flank.
- Fixed BLAST hit-detail alignment display:
  - BLAST identity is computed over the aligned BLAST segment, but the detail
    page previously displayed the full original spacer against BLAST's subject
    segment, which made partial approximate hits look inconsistent;
  - BLAST output now includes `qseq` and `sseq`; spacer-hit exports preserve
    `aligned_spacer_sequence` and `aligned_protospacer_sequence`;
  - the hit-details alignment now renders those aligned BLAST strings when
    available, while keeping the full original spacer for provenance.
- Article-linked release direction:
  - use GitHub as the development repository;
  - archive the frozen `v0.1-beta` source release plus runtime model artifact on
    Zenodo for DOI-backed citation;
  - deploy the live browser demo on Hugging Face Spaces with Docker;
  - configure `SABR_MODEL_URL` to the Zenodo-hosted model URL and
    `SABR_MODEL_SHA256` to
    `ae0e8a5d56d6a2a4eb8206fb916fd9dee51f7fb9276528346b0cb8279b76cd32`.
- Added publication/release metadata files:
  - `CITATION.cff`;
  - `.zenodo.json`;
  - `MODEL_CARD.md`;
  - `RELEASE_RUNBOOK.md`.
- Added an in-app public-demo privacy note warning users not to upload sensitive
  unpublished genomes to hosted demos and recommending local Docker for private
  lab analyses.

Current publication/deployment state as of 2026-06-09:

- Recommended public/article release architecture is settled:
  - GitHub for active development and issue tracking.
  - Zenodo for the frozen citable `v0.1-beta` source/model release DOI.
  - Hugging Face Spaces for the live Docker/Streamlit browser demo.
- Release metadata files now exist:
  - `CITATION.cff`;
  - `.zenodo.json`;
  - `MODEL_CARD.md`;
  - `RELEASE_RUNBOOK.md`;
  - `deployment/huggingface/README.md`.
- Runtime model artifact selected for release:
  - path: `models/cas_subtype_extratrees.joblib`;
  - size: `410202588` bytes;
  - SHA-256:
    `ae0e8a5d56d6a2a4eb8206fb916fd9dee51f7fb9276528346b0cb8279b76cd32`;
  - method: ExtraTrees;
  - training rows: `23478`;
  - classes: `I-A`, `I-B`, `I-C`, `I-D`, `I-E`, `I-F`, `I-G`,
    `II-A`, `II-B`, `II-C`, `III-A`, `III-B`, `III-C`, `III-D`,
    `V-A`, `V-K`, `VI-B1`.
- Account/release metadata status:
  - GitHub repository URL is set to `https://github.com/ESA-zoph/SABR`;
  - Hugging Face Space URL is set to `https://huggingface.co/spaces/EsberSaba/SABR`;
  - license identifier is set to `MIT`;
  - Zenodo DOI is set to `10.5281/zenodo.20737961`;
  - DOI-backed `SABR_MODEL_URL` is set to
    `https://zenodo.org/records/20737961/files/cas_subtype_extratrees.joblib?download=1`.
- Latest release-prep verification:
  - `python -m py_compile app.py` passed;
  - focused tests for BLAST alignment/PAM/output paths passed: `21 passed`;
  - local Streamlit app restarted and returned HTTP `200 OK` at
    `http://127.0.0.1:8502`.
- Current live local app process was restarted after adding release-prep UI/docs
  changes. Refresh the browser before testing.

Current functionality/update state as of 2026-06-11:

- Product direction remains a two-track plan:
  - keep the current app stabilizable as a `v0.1-beta` public evidence-mapping
    release;
  - continue adding functionality on top without changing the claim to
    phenotype-level resistance prediction.
- Added first pass of post-beta functionality focused on usability and
  interpretation rather than performance:
  - sidebar switch to load the built-in demo FASTA panel;
  - saved-run selector for reopening previous `outputs/runs/` analyses;
  - highest-evidence pair summary with no/weak/moderate/strong targeting
    evidence labels;
  - heatmap-cell detail page now shows the score label, plain-language
    interpretation, and score-component table alongside spacer/protospacer,
    PAM/PFS, and seed details;
  - optional known-interaction CSV/TSV upload with `bacterium` and `phage`
    columns, plus optional `expected_label` or `susceptibility_label`, compared
    against SABR targeting evidence;
  - input-quality warning table for fragmented uploads, high ambiguous-base
    content, and many short records;
  - downloadable Markdown report for current/saved runs, also written to each
    run folder as `sabr_report.md` when possible.
- Updated `README.md` current-status bullets to include the new demo,
  interpretation, saved-run, known-interaction comparison, and report features.
- Verification for this pass:
  - `.venv\Scripts\python.exe -m py_compile app.py` passed.
  - focused app/output tests passed: `14 passed`.
  - full pytest suite passed with cache disabled and `.pytest-tmp` basetemp:
    `109 passed`.
- Still recommended before public beta:
  - manually smoke-test the built-in demo through Streamlit;
  - Docker-build and browser-test the same demo;
  - finalize repository/Space/Zenodo/license placeholders before tagging.

UI simplification decision on 2026-06-11:

- Removed the public-facing saved-run selector from the sidebar because it made
  the first-use workflow feel more complicated than necessary.
- Removed the public-facing known-interaction table upload/comparison from the
  Streamlit UI for the same reason.
- Kept automatic timestamped run outputs under `outputs/runs/` and kept the
  Markdown report download/write path, because those support reproducibility
  without adding extra setup decisions for a first-time user.
- Current intended public workflow is now:
  - load the built-in demo or upload bacterial/phage FASTA files;
  - choose methods only if needed;
  - run SABR analysis;
  - inspect heatmap, highest-evidence pairs, hit details, score components, and
    download the report.

Built-in demo update on 2026-06-11:

- Replaced the app's built-in demo target with a lightweight real-data panel
  under `data/examples/real_demo/`.
- The demo bacterial input is `real_demo/bacteria/PA14_CRISPR_region.fasta`, a
  PA14 chromosome excerpt from the existing full PA14 example genome spanning
  `NC_008463.1:2926400-2937300` and covering the two PA14 CRISPR arrays.
- The demo phage inputs are:
  - `real_demo/phages/JBD18_positive_control.fasta`, a real JBD18 phage genome
    expected to show PA14 spacer-targeting evidence;
  - `real_demo/phages/Lambda_negative_control.fasta`, an unrelated real phage
    genome expected to show no PA14 spacer-targeting evidence.
- Internal-backend validation of the real demo produced:
  - `2` candidate arrays;
  - `33` extracted spacers;
  - `4` spacer hits;
  - PA14 CRISPR region vs JBD18 after reverse-strand seed/PAM display fix:
    `4` unique spacers, `2/4` PAM-supported hits, best seed mismatches `0`,
    targeting score `99.97`, strong candidate CRISPR targeting evidence when
    the subtype model artifact is available;
  - PA14 CRISPR region vs Lambda: `0` unique spacers, score `0.0`, no
    spacer-match evidence.

UI interpretation clarification on 2026-06-12:

- Added a plain-language note under each selected pair's evidence label:
  - if PAM/PFS was not evaluated, the app now states that the score is driven by
    spacer count plus identity/coverage and that PAM/PFS/seed layers did not
    contribute;
  - if PAM/PFS was evaluated, the app states whether evaluated hits supported
    the selected rule.
- Added a per-hit warning when `PAM/PFS support` is `not_evaluated` because no
  PAM/PFS rule was selected for that hit. The warning points users to
  `Auto from predicted subtype` or a manual expert override before running the
  analysis.
- Interpretation decision: a score such as `60.00` can be moderate candidate
  CRISPR targeting evidence from multiple high-identity spacer/protospacer hits
  even when PAM/PFS was not evaluated. It should be read as spacer evidence, not
  PAM-supported evidence, unless the PAM/PFS fields are populated.

Demo-default polish on 2026-06-12:

- Changed first-run defaults to the internal exact-repeat CRISPR detector and
  internal exact spacer matcher. This makes the built-in real demo portable and
  predictable even when MinCED/BLAST are not installed.
- Kept PAM/PFS evaluation defaulted to `Auto from predicted subtype`.
- Added an in-app expected-result note when `Use built-in real demo` is checked:
  PA14/JBD18 should show spacer-targeting evidence, while PA14/Lambda should
  show no spacer-match evidence.

Reverse-strand hit display fix on 2026-06-12:

- Hit details previously displayed the raw phage genomic protospacer segment in
  the Spacer/Protospacer alignment block. For reverse-strand hits, this made a
  true 100% reverse-complement match appear mismatched on screen.
- The detail page now reverse-complements internal exact-match protospacers for
  display when `strand == "-"`, so the visual alignment agrees with the
  identity value.
- The phage-reference context block remains in raw genomic orientation and now
  has a caption explaining the distinction for reverse-strand hits.
- The PAM/seed annotation path now also uses a strand-oriented protospacer for
  seed-mismatch summaries. Before this fix, reverse-strand seed mismatch counts
  could be computed against the raw genomic segment instead of the matching
  protospacer orientation.

Saved-output hardening on 2026-06-12:

- Empty summary tables now preserve CSV headers for candidate arrays, spacers,
  spacer hits, PAM subtype support, and evidence matrices.
- This prevents no-hit/no-array runs from writing zero-column CSV files that
  later trigger `pandas.errors.EmptyDataError` when the app reloads saved
  results or reports.
- Removed silent auto-loading of the latest saved run from disk in the main app
  view. This prevents stale old runs from appearing after a refresh and
  confusing users about whether they are looking at the newly selected demo or
  an older analysis. Saved outputs remain on disk, and hit-details recovery can
  still load a specific run ID when one is present.
- Verification after this hardening pass:
  - `.venv\Scripts\python.exe -m compileall -q app.py crispr_phage_predictor tests`
    passed;
  - `.venv\Scripts\python.exe -m pytest --basetemp=.pytest-tmp -p no:cacheprovider -q`
    passed with `111 passed`;
  - local Streamlit endpoint returned HTTP `200 OK` at
    `http://127.0.0.1:8501`.

Release metadata update on 2026-06-12:

- User selected MIT as the public-beta open-source license.
- Updated `CITATION.cff` and `.zenodo.json` license fields to `MIT`.
- Added repository `LICENSE` with MIT text and copyright holder:
  `Copyright (c) 2026 Esber Saba, The Phage Lab, Faculty of Medicine, American University of Beirut`.
- Updated `RELEASE_CHECKLIST.md` and `RELEASE_RUNBOOK.md` so the remaining
  license task is only to confirm the copyright holder text before tagging.
- User provided public support/contact details:
  - Esber Saba, Ph.D.
  - The Phage Lab, Faculty of Medicine, American University of Beirut
  - `es60@aub.edu.lb`
- Added this contact to `CITATION.cff`, `.zenodo.json`, `README.md`,
  `DEPLOYMENT.md`, `deployment/huggingface/README.md`,
  `RELEASE_CHECKLIST.md`, and `RELEASE_RUNBOOK.md`.

Release-readiness checkpoint on 2026-06-17:

- Verification passed:
  - `.venv\Scripts\python.exe -m compileall -q app.py crispr_phage_predictor scripts docs`;
  - `.venv\Scripts\python.exe scripts\ensure_model_artifact.py`;
  - `.venv\Scripts\python.exe scripts\audit_phage_host_database.py`
    reported `audit_status clean` and `issues 0`;
  - `.venv\Scripts\python.exe -m pytest --basetemp=.pytest-tmp -p no:cacheprovider -q`
    passed with `111 passed` and two expected Biopython warnings from malformed
    synthetic GenBank test fixtures.
- Foreground Streamlit startup reached the serving state on
  `http://localhost:8501`; detached background launch still exits in this shell
  environment, so browser testing should use the foreground command or Docker.
- Built-in real demo pipeline validation with internal detector/matcher and
  automatic subtype-derived PAM/PFS annotation produced:
  - `2` candidate arrays;
  - `33` extracted spacers;
  - `4` spacer hits;
  - PA14 CRISPR region vs JBD18: score `99.97`, `4` unique matching spacers,
    `2/4` PAM-supported hits, strong candidate CRISPR targeting evidence;
  - PA14 CRISPR region vs Lambda: score `0.00`, no spacer-match evidence.
- Docker Desktop initially needed time to expose the Linux engine, then worked.
- Docker build passed for `sabr:local`.
- Updated the Dockerfile `CMD` to JSON shell form to avoid Docker's
  `JSONArgsRecommended` warning while preserving environment-variable expansion.
- Fresh Docker container `sabr-local-test` returned HTTP `200` at
  `http://127.0.0.1:7860`.
- Docker startup without `SABR_MODEL_URL` reports the documented missing-model
  fallback and still serves the app.
- Docker model-mounted browser test:
  - the first Docker test container still showed the sidebar warning
    `ExtraTrees subtype model artifact is missing` because it was the old
    container created before mounting the model;
  - recreated `sabr-local-test` with the local artifact mounted to
    `/home/user/app/models/cas_subtype_extratrees.joblib`;
  - verified inside the container that the model file exists at about `392M`;
  - Docker startup log changed to
    `Using model artifact at models/cas_subtype_extratrees.joblib.`;
  - HTTP check still returned `200` at `http://127.0.0.1:7860`;
  - browser refresh confirmed the subtype-model warning was gone.
- External public-release status:
  - Hugging Face Space URL:
    `https://huggingface.co/spaces/EsberSaba/SABR`;
  - Zenodo DOI:
    `10.5281/zenodo.20737961`
    (`https://doi.org/10.5281/zenodo.20737961`);
  - DOI-backed `SABR_MODEL_URL`:
    `https://zenodo.org/records/20737961/files/cas_subtype_extratrees.joblib?download=1`;
  - optional fresh-venv manual Streamlit smoke test remains useful before
    tagging a final GitHub release.

Hugging Face deployment checkpoint on 2026-06-18:

- Pushing the full GitHub repository history to Hugging Face was rejected
  because the history contains binary manuscript/research assets such as
  `.docx`, `.png`, and `.tif` files.
- Created a clean Space-only source tree under
  `release/hf-space-source/` with:
  - runtime app source;
  - `crispr_phage_predictor/`;
  - small `basic_demo` and `real_demo` FASTA inputs;
  - Dockerfile and requirements;
  - public metadata files;
  - no model artifact, no manuscript assets, no generated figures, and no
    research/provenance binaries.
- Removed optional Hugging Face YAML `emoji` metadata because the original
  `emoji: DNA` value was rejected and PowerShell encoded the emoji character as
  `??`.
- Clean Space source was force-pushed successfully to
  `https://huggingface.co/spaces/EsberSaba/SABR`.
- The first Space runtime attempt failed because
  `SABR_MODEL_URL` returned `HTTP Error 404: NOT FOUND`; cause was that the
  Zenodo model file was not yet publicly downloadable.
- After the Zenodo record/model file became public, the Hugging Face Space
  started successfully and the app is working.
- Required Space variables:
  - `SABR_MODEL_URL=https://zenodo.org/records/20737961/files/cas_subtype_extratrees.joblib?download=1`
  - `SABR_MODEL_SHA256=ae0e8a5d56d6a2a4eb8206fb916fd9dee51f7fb9276528346b0cb8279b76cd32`
- Main GitHub repo also received the matching Space README template fix:
  removed invalid optional `emoji` metadata from
  `deployment/huggingface/README.md`.
- Logo handling after Hugging Face binary rejection:
  - The original PNG logo remains supported locally as `assets/aub-logo.png`.
  - Added text sidecar support for `assets/aub-logo.png.b64`; when the PNG is
    absent, the app renders the same PNG from the base64 text file.
  - The Hugging Face Space source can therefore include the logo without pushing
    binary image files.

End-of-day handoff on 2026-06-18:

- Main GitHub repository is clean except for untracked local `release/`.
- Main GitHub repository has been pushed through commit:
  `a177265 Support text-encoded logo for Space deployment`.
- Temporary Hugging Face source repository at `release/hf-space-source/` is
  clean locally and includes commit:
  `01623dc Restore logo without binary asset`.
- If not already done before stopping, push the temporary Space source:

```powershell
cd C:\Users\es60\Desktop\Codex\release\hf-space-source
git push https://huggingface.co/spaces/EsberSaba/SABR master:main --force
```

- Tomorrow's next checks:
  1. Open `https://huggingface.co/spaces/EsberSaba/SABR`.
  2. Confirm the logo is visible.
  3. Confirm there is no subtype-model warning.
  4. Run the built-in real demo.
  5. Confirm PA14/JBD18 shows strong candidate CRISPR targeting evidence.
  6. Confirm PA14/Lambda shows no spacer-match evidence.
  7. If all checks pass, make the Space public if desired.
  8. Tag GitHub:

```powershell
cd C:\Users\es60\Desktop\Codex
git tag v0.1-beta
git push origin v0.1-beta
```
