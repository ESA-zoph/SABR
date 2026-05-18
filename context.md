# Project Context

## Working Name
CRISPR-Phage Resistance Predictor

## Goal
Build an easy-to-deploy bioinformatics tool for the scientific community that estimates hypothetical CRISPR-mediated bacterial resistance against bacteriophages.

The tool should accept multiple bacterial FASTA files and multiple phage FASTA files, identify CRISPR arrays, extract spacers and repeats, compare spacers against phage genomes, infer possible PAM compatibility, and produce an evidence-based bacteria-by-phage resistance likelihood matrix.

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

- predicted CRISPR-mediated resistance likelihood
- hypothetical resistance score
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

- predicted resistance likelihood score
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
- nearby cas gene profile when available
- source database labels

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

1. Bootstrap from CRISPRCasTyper/CCTyper-derived outputs because they expose consensus repeats, repeat length, spacer statistics, subtype predictions, subtype probabilities, trusted-array flags, and CRISPR-near-Cas tables.
2. Use nearby Cas-operon-supported subtype calls as the first training labels.
3. Use repeat-based predictions as a comparison baseline, not as ground truth labels.
4. Exclude ambiguous, hybrid, orphan, low-confidence, and incomplete systems from the first training set.
5. Later, build a reproducible RefSeq/GenBank complete-genome annotation pipeline using fixed tool versions.

The first local schema is tracked in `crispr_phage_predictor/ml/dataset.py`, with documentation in `data/training/README.md` and technical notes in `handout.md`.

Implemented dataset scaffolding now includes:

- repeat/Cas training table schema
- schema validation
- high-confidence label filtering
- CCTyper `crisprs_near_cas.tab`-style importer into the local schema
- tests for schema validation and CCTyper conversion

The first model should use repeat and array features only:

- repeat sequence
- repeat length
- repeat GC percent
- repeat k-mer counts
- spacer count
- mean spacer length
- spacer length variability where available

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
