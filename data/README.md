# Data Organization

This directory is intended for local example and training data. Large datasets should not be committed to version control without a clear reason.

Suggested local layout:

```text
data/
  examples/
    bacteria/
    phages/
  training/
    cas_repeat_labels.csv
  outputs/
```

## Example Inputs
Place bacterial genome FASTA files in `data/examples/bacteria/` and phage genome FASTA files in `data/examples/phages/`.

## Cas-Type Training Dataset
The first curated classifier dataset should use a table like:

```csv
genome_id,organism,contig_id,array_start,array_end,repeat_sequence,repeat_length,spacer_count,mean_spacer_length,cas_type,cas_subtype,pam,source
```

Required fields for the first classifier:

- `repeat_sequence`
- `cas_type`
- `source`

Useful optional fields:

- `cas_subtype`
- `pam`
- `organism`
- `spacer_count`
- `mean_spacer_length`
- `array_start`
- `array_end`

## Curation Notes
Keep notes on where each label came from. For a publishable tool, every training label should be traceable to a database, software output, or paper.
