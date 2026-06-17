# SABR Example Inputs

The Streamlit built-in demo uses `real_demo/`, a small panel derived from real
FASTA files already present in this repository:

- `real_demo/bacteria/PA14_CRISPR_region.fasta`: PA14 chromosome excerpt covering
  the two CRISPR arrays previously detected in the full PA14 example genome.
- `real_demo/phages/JBD18_positive_control.fasta`: real JBD18 phage genome,
  expected to show PA14 spacer-targeting evidence.
- `real_demo/phages/Lambda_negative_control.fasta`: unrelated real phage genome,
  expected to show no PA14 spacer-targeting evidence.

The older `basic_demo/` files are tiny synthetic smoke-test inputs. They remain
useful for parser and export testing, but they are not biological validation
data and should not be used for scientific claims.
