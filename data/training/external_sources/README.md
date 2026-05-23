# External Training Sources

This folder is for raw third-party datasets used to build or audit SABR
training tables.

Raw archives, database dumps, and downloaded source files should stay out of
git. Keep a small README in each source-specific subfolder that records:

- source name and release/version
- original files present locally
- approximate record counts
- expected SABR use
- caveats before using the data for model claims

Derived, filtered, or curated tables should be written to `data/training/` with
provenance fields that point back to the source.
