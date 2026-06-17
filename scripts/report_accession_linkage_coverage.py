from __future__ import annotations

import argparse
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from crispr_phage_predictor.accession_linkage import (
    accession_coverage,
    load_accession_linkage_table,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Report genome-linkage coverage for SABR phage-host interactions."
    )
    parser.add_argument(
        "--interactions",
        type=Path,
        default=Path("data/curation/phage_host_interactions.tsv"),
    )
    parser.add_argument(
        "--linkage",
        type=Path,
        default=Path("data/curation/accession_linkage.tsv"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/curation/accession_linkage_coverage.tsv"),
    )
    args = parser.parse_args()

    interactions = pd.read_csv(args.interactions, sep="\t", dtype=str).fillna("")
    linkage = load_accession_linkage_table(args.linkage)
    coverage = accession_coverage(interactions, linkage)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    coverage.to_csv(args.output, sep="\t", index=False)
    _print_summary(coverage)


def _print_summary(coverage: pd.DataFrame) -> None:
    total = len(coverage)
    bacterium_ready = int(coverage["bacterium_genome_linked"].sum())
    bacterium_hybrid_ready = int(coverage["bacterium_reference_or_genome_linked"].sum())
    phage_ready = int(coverage["phage_genome_linked"].sum())
    pair_ready = int(coverage["pair_genome_ready"].sum())
    pair_hybrid_ready = int(coverage["pair_hybrid_ready"].sum())
    print(f"rows\t{total}")
    print(f"bacterium_genome_linked\t{bacterium_ready}")
    print(f"bacterium_reference_or_genome_linked\t{bacterium_hybrid_ready}")
    print(f"phage_genome_linked\t{phage_ready}")
    print(f"pair_genome_ready\t{pair_ready}")
    print(f"pair_hybrid_ready\t{pair_hybrid_ready}")
    print("pair_ready_by_source")
    print(
        coverage.groupby("source_key")["pair_genome_ready"]
        .agg(["count", "sum"])
        .to_string()
    )
    print("pair_hybrid_ready_by_source")
    print(
        coverage.groupby("source_key")["pair_hybrid_ready"]
        .agg(["count", "sum"])
        .to_string()
    )
    if "dataset_tier" in coverage.columns:
        print("dataset_tier")
        print(coverage.groupby("dataset_tier").size().to_string())


if __name__ == "__main__":
    main()
