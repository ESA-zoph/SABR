from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from crispr_phage_predictor.interaction_features import (
    build_hybrid_interaction_feature_table,
    load_inputs,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build a model-ready feature table for hybrid-ready SABR interactions."
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
        "--coverage",
        type=Path,
        default=Path("data/curation/accession_linkage_coverage.tsv"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/training/phage_host_interaction_features.tsv"),
    )
    parser.add_argument("--kmer-size", type=int, default=3)
    args = parser.parse_args()

    interactions, linkage, coverage = load_inputs(
        args.interactions,
        args.linkage,
        args.coverage,
    )
    features = build_hybrid_interaction_feature_table(
        interactions,
        linkage,
        coverage,
        k=args.kmer_size,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    features.to_csv(args.output, sep="\t", index=False)
    print(f"rows\t{len(features)}")
    print(f"columns\t{len(features.columns)}")
    if not features.empty:
        print("eop_class")
        print(features.groupby("eop_class").size().to_string())
        print("source_key")
        print(features.groupby("source_key").size().to_string())


if __name__ == "__main__":
    main()
