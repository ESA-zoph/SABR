from __future__ import annotations

import argparse
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from crispr_phage_predictor.interaction_targeting_features import add_targeting_features


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Add SABR CRISPR spacer-targeting features to interaction features."
    )
    parser.add_argument(
        "--features",
        type=Path,
        default=Path("data/training/phage_host_interaction_features.tsv"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/training/phage_host_interaction_features_with_targeting.tsv"),
    )
    args = parser.parse_args()

    table = pd.read_csv(args.features, sep="\t")
    augmented = add_targeting_features(table)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    augmented.to_csv(args.output, sep="\t", index=False)
    print(f"rows\t{len(augmented)}")
    print(f"columns\t{len(augmented.columns)}")
    print("targeting_score_summary")
    print(augmented["crispr_targeting_score"].describe().to_string())
    print("rows_with_spacer_hits\t" + str(int((augmented["spacer_hit_count"] > 0).sum())))
    print("graded_crispr_interference_score_summary")
    print(augmented["graded_crispr_interference_score"].describe().to_string())
    print(
        "rows_with_fuzzy_spacer_candidates\t"
        + str(int((augmented["fuzzy_spacer_candidate_count"] > 0).sum()))
    )
    print(
        "rows_with_fuzzy_high_confidence_hits\t"
        + str(int((augmented["fuzzy_high_confidence_hit_count"] > 0).sum()))
    )


if __name__ == "__main__":
    main()
