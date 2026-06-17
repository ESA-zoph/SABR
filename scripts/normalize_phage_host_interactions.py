from __future__ import annotations

import argparse
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from crispr_phage_predictor.interactions import (
    normalize_interaction_table,
    validate_interaction_table,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Normalize SABR phage-host interaction EOP curation tables."
    )
    parser.add_argument("input_tsv", type=Path)
    parser.add_argument("output_tsv", type=Path)
    args = parser.parse_args()

    table = pd.read_csv(args.input_tsv, sep="\t", dtype=str).fillna("")
    normalized = normalize_interaction_table(table)
    validate_interaction_table(normalized)
    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    normalized.to_csv(args.output_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
