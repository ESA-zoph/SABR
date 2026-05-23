from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from crispr_phage_predictor.io import sequence_hash
from crispr_phage_predictor.ml.dataset import load_repeat_cas_training_table


AUDIT_COLUMNS = [
    "metric",
    "value",
]


def audit_crisprcasdb_candidates(
    current_training_csv: str | Path,
    candidate_csv: str | Path,
    max_per_subtype: int = 500,
) -> dict[str, pd.DataFrame]:
    """Compare current repeat/Cas training data with CRISPRCasdb candidates."""
    current = load_repeat_cas_training_table(current_training_csv).copy()
    candidates = load_repeat_cas_training_table(candidate_csv).copy()
    current["_repeat_hash"] = current["repeat_sequence"].map(sequence_hash)
    candidates["_repeat_hash"] = candidates["repeat_sequence"].map(sequence_hash)

    current_hashes = set(current["_repeat_hash"])
    candidate_hashes = set(candidates["_repeat_hash"])
    conflict_hashes = _conflicting_repeat_hashes(candidates)
    nonconflicting = candidates[~candidates["_repeat_hash"].isin(conflict_hashes)].copy()
    novel_nonconflicting = nonconflicting[~nonconflicting["_repeat_hash"].isin(current_hashes)].copy()
    balanced = _balanced_candidate_subset(novel_nonconflicting, max_per_subtype=max_per_subtype)

    summary = pd.DataFrame(
        [
            ("current_rows", len(current)),
            ("candidate_rows", len(candidates)),
            ("current_unique_repeat_hashes", len(current_hashes)),
            ("candidate_unique_repeat_hashes", len(candidate_hashes)),
            ("overlapping_repeat_hashes", len(current_hashes & candidate_hashes)),
            ("candidate_conflicting_repeat_hashes", len(conflict_hashes)),
            ("candidate_rows_after_conflict_filter", len(nonconflicting)),
            ("novel_nonconflicting_candidate_rows", len(novel_nonconflicting)),
            ("balanced_candidate_rows", len(balanced)),
            ("max_per_subtype", max_per_subtype),
        ],
        columns=AUDIT_COLUMNS,
    )

    return {
        "summary": summary,
        "current_subtype_counts": _count_by(current, "cas_subtype"),
        "candidate_subtype_counts": _count_by(candidates, "cas_subtype"),
        "novel_nonconflicting_subtype_counts": _count_by(novel_nonconflicting, "cas_subtype"),
        "balanced_candidate_subtype_counts": _count_by(balanced, "cas_subtype"),
        "candidate_repeat_conflicts": _candidate_repeat_conflicts(candidates),
        "subtype_delta": _subtype_delta(current, candidates, novel_nonconflicting, balanced),
    }


def _conflicting_repeat_hashes(table: pd.DataFrame) -> set[str]:
    subtype_counts = table.groupby("_repeat_hash")["cas_subtype"].nunique()
    return set(subtype_counts[subtype_counts > 1].index)


def _candidate_repeat_conflicts(table: pd.DataFrame) -> pd.DataFrame:
    conflict_hashes = _conflicting_repeat_hashes(table)
    if not conflict_hashes:
        return pd.DataFrame(columns=["repeat_hash", "subtypes", "row_count", "example_repeat_sequence"])
    conflicts = table[table["_repeat_hash"].isin(conflict_hashes)].copy()
    grouped = (
        conflicts.groupby("_repeat_hash")
        .agg(
            subtypes=("cas_subtype", lambda values: ";".join(sorted(set(map(str, values))))),
            row_count=("cas_subtype", "size"),
            example_repeat_sequence=("repeat_sequence", "first"),
        )
        .reset_index()
        .rename(columns={"_repeat_hash": "repeat_hash"})
        .sort_values(["row_count", "repeat_hash"], ascending=[False, True])
    )
    return grouped


def _balanced_candidate_subset(table: pd.DataFrame, max_per_subtype: int) -> pd.DataFrame:
    if max_per_subtype <= 0 or table.empty:
        return table.copy()
    ranked = table.sort_values(
        ["cas_subtype", "spacer_count", "genome_id", "repeat_sequence"],
        ascending=[True, False, True, True],
    )
    return ranked.groupby("cas_subtype", group_keys=False).head(max_per_subtype).copy()


def _count_by(table: pd.DataFrame, column: str) -> pd.DataFrame:
    if table.empty:
        return pd.DataFrame(columns=[column, "count"])
    return (
        table.groupby(column)
        .size()
        .reset_index(name="count")
        .sort_values(["count", column], ascending=[False, True])
    )


def _subtype_delta(
    current: pd.DataFrame,
    candidates: pd.DataFrame,
    novel_nonconflicting: pd.DataFrame,
    balanced: pd.DataFrame,
) -> pd.DataFrame:
    current_counts = current.groupby("cas_subtype").size()
    candidate_counts = candidates.groupby("cas_subtype").size()
    novel_counts = novel_nonconflicting.groupby("cas_subtype").size()
    balanced_counts = balanced.groupby("cas_subtype").size()
    subtypes = sorted(
        set(current_counts.index)
        | set(candidate_counts.index)
        | set(novel_counts.index)
        | set(balanced_counts.index)
    )
    return pd.DataFrame(
        [
            {
                "cas_subtype": subtype,
                "current_count": int(current_counts.get(subtype, 0)),
                "candidate_count": int(candidate_counts.get(subtype, 0)),
                "novel_nonconflicting_count": int(novel_counts.get(subtype, 0)),
                "balanced_candidate_count": int(balanced_counts.get(subtype, 0)),
            }
            for subtype in subtypes
        ]
    ).sort_values(["balanced_candidate_count", "candidate_count", "cas_subtype"], ascending=[False, False, True])


def write_audit_outputs(outputs: dict[str, pd.DataFrame], output_dir: str | Path) -> None:
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    for name, table in outputs.items():
        table.to_csv(output_path / f"{name}.csv", index=False)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Audit CRISPRCasdb SQL candidate labels against the current SABR training table."
    )
    parser.add_argument("current_training_csv", help="Current SABR repeat/Cas training CSV.")
    parser.add_argument("candidate_csv", help="CRISPRCasdb SQL candidate repeat/Cas CSV.")
    parser.add_argument(
        "--output-dir",
        default="data/training/audits/crisprcasdb_sql_candidate",
        help="Directory for audit CSV outputs.",
    )
    parser.add_argument(
        "--max-per-subtype",
        type=int,
        default=500,
        help="Candidate cap per subtype for the proposed balanced subset.",
    )
    args = parser.parse_args()

    outputs = audit_crisprcasdb_candidates(
        args.current_training_csv,
        args.candidate_csv,
        max_per_subtype=args.max_per_subtype,
    )
    write_audit_outputs(outputs, args.output_dir)
    summary = outputs["summary"]
    print(summary.to_string(index=False))
    print(f"Wrote audit outputs to {args.output_dir}")


if __name__ == "__main__":
    main()
