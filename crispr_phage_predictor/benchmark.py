from __future__ import annotations

from pathlib import Path

import pandas as pd


BENCHMARK_LABEL_COLUMNS = [
    "pair_id",
    "label_version",
    "label_status",
    "benchmark_split",
    "bacterium",
    "strain",
    "bacterial_accession",
    "phage",
    "phage_accession",
    "local_bacterium_file",
    "local_phage_file",
    "phenotype_label",
    "crispr_resistance_label",
    "crispr_evidence_level",
    "pam_evidence_level",
    "anti_crispr_status",
    "host_range_status",
    "expected_sabr_behavior",
    "source_keys",
    "curation_confidence",
    "notes",
]

REQUIRED_BENCHMARK_COLUMNS = [
    "pair_id",
    "label_status",
    "bacterium",
    "strain",
    "phage",
    "phenotype_label",
    "crispr_resistance_label",
    "source_keys",
    "curation_confidence",
]

ALLOWED_VALUES = {
    "label_status": {
        "validated",
        "candidate",
        "challenge",
        "exclude_until_verified",
        "needs_literature_review",
    },
    "benchmark_split": {"train_calibration", "validation", "test", "holdout", ""},
    "phenotype_label": {
        "resistant",
        "susceptible",
        "partially_resistant",
        "nonhost",
        "unknown",
    },
    "crispr_resistance_label": {
        "crispr_resistant",
        "not_crispr_resistant",
        "spacer_only",
        "unknown",
        "exclude",
    },
    "crispr_evidence_level": {
        "experimental",
        "spacer_pam_seed",
        "spacer_pam",
        "spacer_only",
        "none",
        "unknown",
    },
    "pam_evidence_level": {
        "validated_pam",
        "candidate_pam",
        "pam_absent",
        "not_evaluated",
        "not_applicable",
        "unknown",
    },
    "anti_crispr_status": {"present", "absent", "unknown", "not_applicable"},
    "host_range_status": {"host", "nonhost", "unknown"},
    "expected_sabr_behavior": {
        "high_score_expected",
        "moderate_score_expected",
        "low_score_expected",
        "no_score_expected",
        "do_not_score",
        "unknown",
    },
    "curation_confidence": {"high", "medium", "low", "exclude"},
}


def empty_benchmark_label_table() -> pd.DataFrame:
    return pd.DataFrame(columns=BENCHMARK_LABEL_COLUMNS)


def load_benchmark_label_table(path: str | Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    validate_benchmark_label_table(table)
    return table


def validate_benchmark_label_table(table: pd.DataFrame) -> None:
    missing_columns = [column for column in REQUIRED_BENCHMARK_COLUMNS if column not in table.columns]
    if missing_columns:
        raise ValueError("Missing required benchmark columns: " + ", ".join(missing_columns))

    invalid_rows: list[str] = []
    seen_pair_ids: set[str] = set()
    for index, row in table.iterrows():
        row_number = index + 2
        pair_id = str(row.get("pair_id", "")).strip()
        if not pair_id:
            invalid_rows.append(f"row {row_number}: pair_id is empty")
        elif pair_id in seen_pair_ids:
            invalid_rows.append(f"row {row_number}: duplicate pair_id {pair_id}")
        seen_pair_ids.add(pair_id)

        for column in REQUIRED_BENCHMARK_COLUMNS:
            if not str(row.get(column, "")).strip():
                invalid_rows.append(f"row {row_number}: {column} is empty")

        for column, allowed in ALLOWED_VALUES.items():
            if column not in table.columns:
                continue
            value = str(row.get(column, "")).strip()
            if value not in allowed:
                invalid_rows.append(
                    f"row {row_number}: invalid {column} '{value}'"
                )

    if invalid_rows:
        raise ValueError("Invalid benchmark label table:\n" + "\n".join(invalid_rows))


def evaluate_benchmark_run(
    run_dir: str | Path,
    benchmark_path: str | Path,
) -> pd.DataFrame:
    run_path = Path(run_dir)
    benchmark = load_benchmark_label_table(benchmark_path)
    evidence = _load_run_table(run_path / "evidence_matrix.csv")
    bacteria = _load_optional_run_table(run_path / "bacterial_records.csv")
    phages = _load_optional_run_table(run_path / "phage_records.csv")

    rows = []
    for _, label in benchmark.iterrows():
        matched = _matching_evidence_rows(label, evidence)
        if matched.empty:
            rows.append(_evaluation_row(label, None, bacteria, phages))
            continue
        for _, evidence_row in matched.iterrows():
            rows.append(_evaluation_row(label, evidence_row, bacteria, phages))
    return pd.DataFrame(rows)


def summarize_benchmark_evaluation(evaluation: pd.DataFrame) -> pd.DataFrame:
    if evaluation.empty:
        return pd.DataFrame(
            columns=["label_status", "expected_sabr_behavior", "rows", "matched_rows"]
        )
    grouped = (
        evaluation.groupby(["label_status", "expected_sabr_behavior"], dropna=False)
        .agg(
            rows=("pair_id", "count"),
            matched_rows=("run_match_status", lambda values: int((values == "matched").sum())),
            mean_score=("crispr_targeting_score", _mean_numeric),
            mean_unique_spacers=("unique_matching_spacers", _mean_numeric),
        )
        .reset_index()
    )
    return grouped


def _load_run_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Required run output not found: {path}")
    return pd.read_csv(path, dtype=str).fillna("")


def _load_optional_run_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, dtype=str).fillna("")


def _matching_evidence_rows(label: pd.Series, evidence: pd.DataFrame) -> pd.DataFrame:
    if evidence.empty:
        return evidence
    bacterium_files = _split_file_list(str(label.get("local_bacterium_file", "")))
    phage_files = _split_file_list(str(label.get("local_phage_file", "")))
    if not bacterium_files or not phage_files:
        return evidence.iloc[0:0]
    return evidence[
        evidence["bacterium"].map(_basename).isin(bacterium_files)
        & evidence["phage"].map(_basename).isin(phage_files)
    ].copy()


def _evaluation_row(
    label: pd.Series,
    evidence_row: pd.Series | None,
    bacteria: pd.DataFrame,
    phages: pd.DataFrame,
) -> dict[str, object]:
    bacterium_files = _split_file_list(str(label.get("local_bacterium_file", "")))
    phage_files = _split_file_list(str(label.get("local_phage_file", "")))
    row = {
        "pair_id": label.get("pair_id", ""),
        "label_status": label.get("label_status", ""),
        "benchmark_split": label.get("benchmark_split", ""),
        "expected_sabr_behavior": label.get("expected_sabr_behavior", ""),
        "phenotype_label": label.get("phenotype_label", ""),
        "crispr_resistance_label": label.get("crispr_resistance_label", ""),
        "crispr_evidence_level": label.get("crispr_evidence_level", ""),
        "pam_evidence_level": label.get("pam_evidence_level", ""),
        "curation_confidence": label.get("curation_confidence", ""),
        "benchmark_bacterium_file": ";".join(bacterium_files),
        "benchmark_phage_file": ";".join(phage_files),
        "benchmark_bacterial_accession": label.get("bacterial_accession", ""),
        "benchmark_phage_accession": label.get("phage_accession", ""),
        "run_bacterial_accession": _lookup_accessions(bacteria, bacterium_files),
        "run_phage_accession": _lookup_accessions(phages, phage_files),
        "run_match_status": "matched" if evidence_row is not None else "missing_from_run",
    }
    for column in _EVIDENCE_COLUMNS:
        row[column] = _evidence_value(evidence_row, column)
    row["score_expectation_result"] = _score_expectation_result(row)
    return row


_EVIDENCE_COLUMNS = [
    "crispr_targeting_score",
    "experimental_pam_weighted_score",
    "hypothetical_resistance_score",
    "spacer_hits",
    "unique_matching_spacers",
    "best_identity_percent",
    "best_coverage_percent",
    "pam_support_level",
    "best_pam_compatibility_score",
    "mean_pam_compatibility_score",
    "seed_evaluated_hits",
    "best_seed_mismatches",
    "current_evidence_level",
    "interpretation",
]


def _evidence_value(evidence_row: pd.Series | None, column: str) -> str:
    if evidence_row is None:
        return ""
    if column == "crispr_targeting_score" and column not in evidence_row:
        return str(evidence_row.get("hypothetical_resistance_score", ""))
    if column == "hypothetical_resistance_score" and column not in evidence_row:
        return str(evidence_row.get("crispr_targeting_score", ""))
    return str(evidence_row.get(column, ""))


def _lookup_accessions(records: pd.DataFrame, file_names: list[str]) -> str:
    if records.empty or "source_file" not in records.columns or "accession" not in records.columns:
        return ""
    matches = records[records["source_file"].map(_basename).isin(file_names)]
    accessions = [
        str(accession)
        for accession in matches["accession"].tolist()
        if str(accession).strip()
    ]
    return ";".join(dict.fromkeys(accessions))


def _score_expectation_result(row: dict[str, object]) -> str:
    expected = str(row.get("expected_sabr_behavior", ""))
    if str(row.get("run_match_status")) != "matched":
        return "not_evaluated"
    score = _float_or_none(row.get("crispr_targeting_score"))
    if score is None:
        score = _float_or_none(row.get("hypothetical_resistance_score"))
    if score is None:
        unique_spacers = _float_or_none(row.get("unique_matching_spacers"))
        if unique_spacers is None:
            return "not_evaluated"
        score = _legacy_proxy_score(unique_spacers)
    thresholds = {
        "high_score_expected": (70.0, None),
        "moderate_score_expected": (35.0, 85.0),
        "low_score_expected": (None, 40.0),
        "no_score_expected": (None, 0.0),
    }
    if expected == "do_not_score":
        return "excluded"
    if expected not in thresholds:
        return "not_evaluated"
    minimum, maximum = thresholds[expected]
    if minimum is not None and score < minimum:
        return "below_expected"
    if maximum is not None and score > maximum:
        return "above_expected"
    return "meets_expectation"


def _legacy_proxy_score(unique_spacers: float) -> float:
    if unique_spacers >= 3:
        return 75.0
    if unique_spacers == 2:
        return 55.0
    if unique_spacers == 1:
        return 35.0
    return 0.0


def _split_file_list(value: str) -> list[str]:
    return [_basename(item) for item in value.split(";") if item.strip()]


def _basename(value: object) -> str:
    return Path(str(value).strip()).name


def _float_or_none(value: object) -> float | None:
    try:
        if value == "":
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def _mean_numeric(values: pd.Series) -> float:
    numeric = pd.to_numeric(values, errors="coerce").dropna()
    if numeric.empty:
        return 0.0
    return round(float(numeric.mean()), 3)
