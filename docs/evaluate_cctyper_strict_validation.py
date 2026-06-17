from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from crispr_phage_predictor.ml.dataset import (
    REPEAT_CAS_DATASET_COLUMNS,
    cas_type_from_subtype,
    load_cctyper_crisprs_near_cas,
)
from crispr_phage_predictor.ml.evaluate_external_dataset import evaluate_external_dataset


PACKAGE_ROOT = Path("data") / "validation" / "cctyper_full_balanced_validation_package"
SUMMARY_CSV = PACKAGE_ROOT / "cctyper_genome_level_summary.csv"
RESULTS_ROOT = PACKAGE_ROOT / "cctyper_results_balanced"
TRAINING_TABLE = Path("data") / "training" / "repeats_cas_types_crisprcasdb_sql_candidate.csv"
OUTPUT_DIR = Path("data") / "validation" / "cctyper_strict_validation"


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    summary = pd.read_csv(SUMMARY_CSV)
    strict = summary[summary["StrictMatch"].astype(str).str.upper() == "YES"].copy()

    manifest = _build_manifest(strict)
    manifest_path = OUTPUT_DIR / "strict_manifest.csv"
    manifest.to_csv(manifest_path, index=False)

    full_table = _load_full_cctyper_table(manifest)
    full_table_path = OUTPUT_DIR / "strict_cctyper_all_near_cas_rows.csv"
    full_table.to_csv(full_table_path, index=False)

    expected_table = full_table[
        full_table["cas_subtype"].astype(str) == full_table["expected_subtype"].astype(str)
    ].copy()
    expected_table = expected_table.drop(columns=["expected_subtype"])
    expected_table_path = OUTPUT_DIR / "strict_cctyper_expected_subtype_rows.csv"
    expected_table.to_csv(expected_table_path, index=False)

    result = evaluate_external_dataset(TRAINING_TABLE, expected_table_path)
    _write_result_outputs(result)
    _write_summary(strict, full_table, expected_table, result)


def _build_manifest(strict: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in strict.iterrows():
        expected = str(row["ExpectedFolder"])
        genome = str(row["Genome"])
        output_dir = RESULTS_ROOT / expected / genome
        rows.append(
            {
                "expected_subtype": expected,
                "cctyper_output_dir": str(output_dir),
                "genome_id": genome,
                "organism": _organism_from_genome_name(genome, expected),
                "taxonomy": "",
                "assembly_level": "",
            }
        )
    return pd.DataFrame(rows)


def _load_full_cctyper_table(manifest: pd.DataFrame) -> pd.DataFrame:
    tables = []
    for _, row in manifest.iterrows():
        output_dir = Path(row["cctyper_output_dir"])
        near_cas = output_dir / "crisprs_near_cas.tab"
        if near_cas.exists():
            table = load_cctyper_crisprs_near_cas(
                near_cas,
                genome_id=str(row["genome_id"]),
                organism=str(row["organism"]),
                taxonomy=str(row["taxonomy"]),
                assembly_level=str(row["assembly_level"]),
            )
        else:
            table = _load_genome_level_crisprs_all(
                output_dir / "crisprs_all.tab",
                expected_subtype=str(row["expected_subtype"]),
                genome_id=str(row["genome_id"]),
                organism=str(row["organism"]),
                taxonomy=str(row["taxonomy"]),
                assembly_level=str(row["assembly_level"]),
            )
        table["expected_subtype"] = str(row["expected_subtype"])
        tables.append(table)
    return pd.concat(tables, ignore_index=True) if tables else pd.DataFrame()


def _load_genome_level_crisprs_all(
    crisprs_all: Path,
    expected_subtype: str,
    genome_id: str,
    organism: str,
    taxonomy: str,
    assembly_level: str,
) -> pd.DataFrame:
    cctyper_table = pd.read_csv(crisprs_all, sep="\t")
    rows = []
    for _, row in cctyper_table.iterrows():
        repeat = str(row["Consensus_repeat"]).upper().strip()
        if not repeat or not set(repeat).issubset({"A", "C", "G", "T", "N"}):
            continue
        rows.append(
            {
                "source": "cctyper_strict_genome_level",
                "genome_id": genome_id,
                "organism": organism,
                "taxonomy": taxonomy,
                "assembly_level": assembly_level,
                "contig_id": str(row["Contig"]),
                "array_start": _coerce_int(row.get("Start"), None),
                "array_end": _coerce_int(row.get("End"), None),
                "repeat_sequence": repeat,
                "repeat_length": len(repeat),
                "spacer_count": max(0, (_coerce_int(row.get("N_repeats"), 0) or 0) - 1),
                "mean_spacer_length": _coerce_float(row.get("Spacer_len_avg"), 0.0),
                "cas_type": cas_type_from_subtype(expected_subtype),
                "cas_subtype": expected_subtype,
                "label_source": "strict_cctyper_genome_operon",
                "label_confidence": "cas_operon_supported",
                "pam_rule": "",
            }
        )
    return pd.DataFrame(rows, columns=REPEAT_CAS_DATASET_COLUMNS)


def _write_result_outputs(result: dict[str, object]) -> None:
    (OUTPUT_DIR / "classification_report.txt").write_text(str(result["report"]), encoding="utf-8")
    result["confusion"].to_csv(OUTPUT_DIR / "confusion_matrix.csv")


def _write_summary(
    strict: pd.DataFrame,
    full_table: pd.DataFrame,
    expected_table: pd.DataFrame,
    result: dict[str, object],
) -> None:
    subtype_counts = (
        expected_table["cas_subtype"].value_counts().rename_axis("cas_subtype").reset_index(name="array_rows")
    )
    genome_counts = (
        strict["ExpectedFolder"].value_counts().rename_axis("cas_subtype").reset_index(name="strict_genomes")
    )
    counts = genome_counts.merge(subtype_counts, on="cas_subtype", how="outer").fillna(0)
    counts.to_csv(OUTPUT_DIR / "strict_validation_counts.csv", index=False)

    summary_rows = [
        {"metric": "strict_genomes", "value": len(strict)},
        {"metric": "full_imported_cctyper_rows", "value": len(full_table)},
        {"metric": "expected_subtype_rows", "value": len(expected_table)},
        {"metric": "train_rows", "value": result["train_rows"]},
        {"metric": "raw_test_rows", "value": result["raw_test_rows"]},
        {"metric": "evaluated_test_rows", "value": result["evaluated_test_rows"]},
        {"metric": "excluded_test_rows", "value": result["excluded_test_rows"]},
        {"metric": "accuracy", "value": result["accuracy"]},
        {
            "metric": "excluded_test_subtypes",
            "value": ";".join(result["excluded_test_subtypes"]),
        },
    ]
    pd.DataFrame(summary_rows).to_csv(OUTPUT_DIR / "summary.csv", index=False)


def _organism_from_genome_name(genome: str, expected: str) -> str:
    prefix = f"{expected}_"
    text = genome[len(prefix) :] if genome.startswith(prefix) else genome
    gcf_index = text.find("_GCF_")
    if gcf_index >= 0:
        text = text[:gcf_index]
    return text.replace("_", " ")


def _coerce_int(value: object, default: int | None) -> int | None:
    if pd.isna(value):
        return default
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def _coerce_float(value: object, default: float | None) -> float | None:
    if pd.isna(value):
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


if __name__ == "__main__":
    main()
