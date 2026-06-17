from __future__ import annotations

import argparse
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from crispr_phage_predictor.accession_linkage import (
    load_accession_linkage_table,
    validate_accession_linkage_table,
)
from crispr_phage_predictor.interactions import (
    eop_class_from_value,
    load_interaction_table,
    susceptibility_from_eop_class,
)
from scripts.download_accession_linkage_genomes import _fasta_stats


REPORT_COLUMNS = ["severity", "check", "row_id", "message"]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Audit curated phage-host interaction database consistency."
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
        "--features",
        type=Path,
        default=Path("data/training/phage_host_interaction_features_with_annotations.tsv"),
    )
    parser.add_argument(
        "--downloaded-records",
        type=Path,
        default=Path("data/curation/downloaded_records.tsv"),
    )
    parser.add_argument(
        "--assembly-candidates",
        type=Path,
        default=Path("data/curation/host_assembly_candidates.tsv"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/curation/phage_host_database_audit.tsv"),
    )
    args = parser.parse_args()

    issues: list[dict[str, str]] = []
    interactions = _load_checked(args.interactions, load_interaction_table, issues, "interactions")
    linkage = _load_checked(args.linkage, load_accession_linkage_table, issues, "accession_linkage")
    coverage = _load_table(args.coverage, issues, "coverage")
    features = _load_table(args.features, issues, "features")
    downloaded = _load_table(args.downloaded_records, issues, "downloaded_records")
    assemblies = _load_table(args.assembly_candidates, issues, "host_assembly_candidates")

    if interactions is not None:
        _audit_interactions(interactions, issues)
    if linkage is not None:
        _audit_linkage(linkage, issues)
    if downloaded is not None:
        _audit_downloaded_records(downloaded, issues)
    if interactions is not None and coverage is not None:
        _audit_coverage(interactions, coverage, issues)
    if coverage is not None and features is not None:
        _audit_features(coverage, features, issues)
    if linkage is not None and assemblies is not None:
        _audit_promoted_assemblies(linkage, assemblies, issues)

    report = pd.DataFrame(issues, columns=REPORT_COLUMNS)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(args.output, sep="\t", index=False)

    if report.empty:
        print("audit_status\tclean")
        print("issues\t0")
        return

    counts = report.groupby("severity").size().to_dict()
    print("audit_status\tissues_found")
    for severity in ["error", "warning"]:
        print(f"{severity}s\t{counts.get(severity, 0)}")
    print(f"report\t{args.output.as_posix()}")
    if counts.get("error", 0):
        raise SystemExit(1)


def _load_checked(path: Path, loader, issues: list[dict[str, str]], check: str) -> pd.DataFrame | None:
    try:
        return loader(path)
    except Exception as exc:
        _issue(issues, "error", check, path.as_posix(), str(exc))
        return None


def _load_table(path: Path, issues: list[dict[str, str]], check: str) -> pd.DataFrame | None:
    if not path.exists():
        _issue(issues, "error", check, path.as_posix(), "File does not exist.")
        return None
    try:
        return pd.read_csv(path, sep="\t", dtype=str).fillna("")
    except Exception as exc:
        _issue(issues, "error", check, path.as_posix(), str(exc))
        return None


def _audit_interactions(table: pd.DataFrame, issues: list[dict[str, str]]) -> None:
    for _, row in table.iterrows():
        row_id = str(row["interaction_id"])
        eop_value = str(row.get("eop_value", "")).strip()
        eop_class = str(row.get("eop_class", "")).strip()
        label = str(row.get("susceptibility_label", "")).strip()
        if eop_value:
            value = float(eop_value)
            expected_class = eop_class_from_value(value, str(row.get("eop_relation", "=")).strip())
            if eop_class != expected_class:
                _issue(
                    issues,
                    "error",
                    "eop_class_consistency",
                    row_id,
                    f"eop_value={eop_value} implies {expected_class}, found {eop_class}.",
                )
        expected_label = susceptibility_from_eop_class(eop_class)
        if expected_label != "unknown" and label not in {expected_label, "nonhost"}:
            _issue(
                issues,
                "warning",
                "susceptibility_label_consistency",
                row_id,
                f"eop_class={eop_class} usually maps to {expected_label}, found {label}.",
            )
        if str(row.get("curation_status", "")) == "curated" and not str(row.get("notes", "")).strip():
            _issue(issues, "warning", "curation_notes", row_id, "Curated row has no notes.")


def _audit_linkage(table: pd.DataFrame, issues: list[dict[str, str]]) -> None:
    validate_accession_linkage_table(table)
    linked = table[
        table["linkage_status"].isin(["exact", "strain_alias", "reference_proxy"])
        & table["sequence_status"].eq("available")
    ]
    for _, row in linked.iterrows():
        row_id = str(row["linkage_id"])
        local_path = Path(str(row.get("local_path", "")))
        if not str(row.get("local_path", "")).strip():
            _issue(issues, "error", "linkage_local_path", row_id, "Available linked row has no local_path.")
        elif not local_path.exists():
            _issue(
                issues,
                "error",
                "linkage_local_path",
                row_id,
                f"local_path does not exist: {local_path.as_posix()}",
            )
    exact_keys = table[
        table["linkage_status"].isin(["exact", "strain_alias"])
        & table["accession"].astype(str).str.strip().astype(bool)
    ][["entity_type", "source_key", "display_name", "strain_or_isolate"]]
    duplicated = exact_keys[exact_keys.duplicated(keep=False)]
    for index in duplicated.index:
        _issue(
            issues,
            "warning",
            "duplicate_exact_linkage_target",
            str(table.loc[index, "linkage_id"]),
            "Multiple exact/strain-alias links share the same entity/source/name/strain target.",
        )


def _audit_downloaded_records(table: pd.DataFrame, issues: list[dict[str, str]]) -> None:
    required = {"record_type", "name", "accession", "local_path", "record_count", "total_bp"}
    missing = required - set(table.columns)
    if missing:
        _issue(
            issues,
            "error",
            "downloaded_records_schema",
            "downloaded_records",
            "Missing columns: " + ", ".join(sorted(missing)),
        )
        return
    duplicate_accessions = table[table["accession"].astype(str).str.strip().duplicated(keep=False)]
    for _, row in duplicate_accessions.iterrows():
        _issue(
            issues,
            "warning",
            "duplicate_downloaded_accession",
            str(row["accession"]),
            "Downloaded records table contains duplicate accession entries.",
        )
    for _, row in table.iterrows():
        path = Path(str(row["local_path"]))
        if not path.exists():
            _issue(
                issues,
                "error",
                "downloaded_record_path",
                str(row["accession"]),
                f"Downloaded record path does not exist: {path.as_posix()}",
            )
            continue
        if path.suffix.lower() in {".fasta", ".fa", ".fna"}:
            stats = _fasta_stats(path)
            expected_count = _int_or_none(row.get("record_count", ""))
            expected_bp = _int_or_none(row.get("total_bp", ""))
            if expected_count is not None and stats["record_count"] != expected_count:
                _issue(
                    issues,
                    "error",
                    "downloaded_record_count",
                    str(row["accession"]),
                    f"record_count={expected_count} but FASTA has {stats['record_count']}.",
                )
            if expected_bp is not None and stats["total_bp"] != expected_bp:
                _issue(
                    issues,
                    "error",
                    "downloaded_record_bp",
                    str(row["accession"]),
                    f"total_bp={expected_bp} but FASTA has {stats['total_bp']}.",
                )


def _audit_coverage(
    interactions: pd.DataFrame,
    coverage: pd.DataFrame,
    issues: list[dict[str, str]],
) -> None:
    if len(interactions) != len(coverage):
        _issue(
            issues,
            "error",
            "coverage_row_count",
            "coverage",
            f"coverage rows={len(coverage)} but interaction rows={len(interactions)}.",
        )
    interaction_ids = set(interactions["interaction_id"].astype(str))
    coverage_ids = set(coverage["interaction_id"].astype(str))
    missing = sorted(interaction_ids - coverage_ids)
    extra = sorted(coverage_ids - interaction_ids)
    for row_id in missing[:20]:
        _issue(issues, "error", "coverage_missing_interaction", row_id, "Missing from coverage.")
    for row_id in extra[:20]:
        _issue(issues, "error", "coverage_extra_interaction", row_id, "Not present in interactions.")


def _audit_features(
    coverage: pd.DataFrame,
    features: pd.DataFrame,
    issues: list[dict[str, str]],
) -> None:
    if "interaction_id" not in features.columns:
        _issue(issues, "error", "features_schema", "features", "Missing interaction_id column.")
        return
    duplicated = features[features["interaction_id"].astype(str).duplicated(keep=False)]
    for _, row in duplicated.iterrows():
        _issue(
            issues,
            "error",
            "duplicate_feature_interaction",
            str(row["interaction_id"]),
            "Feature table contains duplicate interaction_id.",
        )
    expected_ids = set(
        coverage.loc[
            coverage["pair_hybrid_ready"].astype(str).str.lower().isin(["true", "1", "yes"]),
            "interaction_id",
        ].astype(str)
    )
    feature_ids = set(features["interaction_id"].astype(str))
    missing = sorted(expected_ids - feature_ids)
    extra = sorted(feature_ids - expected_ids)
    for row_id in missing[:20]:
        _issue(issues, "error", "features_missing_hybrid_ready_row", row_id, "Hybrid-ready row missing from features.")
    for row_id in extra[:20]:
        _issue(issues, "error", "features_extra_row", row_id, "Feature row is not hybrid-ready.")


def _audit_promoted_assemblies(
    linkage: pd.DataFrame,
    assemblies: pd.DataFrame,
    issues: list[dict[str, str]],
) -> None:
    promoted = assemblies[assemblies["review_status"].eq("promoted_exact_assembly")]
    promoted_accessions = set(promoted["assembly_accession"].astype(str))
    linked_promoted = linkage[
        linkage["linkage_id"].astype(str).str.startswith("host_assembly_")
        & linkage["linkage_status"].eq("exact")
    ]
    linked_accessions = set(linked_promoted["accession"].astype(str))
    for accession in sorted(promoted_accessions - linked_accessions):
        _issue(issues, "error", "promoted_assembly_missing_linkage", accession, "Promoted candidate not in linkage.")
    for accession in sorted(linked_accessions - promoted_accessions):
        _issue(issues, "error", "linked_assembly_missing_candidate", accession, "Linked host assembly lacks promoted candidate evidence.")
    for _, row in linked_promoted.iterrows():
        notes = str(row.get("notes", ""))
        if "full genome representation" not in notes:
            _issue(
                issues,
                "warning",
                "promoted_assembly_evidence",
                str(row["linkage_id"]),
                "Promoted assembly note does not mention full genome representation.",
            )


def _issue(
    issues: list[dict[str, str]],
    severity: str,
    check: str,
    row_id: str,
    message: str,
) -> None:
    issues.append(
        {
            "severity": severity,
            "check": check,
            "row_id": row_id,
            "message": message,
        }
    )


def _int_or_none(value: object) -> int | None:
    try:
        return int(str(value).strip())
    except (TypeError, ValueError):
        return None


if __name__ == "__main__":
    main()
