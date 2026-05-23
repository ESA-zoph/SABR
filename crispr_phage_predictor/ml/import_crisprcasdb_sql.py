from __future__ import annotations

import argparse
from collections.abc import Iterable
from pathlib import Path

import pandas as pd

from crispr_phage_predictor.ml.dataset import (
    REPEAT_CAS_DATASET_COLUMNS,
    cas_type_from_subtype,
    validate_repeat_cas_training_table,
)


VALID_BASES = {"A", "C", "G", "T", "N"}


def import_crisprcasdb_sql_candidate_labels(
    sql_path: str | Path,
    release: str = "34",
    max_cas_distance_bp: int = 20_000,
    min_evidence_level: int = 4,
) -> pd.DataFrame:
    """Build candidate repeat/Cas rows from a CRISPRCasdb PostgreSQL dump.

    Rows are computational candidates, not curated gold labels. A CRISPR locus
    is labeled only when its direct-repeat consensus can be resolved and the
    nearest same-sequence Cas cluster has an unambiguous CAS-Type subtype within
    the requested distance threshold.
    """
    parsed = _load_required_tables(sql_path)
    regions = parsed["region"]
    loci = parsed["crisprlocus"]
    locus_regions = parsed["crisprlocus_region"]
    clusters = parsed["clustercas"]
    sequences = parsed["sequence"]
    strains = parsed["strain"]

    region_by_id = {
        row["id"]: {"sequence": row["sequence"], "category": row["category"]}
        for row in regions
    }
    sequence_by_id = {row["id"]: row for row in sequences}
    strain_by_id = {row["id"]: row for row in strains}

    cluster_by_sequence: dict[str, list[dict[str, object]]] = {}
    for cluster in clusters:
        subtype = _normalize_crisprcasdb_subtype(cluster["class"])
        if not subtype:
            continue
        cluster = {**cluster, "cas_subtype": subtype}
        cluster_by_sequence.setdefault(str(cluster["sequence"]), []).append(cluster)

    spacer_stats = _build_spacer_stats(locus_regions, region_by_id)

    rows = []
    for locus in loci:
        if int(locus["evidencelevel"]) < min_evidence_level:
            continue
        repeat_region = region_by_id.get(str(locus["drconsensus"]))
        if not repeat_region:
            continue
        repeat = str(repeat_region["sequence"]).upper().strip()
        if not _is_usable_repeat(repeat):
            continue

        nearest = _nearest_cluster(
            locus_start=int(locus["start"]),
            locus_end=int(locus["start"]) + int(locus["length"]),
            clusters=cluster_by_sequence.get(str(locus["sequence"]), []),
        )
        if nearest is None or int(nearest["distance_bp"]) > max_cas_distance_bp:
            continue

        sequence = sequence_by_id.get(str(locus["sequence"]), {})
        strain = strain_by_id.get(str(sequence.get("strain", "")), {})
        refseq = _clean_null(str(strain.get("refseq", "")))
        genbank = _clean_null(str(strain.get("genbank", "")))
        genome_id = refseq or genbank or str(locus["sequence"])
        stats = spacer_stats.get(str(locus["id"]), {"spacer_count": 0, "mean_spacer_length": 0.0})
        subtype = str(nearest["cas_subtype"])
        rows.append(
            {
                "source": "crisprcasdb_sql_nearest_cas_cluster_candidate",
                "genome_id": genome_id,
                "organism": "",
                "taxonomy": "",
                "assembly_level": _clean_null(str(strain.get("assembly_status", ""))),
                "contig_id": str(locus["sequence"]),
                "array_start": int(locus["start"]),
                "array_end": int(locus["start"]) + int(locus["length"]),
                "repeat_sequence": repeat,
                "repeat_length": len(repeat),
                "spacer_count": int(stats["spacer_count"]),
                "mean_spacer_length": round(float(stats["mean_spacer_length"]), 6),
                "cas_type": cas_type_from_subtype(subtype),
                "cas_subtype": subtype,
                "label_source": (
                    f"CRISPRCasdb_release_{release}_nearest_same_sequence_cas_cluster"
                    f"_{int(nearest['distance_bp'])}bp"
                ),
                "label_confidence": "computational_nearby_cas_cluster",
                "pam_rule": "",
            }
        )

    table = pd.DataFrame(rows, columns=REPEAT_CAS_DATASET_COLUMNS)
    validate_repeat_cas_training_table(table)
    return table


def _load_required_tables(sql_path: str | Path) -> dict[str, list[dict[str, object]]]:
    columns = {
        "clustercas": ["id", "sequence", "start", "length", "class"],
        "crisprlocus": [
            "id",
            "sequence",
            "start",
            "length",
            "orientation",
            "trusted",
            "evidencelevel",
            "drconsensus",
            "drconservation",
            "spacerconservation",
            "potentialorientation",
            "evidencelevelreeval",
            "blastscore",
        ],
        "crisprlocus_region": ["crisprlocus", "region", "start", "length"],
        "region": ["id", "sequence", "category"],
        "sequence": ["id", "strain", "category", "length", "ncount", "description", "job"],
        "strain": [
            "id",
            "genbank",
            "refseq",
            "taxon",
            "gb_release_date",
            "release_level",
            "assembly_status",
        ],
    }
    parsed = {name: [] for name in columns}
    active_table = ""
    active_columns: list[str] = []
    with Path(sql_path).open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if active_table:
                if line == r"\.":
                    active_table = ""
                    active_columns = []
                    continue
                parsed[active_table].append(_parse_copy_row(line, active_columns))
                continue
            for table_name, table_columns in columns.items():
                if line.startswith(f"COPY public.{table_name} "):
                    active_table = table_name
                    active_columns = table_columns
                    break
    return parsed


def _parse_copy_row(line: str, columns: list[str]) -> dict[str, object]:
    values = [_parse_pg_value(value) for value in line.split("\t")]
    return dict(zip(columns, values, strict=True))


def _parse_pg_value(value: str) -> object:
    if value == r"\N":
        return ""
    return value.replace(r"\n", "").strip()


def _normalize_crisprcasdb_subtype(value: object) -> str:
    text = str(value).strip().upper()
    if not text.startswith("CAS-TYPE"):
        return ""
    suffix = text.replace("CAS-TYPE", "", 1)
    if not suffix or "-" not in suffix:
        return ""
    subtype = suffix
    if subtype.count("-") != 1:
        return ""
    cas_type, subtype_suffix = subtype.split("-", 1)
    if not cas_type or not subtype_suffix:
        return ""
    return f"{cas_type}-{subtype_suffix}"


def _build_spacer_stats(
    locus_regions: Iterable[dict[str, object]],
    region_by_id: dict[str, dict[str, object]],
) -> dict[str, dict[str, float]]:
    lengths_by_locus: dict[str, list[int]] = {}
    for link in locus_regions:
        region = region_by_id.get(str(link["region"]))
        if not region or int(region["category"]) != 3:
            continue
        lengths_by_locus.setdefault(str(link["crisprlocus"]), []).append(int(link["length"]))
    return {
        locus_id: {
            "spacer_count": len(lengths),
            "mean_spacer_length": sum(lengths) / len(lengths) if lengths else 0.0,
        }
        for locus_id, lengths in lengths_by_locus.items()
    }


def _nearest_cluster(
    locus_start: int,
    locus_end: int,
    clusters: list[dict[str, object]],
) -> dict[str, object] | None:
    nearest: dict[str, object] | None = None
    for cluster in clusters:
        cluster_start = int(cluster["start"])
        cluster_end = cluster_start + int(cluster["length"])
        distance = _interval_distance(locus_start, locus_end, cluster_start, cluster_end)
        candidate = {**cluster, "distance_bp": distance}
        if nearest is None or distance < int(nearest["distance_bp"]):
            nearest = candidate
    return nearest


def _interval_distance(left_start: int, left_end: int, right_start: int, right_end: int) -> int:
    if left_end < right_start:
        return right_start - left_end
    if right_end < left_start:
        return left_start - right_end
    return 0


def _is_usable_repeat(repeat: str) -> bool:
    return 23 <= len(repeat) <= 47 and set(repeat).issubset(VALID_BASES)


def _clean_null(value: str) -> str:
    return str(value).replace("\\n", "").strip()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Import candidate repeat/Cas labels from a CRISPRCasdb PostgreSQL dump."
    )
    parser.add_argument("sql_path", help="Extracted CRISPRCasdb PostgreSQL dump path.")
    parser.add_argument(
        "--output",
        default="data/training/repeats_cas_types_crisprcasdb_sql_candidate.csv",
        help="Output SABR repeat/Cas training CSV path.",
    )
    parser.add_argument("--release", default="34", help="CRISPRCasdb release/version label.")
    parser.add_argument(
        "--max-cas-distance-bp",
        type=int,
        default=20_000,
        help="Maximum distance between CRISPR locus and nearest same-sequence Cas cluster.",
    )
    parser.add_argument(
        "--min-evidence-level",
        type=int,
        default=4,
        help="Minimum CRISPRCasdb evidencelevel to keep.",
    )
    args = parser.parse_args()

    table = import_crisprcasdb_sql_candidate_labels(
        args.sql_path,
        release=args.release,
        max_cas_distance_bp=args.max_cas_distance_bp,
        min_evidence_level=args.min_evidence_level,
    )
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_path, index=False)
    print(f"Wrote {len(table)} candidate rows to {output_path}")


if __name__ == "__main__":
    main()
