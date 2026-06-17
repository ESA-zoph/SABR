from __future__ import annotations

from pathlib import Path

import pandas as pd


ACCESSION_LINKAGE_COLUMNS = [
    "linkage_id",
    "entity_type",
    "source_key",
    "display_name",
    "strain_or_isolate",
    "accession",
    "accession_database",
    "assembly_level",
    "sequence_status",
    "linkage_status",
    "confidence",
    "local_path",
    "notes",
]

REQUIRED_ACCESSION_LINKAGE_COLUMNS = [
    "linkage_id",
    "entity_type",
    "source_key",
    "display_name",
    "strain_or_isolate",
    "linkage_status",
    "confidence",
]

ALLOWED_ACCESSION_LINKAGE_VALUES = {
    "entity_type": {"bacterium", "phage", "cocktail"},
    "accession_database": {"GenBank", "RefSeq", "ENA", "DDBJ", "SRA", "unknown", ""},
    "assembly_level": {
        "complete_genome",
        "chromosome",
        "contig",
        "scaffold",
        "metagenome",
        "not_applicable",
        "unknown",
        "",
    },
    "sequence_status": {
        "available",
        "local",
        "not_found",
        "not_applicable",
        "needs_lookup",
        "unknown",
        "",
    },
    "linkage_status": {
        "exact",
        "strain_alias",
        "species_only",
        "source_panel_only",
        "reference_proxy",
        "not_applicable",
        "needs_lookup",
        "unresolved",
    },
    "confidence": {"high", "medium", "low", "exclude"},
}


def empty_accession_linkage_table() -> pd.DataFrame:
    return pd.DataFrame(columns=ACCESSION_LINKAGE_COLUMNS)


def load_accession_linkage_table(path: str | Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    validate_accession_linkage_table(table)
    return table


def validate_accession_linkage_table(table: pd.DataFrame) -> None:
    missing_columns = [
        column for column in REQUIRED_ACCESSION_LINKAGE_COLUMNS if column not in table.columns
    ]
    if missing_columns:
        raise ValueError(
            "Missing required accession-linkage columns: " + ", ".join(missing_columns)
        )

    invalid_rows: list[str] = []
    seen_ids: set[str] = set()
    for index, row in table.iterrows():
        row_number = index + 2
        linkage_id = str(row.get("linkage_id", "")).strip()
        if not linkage_id:
            invalid_rows.append(f"row {row_number}: linkage_id is empty")
        elif linkage_id in seen_ids:
            invalid_rows.append(f"row {row_number}: duplicate linkage_id {linkage_id}")
        seen_ids.add(linkage_id)

        for column in REQUIRED_ACCESSION_LINKAGE_COLUMNS:
            if not str(row.get(column, "")).strip():
                invalid_rows.append(f"row {row_number}: {column} is empty")

        for column, allowed in ALLOWED_ACCESSION_LINKAGE_VALUES.items():
            if column not in table.columns:
                continue
            value = str(row.get(column, "")).strip()
            if value not in allowed:
                invalid_rows.append(f"row {row_number}: invalid {column} '{value}'")

        accession = str(row.get("accession", "")).strip()
        sequence_status = str(row.get("sequence_status", "")).strip()
        linkage_status = str(row.get("linkage_status", "")).strip()
        if linkage_status in {"exact", "strain_alias", "reference_proxy"} and not accession:
            invalid_rows.append(
                f"row {row_number}: exact/strain_alias/reference_proxy linkage requires accession"
            )
        if accession and sequence_status in {"not_found", "not_applicable"}:
            invalid_rows.append(
                f"row {row_number}: accession conflicts with sequence_status {sequence_status}"
            )

    if invalid_rows:
        raise ValueError("Invalid accession-linkage table:\n" + "\n".join(invalid_rows))


def accession_coverage(
    interactions: pd.DataFrame,
    linkage: pd.DataFrame,
) -> pd.DataFrame:
    validate_accession_linkage_table(linkage)
    rows = []
    for _, row in interactions.iterrows():
        bacterium_exact_ready = _has_linked_accession(
            linkage,
            entity_type="bacterium",
            source_key=row.get("source_key", ""),
            display_name=row.get("bacterium", ""),
            strain_or_isolate=row.get("strain", ""),
            allowed_statuses=("exact", "strain_alias"),
        ) or bool(str(row.get("bacterial_accession", "")).strip())
        bacterium_hybrid_ready = bacterium_exact_ready or _has_linked_accession(
            linkage,
            entity_type="bacterium",
            source_key=row.get("source_key", ""),
            display_name=row.get("bacterium", ""),
            strain_or_isolate=row.get("strain", ""),
            allowed_statuses=("reference_proxy",),
        )
        phage_ready = _has_linked_accession(
            linkage,
            entity_type="phage",
            source_key=row.get("source_key", ""),
            display_name=row.get("phage", ""),
            strain_or_isolate=row.get("phage", ""),
            allowed_statuses=("exact", "strain_alias"),
        ) or bool(str(row.get("phage_accession", "")).strip())
        pair_genome_ready = bacterium_exact_ready and phage_ready
        pair_hybrid_ready = bacterium_hybrid_ready and phage_ready
        rows.append(
            {
                "interaction_id": row.get("interaction_id", ""),
                "source_key": row.get("source_key", ""),
                "bacterium": row.get("bacterium", ""),
                "strain": row.get("strain", ""),
                "phage": row.get("phage", ""),
                "bacterium_genome_linked": bacterium_exact_ready,
                "bacterium_reference_or_genome_linked": bacterium_hybrid_ready,
                "phage_genome_linked": phage_ready,
                "pair_genome_ready": pair_genome_ready,
                "pair_hybrid_ready": pair_hybrid_ready,
                "dataset_tier": _dataset_tier(
                    pair_genome_ready=pair_genome_ready,
                    pair_hybrid_ready=pair_hybrid_ready,
                    bacterium_hybrid_ready=bacterium_hybrid_ready,
                    phage_ready=phage_ready,
                ),
            }
        )
    return pd.DataFrame(rows)


def _has_linked_accession(
    linkage: pd.DataFrame,
    entity_type: str,
    source_key: object,
    display_name: object,
    strain_or_isolate: object,
    allowed_statuses: tuple[str, ...],
) -> bool:
    exact = linkage[
        (linkage["entity_type"] == entity_type)
        & (linkage["source_key"] == str(source_key))
        & (linkage["display_name"] == str(display_name))
        & (linkage["strain_or_isolate"] == str(strain_or_isolate))
    ]
    if exact.empty:
        exact = linkage[
            (linkage["entity_type"] == entity_type)
            & (linkage["source_key"] == str(source_key))
            & (linkage["display_name"] == str(display_name))
            & (linkage["strain_or_isolate"] == "*")
        ]
    if exact.empty:
        return False
    linked = exact[
        exact["linkage_status"].isin(allowed_statuses)
        & exact["accession"].astype(str).str.strip().astype(bool)
    ]
    return not linked.empty


def _dataset_tier(
    pair_genome_ready: bool,
    pair_hybrid_ready: bool,
    bacterium_hybrid_ready: bool,
    phage_ready: bool,
) -> str:
    if pair_genome_ready:
        return "tier1_exact_pair"
    if pair_hybrid_ready:
        return "tier2_proxy_host_exact_phage"
    if phage_ready:
        return "tier3_phage_only"
    if bacterium_hybrid_ready:
        return "tier4_host_only_or_proxy"
    return "tier5_phenotype_only"
