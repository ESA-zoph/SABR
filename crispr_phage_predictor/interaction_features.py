from __future__ import annotations

from collections import Counter
from itertools import product
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from crispr_phage_predictor.accession_linkage import load_accession_linkage_table


DNA_ALPHABET = "ACGT"
DEFAULT_KMER_SIZE = 3


def build_hybrid_interaction_feature_table(
    interactions: pd.DataFrame,
    linkage: pd.DataFrame,
    coverage: pd.DataFrame,
    k: int = DEFAULT_KMER_SIZE,
) -> pd.DataFrame:
    rows = []
    fasta_cache: dict[str, dict[str, object]] = {}
    hybrid_ids = set(
        coverage.loc[coverage["pair_hybrid_ready"].astype(bool), "interaction_id"].astype(str)
    )
    tiers_by_id = dict(zip(coverage["interaction_id"].astype(str), coverage.get("dataset_tier", "")))
    for _, interaction in interactions.iterrows():
        interaction_id = str(interaction["interaction_id"])
        if interaction_id not in hybrid_ids:
            continue
        host_link = _select_linkage(
            linkage,
            entity_type="bacterium",
            source_key=interaction["source_key"],
            display_name=interaction["bacterium"],
            strain_or_isolate=interaction["strain"],
            allow_proxy=True,
        )
        phage_link = _select_linkage(
            linkage,
            entity_type="phage",
            source_key=interaction["source_key"],
            display_name=interaction["phage"],
            strain_or_isolate=interaction["phage"],
            allow_proxy=False,
        )
        if host_link is None or phage_link is None:
            continue
        host_stats = _fasta_stats_cached(str(host_link["local_path"]), fasta_cache, prefix="host", k=k)
        phage_stats = _fasta_stats_cached(str(phage_link["local_path"]), fasta_cache, prefix="phage", k=k)
        if not host_stats or not phage_stats:
            continue
        eop_value = _float_or_none(interaction.get("eop_value", ""))
        row = {
            "interaction_id": interaction_id,
            "source_key": interaction["source_key"],
            "bacterium": interaction["bacterium"],
            "strain": interaction["strain"],
            "phage": interaction["phage"],
            "eop_class": interaction["eop_class"],
            "susceptibility_label": interaction["susceptibility_label"],
            "binary_susceptibility": _binary_susceptibility(interaction["susceptibility_label"]),
            "eop_value": "" if eop_value is None else eop_value,
            "host_linkage_status": host_link["linkage_status"],
            "host_accession": host_link["accession"],
            "host_local_path": host_link["local_path"],
            "phage_linkage_status": phage_link["linkage_status"],
            "phage_accession": phage_link["accession"],
            "phage_local_path": phage_link["local_path"],
            "uses_reference_proxy_host": host_link["linkage_status"] == "reference_proxy",
            "dataset_tier": tiers_by_id.get(interaction_id, ""),
        }
        row.update(host_stats)
        row.update(phage_stats)
        row["phage_to_host_length_ratio"] = _safe_ratio(
            row["phage_total_bp"], row["host_total_bp"]
        )
        row["phage_host_gc_delta"] = round(row["phage_gc_percent"] - row["host_gc_percent"], 6)
        rows.append(row)
    return pd.DataFrame(rows)


def load_inputs(
    interactions_path: str | Path,
    linkage_path: str | Path,
    coverage_path: str | Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    interactions = pd.read_csv(interactions_path, sep="\t", dtype=str).fillna("")
    linkage = load_accession_linkage_table(linkage_path)
    coverage = pd.read_csv(coverage_path, sep="\t", dtype=str).fillna("")
    for column in [
        "bacterium_genome_linked",
        "bacterium_reference_or_genome_linked",
        "phage_genome_linked",
        "pair_genome_ready",
        "pair_hybrid_ready",
    ]:
        if column in coverage.columns:
            coverage[column] = coverage[column].map(_truthy)
    return interactions, linkage, coverage


def _select_linkage(
    linkage: pd.DataFrame,
    entity_type: str,
    source_key: object,
    display_name: object,
    strain_or_isolate: object,
    allow_proxy: bool,
) -> pd.Series | None:
    candidates = linkage[
        (linkage["entity_type"] == entity_type)
        & (linkage["source_key"] == str(source_key))
        & (linkage["display_name"] == str(display_name))
        & (linkage["strain_or_isolate"].isin([str(strain_or_isolate), "*"]))
        & (linkage["local_path"].astype(str).str.strip() != "")
    ].copy()
    if candidates.empty:
        return None
    allowed = ["exact", "strain_alias"]
    if allow_proxy:
        allowed.append("reference_proxy")
    candidates = candidates[candidates["linkage_status"].isin(allowed)]
    if candidates.empty:
        return None
    priority = {"exact": 0, "strain_alias": 1, "reference_proxy": 2}
    candidates["_priority"] = candidates["linkage_status"].map(priority)
    return candidates.sort_values("_priority").iloc[0]


def _fasta_stats_cached(
    path: str,
    cache: dict[str, dict[str, object]],
    prefix: str,
    k: int,
) -> dict[str, object]:
    if not path:
        return {}
    if path not in cache:
        cache[path] = _fasta_stats(Path(path), k=k)
    return {f"{prefix}_{key}": value for key, value in cache[path].items()}


def _fasta_stats(path: Path, k: int) -> dict[str, object]:
    sequences = [str(record.seq).upper() for record in SeqIO.parse(path, "fasta")]
    sequence = "".join(sequences)
    valid_sequence = "".join(base for base in sequence if base in DNA_ALPHABET)
    length = len(valid_sequence)
    kmer_features = _kmer_frequencies(valid_sequence, k=k)
    return {
        "record_count": len(sequences),
        "total_bp": length,
        "gc_percent": _gc_percent(valid_sequence),
        "ambiguous_fraction": _safe_ratio(len(sequence) - length, len(sequence)),
        **kmer_features,
    }


def _kmer_frequencies(sequence: str, k: int) -> dict[str, float]:
    observed = Counter(sequence[index : index + k] for index in range(max(0, len(sequence) - k + 1)))
    total = sum(observed.values())
    features: dict[str, float] = {}
    for kmer in ("".join(parts) for parts in product(DNA_ALPHABET, repeat=k)):
        features[f"kmer_{k}_{kmer}_fraction"] = _safe_ratio(observed.get(kmer, 0), total)
    return features


def _gc_percent(sequence: str) -> float:
    if not sequence:
        return 0.0
    gc = sequence.count("G") + sequence.count("C")
    return round((gc / len(sequence)) * 100, 6)


def _safe_ratio(numerator: float, denominator: float) -> float:
    denominator = float(denominator)
    if denominator == 0:
        return 0.0
    return round(float(numerator) / denominator, 8)


def _float_or_none(value: object) -> float | None:
    try:
        if value == "":
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def _binary_susceptibility(label: object) -> str:
    normalized = str(label).strip()
    if normalized in {"susceptible", "reduced_susceptibility"}:
        return "susceptible"
    if normalized in {"resistant", "nonhost"}:
        return "resistant"
    return "unknown"


def _truthy(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}
