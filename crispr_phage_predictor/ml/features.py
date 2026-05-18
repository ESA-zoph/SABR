from __future__ import annotations

from collections import Counter
from itertools import product

import pandas as pd


DEFAULT_KMER_SIZES = (2, 3, 4)
DNA_ALPHABET = "ACGT"


def kmer_counts(sequence: str, k: int = 3) -> dict[str, int]:
    normalized = sequence.upper()
    if k <= 0:
        raise ValueError("k must be greater than zero")
    if len(normalized) < k:
        return {}
    return Counter(normalized[index : index + k] for index in range(len(normalized) - k + 1))


def build_repeat_feature_table(
    training_table: pd.DataFrame,
    kmer_sizes: tuple[int, ...] = DEFAULT_KMER_SIZES,
    normalize_kmers: bool = True,
) -> pd.DataFrame:
    """Convert repeat/Cas training rows into numeric model features."""
    rows = []
    for _, row in training_table.iterrows():
        repeat = str(row["repeat_sequence"]).upper()
        feature_row = {
            "genome_id": row.get("genome_id", ""),
            "contig_id": row.get("contig_id", ""),
            "cas_type": row.get("cas_type", ""),
            "cas_subtype": row.get("cas_subtype", ""),
            "repeat_length": int(row.get("repeat_length", len(repeat))),
            "repeat_gc_percent": _gc_percent(repeat),
            "spacer_count": int(row.get("spacer_count", 0)),
            "mean_spacer_length": float(row.get("mean_spacer_length", 0.0) or 0.0),
        }
        for k in kmer_sizes:
            feature_row.update(_kmer_feature_values(repeat, k, normalize=normalize_kmers))
        rows.append(feature_row)
    return pd.DataFrame(rows)


def feature_columns(feature_table: pd.DataFrame) -> list[str]:
    excluded = {"genome_id", "contig_id", "cas_type", "cas_subtype"}
    return [column for column in feature_table.columns if column not in excluded]


def _gc_percent(sequence: str) -> float:
    valid_bases = [base for base in sequence.upper() if base in DNA_ALPHABET]
    if not valid_bases:
        return 0.0
    gc_count = valid_bases.count("G") + valid_bases.count("C")
    return round((gc_count / len(valid_bases)) * 100, 6)


def _kmer_feature_values(sequence: str, k: int, normalize: bool) -> dict[str, float]:
    observed = kmer_counts(_strip_ambiguous_bases(sequence), k=k)
    total = sum(observed.values())
    values: dict[str, float] = {}
    for kmer in _all_dna_kmers(k):
        key = f"kmer_{k}_{kmer}"
        count = observed.get(kmer, 0)
        values[key] = (count / total) if normalize and total else float(count)
    return values


def _all_dna_kmers(k: int) -> list[str]:
    if k <= 0:
        raise ValueError("k must be greater than zero")
    return ["".join(parts) for parts in product(DNA_ALPHABET, repeat=k)]


def _strip_ambiguous_bases(sequence: str) -> str:
    return "".join(base for base in sequence.upper() if base in DNA_ALPHABET)
