from __future__ import annotations

from collections import Counter
from itertools import product
from math import sqrt
from statistics import median

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
            "repeat_at_percent": _at_percent(repeat),
            "repeat_gc_skew": _gc_skew(repeat),
            "spacer_count": int(row.get("spacer_count", 0)),
            "mean_spacer_length": float(row.get("mean_spacer_length", 0.0) or 0.0),
            "min_spacer_length": _float_or_default(row.get("min_spacer_length"), row.get("mean_spacer_length", 0.0)),
            "max_spacer_length": _float_or_default(row.get("max_spacer_length"), row.get("mean_spacer_length", 0.0)),
            "median_spacer_length": _float_or_default(row.get("median_spacer_length"), row.get("mean_spacer_length", 0.0)),
            "std_spacer_length": _float_or_default(row.get("std_spacer_length"), 0.0),
            "repeat_self_rc_identity": _self_reverse_complement_identity(repeat),
            "repeat_longest_inverted_stem": _longest_terminal_inverted_stem(repeat),
            "repeat_hairpin_score": _hairpin_score(repeat),
        }
        feature_row["repeat_count"] = feature_row["spacer_count"] + 1
        feature_row["array_length_estimate"] = _array_length_estimate(
            repeat_length=feature_row["repeat_length"],
            repeat_count=feature_row["repeat_count"],
            spacer_count=feature_row["spacer_count"],
            mean_spacer_length=feature_row["mean_spacer_length"],
        )
        feature_row["spacer_repeat_length_ratio"] = _safe_ratio(
            feature_row["mean_spacer_length"],
            feature_row["repeat_length"],
        )
        feature_row.update(_terminal_base_features(repeat, width=4))
        feature_row.update(_terminal_base_features(repeat, width=6))
        for k in kmer_sizes:
            feature_row.update(_kmer_feature_values(repeat, k, normalize=normalize_kmers))
            feature_row.update(
                _terminal_kmer_feature_values(repeat, k, width=8, normalize=normalize_kmers)
            )
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


def _at_percent(sequence: str) -> float:
    valid_bases = [base for base in sequence.upper() if base in DNA_ALPHABET]
    if not valid_bases:
        return 0.0
    at_count = valid_bases.count("A") + valid_bases.count("T")
    return round((at_count / len(valid_bases)) * 100, 6)


def _gc_skew(sequence: str) -> float:
    normalized = sequence.upper()
    g_count = normalized.count("G")
    c_count = normalized.count("C")
    return round(_safe_ratio(g_count - c_count, g_count + c_count), 6)


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


def _float_or_default(value: object, default: object) -> float:
    try:
        if pd.isna(value):
            raise ValueError
        return float(value)
    except (TypeError, ValueError):
        try:
            if pd.isna(default):
                return 0.0
            return float(default)
        except (TypeError, ValueError):
            return 0.0


def _safe_ratio(numerator: float, denominator: float) -> float:
    denominator = float(denominator)
    if denominator == 0.0:
        return 0.0
    return round(float(numerator) / denominator, 6)


def _array_length_estimate(
    repeat_length: int,
    repeat_count: int,
    spacer_count: int,
    mean_spacer_length: float,
) -> float:
    return round((repeat_length * repeat_count) + (spacer_count * mean_spacer_length), 6)


def _terminal_base_features(sequence: str, width: int) -> dict[str, float]:
    normalized = _strip_ambiguous_bases(sequence)
    start = normalized[:width]
    end = normalized[-width:] if normalized else ""
    features: dict[str, float] = {}
    for prefix, segment in [(f"repeat_start_{width}", start), (f"repeat_end_{width}", end)]:
        features[f"{prefix}_gc_percent"] = _gc_percent(segment)
        for base in DNA_ALPHABET:
            features[f"{prefix}_{base}_fraction"] = _safe_ratio(segment.count(base), len(segment))
    return features


def _terminal_kmer_feature_values(
    sequence: str,
    k: int,
    width: int,
    normalize: bool,
) -> dict[str, float]:
    normalized = _strip_ambiguous_bases(sequence)
    terminal = normalized[:width] + normalized[-width:]
    observed = kmer_counts(terminal, k=k)
    total = sum(observed.values())
    values: dict[str, float] = {}
    for kmer in _all_dna_kmers(k):
        key = f"terminal_kmer_{k}_{kmer}"
        count = observed.get(kmer, 0)
        values[key] = (count / total) if normalize and total else float(count)
    return values


def _reverse_complement(sequence: str) -> str:
    complement = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return sequence.translate(complement)[::-1].upper()


def _self_reverse_complement_identity(sequence: str) -> float:
    normalized = _strip_ambiguous_bases(sequence)
    if not normalized:
        return 0.0
    reverse_complement = _reverse_complement(normalized)
    matches = sum(
        1 for left_base, right_base in zip(normalized, reverse_complement) if left_base == right_base
    )
    return round(matches / len(normalized), 6)


def _longest_terminal_inverted_stem(sequence: str, max_stem: int = 12) -> int:
    normalized = _strip_ambiguous_bases(sequence)
    longest = 0
    for stem_length in range(2, min(max_stem, len(normalized) // 2) + 1):
        left = normalized[:stem_length]
        right = normalized[-stem_length:]
        if left == _reverse_complement(right):
            longest = stem_length
    return longest


def _hairpin_score(sequence: str) -> float:
    normalized = _strip_ambiguous_bases(sequence)
    if len(normalized) < 8:
        return 0.0
    stem = _longest_terminal_inverted_stem(normalized)
    self_identity = _self_reverse_complement_identity(normalized)
    length_scale = sqrt(len(normalized))
    return round((stem / length_scale) + self_identity, 6)
