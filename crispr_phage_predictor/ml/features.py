from __future__ import annotations

from collections import Counter


def kmer_counts(sequence: str, k: int = 3) -> dict[str, int]:
    normalized = sequence.upper()
    if k <= 0:
        raise ValueError("k must be greater than zero")
    if len(normalized) < k:
        return {}
    return Counter(normalized[index : index + k] for index in range(len(normalized) - k + 1))
