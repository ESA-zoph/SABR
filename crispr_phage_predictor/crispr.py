from __future__ import annotations

from dataclasses import dataclass
from statistics import mean
from typing import Callable


DEFAULT_MIN_REPEAT_LENGTH = 23
DEFAULT_MAX_REPEAT_LENGTH = 47
DEFAULT_MIN_SPACER_LENGTH = 20
DEFAULT_MAX_SPACER_LENGTH = 72
DEFAULT_MIN_REPEATS = 3


@dataclass(frozen=True)
class CrisprArray:
    array_id: str
    genome_id: str
    contig_id: str
    start: int
    end: int
    repeat_consensus: str
    spacers: list[str]

    @property
    def repeat_length(self) -> int:
        return len(self.repeat_consensus)

    @property
    def repeat_count(self) -> int:
        return len(self.spacers) + 1

    @property
    def spacer_count(self) -> int:
        return len(self.spacers)

    @property
    def mean_spacer_length(self) -> float:
        if not self.spacers:
            return 0.0
        return mean(len(spacer) for spacer in self.spacers)


def detect_crispr_arrays(
    sequence: str,
    genome_id: str,
    contig_id: str,
    min_repeat_length: int = DEFAULT_MIN_REPEAT_LENGTH,
    max_repeat_length: int = DEFAULT_MAX_REPEAT_LENGTH,
    min_spacer_length: int = DEFAULT_MIN_SPACER_LENGTH,
    max_spacer_length: int = DEFAULT_MAX_SPACER_LENGTH,
    min_repeats: int = DEFAULT_MIN_REPEATS,
    scan_progress_callback: Callable[[int, int], None] | None = None,
) -> list[CrisprArray]:
    """Detect candidate CRISPR arrays using exact direct repeats.

    This is an MVP detector intended to be transparent and easy to benchmark.
    It finds repeated sequences of CRISPR-like length separated by plausible
    spacer lengths. Later versions should add approximate repeats, orientation
    prediction, and comparison against established CRISPR callers.
    """
    normalized = _normalize_dna(sequence)
    candidates: list[CrisprArray] = []

    repeat_lengths = list(range(min_repeat_length, max_repeat_length + 1))
    for repeat_index, repeat_length in enumerate(repeat_lengths, start=1):
        repeat_positions = _index_repeat_positions(normalized, repeat_length)
        for repeat, positions in repeat_positions.items():
            if len(positions) < min_repeats or not _is_candidate_repeat(repeat):
                continue
            candidates.extend(
                _chains_to_arrays(
                    positions=positions,
                    repeat=repeat,
                    sequence=normalized,
                    genome_id=genome_id,
                    contig_id=contig_id,
                    min_spacer_length=min_spacer_length,
                    max_spacer_length=max_spacer_length,
                    min_repeats=min_repeats,
                )
            )
        if scan_progress_callback:
            scan_progress_callback(repeat_index, len(repeat_lengths))

    selected = _select_non_overlapping_arrays(candidates)
    return [
        CrisprArray(
            array_id=f"{genome_id}|{contig_id}|array_{index}",
            genome_id=array.genome_id,
            contig_id=array.contig_id,
            start=array.start,
            end=array.end,
            repeat_consensus=array.repeat_consensus,
            spacers=array.spacers,
        )
        for index, array in enumerate(selected, start=1)
    ]


def _normalize_dna(sequence: str) -> str:
    return "".join(base for base in sequence.upper() if not base.isspace())


def _index_repeat_positions(sequence: str, repeat_length: int) -> dict[str, list[int]]:
    positions_by_repeat: dict[str, list[int]] = {}
    for position in range(0, len(sequence) - repeat_length + 1):
        repeat = sequence[position : position + repeat_length]
        if set(repeat).issubset({"A", "C", "G", "T"}):
            positions_by_repeat.setdefault(repeat, []).append(position)
    return positions_by_repeat


def _is_candidate_repeat(repeat: str) -> bool:
    counts = {base: repeat.count(base) for base in "ACGT"}
    dominant_fraction = max(counts.values()) / len(repeat)
    gc_fraction = (counts["G"] + counts["C"]) / len(repeat)
    return dominant_fraction <= 0.8 and 0.2 <= gc_fraction <= 0.8


def _chains_to_arrays(
    positions: list[int],
    repeat: str,
    sequence: str,
    genome_id: str,
    contig_id: str,
    min_spacer_length: int,
    max_spacer_length: int,
    min_repeats: int,
) -> list[CrisprArray]:
    arrays: list[CrisprArray] = []
    repeat_length = len(repeat)
    chain = [positions[0]]

    for position in positions[1:]:
        spacer_length = position - chain[-1] - repeat_length
        if min_spacer_length <= spacer_length <= max_spacer_length:
            chain.append(position)
        else:
            arrays.extend(
                _chain_to_array(
                    chain=chain,
                    repeat=repeat,
                    sequence=sequence,
                    genome_id=genome_id,
                    contig_id=contig_id,
                    min_repeats=min_repeats,
                )
            )
            chain = [position]

    arrays.extend(
        _chain_to_array(
            chain=chain,
            repeat=repeat,
            sequence=sequence,
            genome_id=genome_id,
            contig_id=contig_id,
            min_repeats=min_repeats,
        )
    )
    return arrays


def _chain_to_array(
    chain: list[int],
    repeat: str,
    sequence: str,
    genome_id: str,
    contig_id: str,
    min_repeats: int,
) -> list[CrisprArray]:
    if len(chain) < min_repeats:
        return []

    repeat_length = len(repeat)
    spacers = [
        sequence[left + repeat_length : right]
        for left, right in zip(chain, chain[1:])
    ]
    return [
        CrisprArray(
            array_id="pending",
            genome_id=genome_id,
            contig_id=contig_id,
            start=chain[0] + 1,
            end=chain[-1] + repeat_length,
            repeat_consensus=repeat,
            spacers=spacers,
        )
    ]


def _select_non_overlapping_arrays(candidates: list[CrisprArray]) -> list[CrisprArray]:
    ranked = sorted(
        candidates,
        key=lambda array: (
            array.repeat_count,
            array.spacer_count,
            array.end - array.start,
            array.repeat_length,
        ),
        reverse=True,
    )
    selected: list[CrisprArray] = []
    for candidate in ranked:
        if any(_overlap_fraction(candidate, existing) > 0.5 for existing in selected):
            continue
        selected.append(candidate)
    return sorted(selected, key=lambda array: (array.contig_id, array.start, array.end))


def _overlap_fraction(left: CrisprArray, right: CrisprArray) -> float:
    if left.contig_id != right.contig_id or left.genome_id != right.genome_id:
        return 0.0
    overlap_start = max(left.start, right.start)
    overlap_end = min(left.end, right.end)
    overlap = max(0, overlap_end - overlap_start + 1)
    shorter = min(left.end - left.start + 1, right.end - right.start + 1)
    if shorter <= 0:
        return 0.0
    return overlap / shorter
