from __future__ import annotations

from pathlib import Path
from dataclasses import dataclass
from math import log2

import pandas as pd
from Bio import SeqIO

from crispr_phage_predictor.crispr import CrisprArray
from crispr_phage_predictor.crispr import detect_crispr_arrays
from crispr_phage_predictor.external.minced import detect_arrays_with_minced, minced_available
from crispr_phage_predictor.io import FastaRecord
from crispr_phage_predictor.matching import SpacerHit, find_spacer_hits, reverse_complement
from crispr_phage_predictor.scoring import score_crispr_targeting_evidence
from crispr_phage_predictor.pipeline import annotate_spacer_hits_with_pam


TARGETING_FEATURE_COLUMNS = [
    "crispr_array_count",
    "crispr_spacer_count",
    "crispr_mean_spacers_per_array",
    "spacer_hit_count",
    "unique_matching_spacers",
    "targeted_phage_contig_count",
    "best_spacer_identity",
    "best_spacer_coverage",
    "pam_evaluated_hit_count",
    "pam_supported_hit_count",
    "best_pam_compatibility_score",
    "seed_evaluated_hit_count",
    "best_seed_mismatches",
    "crispr_targeting_score",
    "fuzzy_spacer_candidate_count",
    "fuzzy_unique_spacers",
    "fuzzy_targeted_phage_contig_count",
    "best_fuzzy_spacer_identity",
    "best_fuzzy_spacer_mismatches",
    "best_fuzzy_seed_edge_mismatches",
    "best_fuzzy_distal_mismatches",
    "fuzzy_seed_perfect_hit_count",
    "fuzzy_near_perfect_hit_count",
    "fuzzy_high_confidence_hit_count",
    "graded_crispr_interference_score",
]


def add_targeting_features(
    feature_table: pd.DataFrame,
    internal_fallback_max_bp: int = 200_000,
) -> pd.DataFrame:
    host_array_cache: dict[str, list[CrisprArray]] = {}
    phage_record_cache: dict[str, list[FastaRecord]] = {}
    rows = []
    for _, row in feature_table.iterrows():
        host_path = str(row.get("host_local_path", ""))
        phage_path = str(row.get("phage_local_path", ""))
        arrays = _arrays_for_path(host_path, host_array_cache, internal_fallback_max_bp)
        phage_records = _records_for_path(phage_path, phage_record_cache)
        hits = find_spacer_hits(arrays, phage_records) if arrays and phage_records else []
        hits = annotate_spacer_hits_with_pam(hits)
        fuzzy_hits = (
            find_fuzzy_spacer_hits(arrays, phage_records) if arrays and phage_records else []
        )
        targeting = _targeting_features(arrays, hits, fuzzy_hits)
        merged = row.to_dict()
        merged.update(targeting)
        rows.append(merged)
    return pd.DataFrame(rows)


@dataclass(frozen=True)
class FuzzySpacerHit:
    spacer_id: str
    phage_contig_id: str
    strand: str
    start: int
    end: int
    identity: float
    mismatches: int
    spacer_length: int
    seed_edge_mismatches: int
    distal_mismatches: int


def find_fuzzy_spacer_hits(
    crispr_arrays: list[CrisprArray],
    phage_records: list[FastaRecord],
    max_mismatches: int = 4,
    min_identity: float = 0.85,
    seed_length: int = 8,
    max_hits_per_spacer: int = 5,
) -> list[FuzzySpacerHit]:
    """Find near spacer/protospacer matches with a bounded Hamming scan.

    This deliberately allows imperfect matches because CRISPR interference can
    tolerate distal mismatches while remaining sensitive to seed/PAM-proximal
    mismatches. Without a confident subtype/PAM orientation for every array, the
    seed-edge feature uses the better of the two spacer ends as a conservative
    "possible PAM-proximal seed" signal.
    """
    hits: list[FuzzySpacerHit] = []
    for array in crispr_arrays:
        for spacer_index, spacer in enumerate(array.spacers, start=1):
            spacer_id = f"{array.array_id}|spacer_{spacer_index}"
            spacer_hits: list[FuzzySpacerHit] = []
            query = spacer.upper()
            if len(query) < seed_length or set(query) - {"A", "C", "G", "T"}:
                continue
            for phage_record in phage_records:
                spacer_hits.extend(
                    _scan_fuzzy_one_spacer(
                        spacer_id=spacer_id,
                        spacer_sequence=query,
                        phage_contig_id=phage_record.record_id,
                        phage_sequence=phage_record.sequence.upper(),
                        max_mismatches=max_mismatches,
                        min_identity=min_identity,
                        seed_length=seed_length,
                    )
                )
            hits.extend(
                sorted(
                    spacer_hits,
                    key=lambda hit: (
                        hit.mismatches,
                        hit.seed_edge_mismatches,
                        -hit.identity,
                    ),
                )[:max_hits_per_spacer]
            )
    return hits


def _scan_fuzzy_one_spacer(
    spacer_id: str,
    spacer_sequence: str,
    phage_contig_id: str,
    phage_sequence: str,
    max_mismatches: int,
    min_identity: float,
    seed_length: int,
) -> list[FuzzySpacerHit]:
    hits = []
    for query, strand in [(spacer_sequence, "+"), (reverse_complement(spacer_sequence), "-")]:
        if strand == "-" and query == spacer_sequence:
            continue
        query_length = len(query)
        if len(phage_sequence) < query_length:
            continue
        for offset in range(0, len(phage_sequence) - query_length + 1):
            window = phage_sequence[offset : offset + query_length]
            mismatches = _bounded_hamming_mismatches(query, window, max_mismatches)
            if mismatches > max_mismatches:
                continue
            identity = 1.0 - (mismatches / query_length)
            if identity < min_identity:
                continue
            seed_edge_mismatches = min(
                _hamming_mismatches(query[:seed_length], window[:seed_length]),
                _hamming_mismatches(query[-seed_length:], window[-seed_length:]),
            )
            hits.append(
                FuzzySpacerHit(
                    spacer_id=spacer_id,
                    phage_contig_id=phage_contig_id,
                    strand=strand,
                    start=offset + 1,
                    end=offset + query_length,
                    identity=round(identity, 6),
                    mismatches=mismatches,
                    spacer_length=query_length,
                    seed_edge_mismatches=seed_edge_mismatches,
                    distal_mismatches=mismatches - seed_edge_mismatches,
                )
            )
    return hits


def _bounded_hamming_mismatches(left: str, right: str, max_mismatches: int) -> int:
    mismatches = 0
    for left_base, right_base in zip(left, right):
        if left_base != right_base:
            mismatches += 1
            if mismatches > max_mismatches:
                return mismatches
    return mismatches


def _hamming_mismatches(left: str, right: str) -> int:
    return sum(left_base != right_base for left_base, right_base in zip(left, right))


def _arrays_for_path(
    path: str,
    cache: dict[str, list[CrisprArray]],
    internal_fallback_max_bp: int,
) -> list[CrisprArray]:
    if not path:
        return []
    if path not in cache:
        records = _records_for_path(path, {})
        arrays = detect_arrays_with_minced(records) if minced_available() else []
        total_bp = sum(len(record.sequence) for record in records)
        if not arrays and total_bp <= internal_fallback_max_bp:
            for record in records:
                arrays.extend(
                    detect_crispr_arrays(
                        sequence=record.sequence,
                        genome_id=record.source_file,
                        contig_id=record.record_id,
                    )
                )
        cache[path] = arrays
    return cache[path]


def _records_for_path(path: str, cache: dict[str, list[FastaRecord]]) -> list[FastaRecord]:
    if not path:
        return []
    if path not in cache:
        fasta_path = Path(path)
        records = []
        for record in SeqIO.parse(fasta_path, "fasta"):
            records.append(
                FastaRecord(
                    source_file=fasta_path.name,
                    record_id=record.id,
                    description=record.description,
                    sequence=str(record.seq).upper(),
                )
            )
        cache[path] = records
    return cache[path]


def _targeting_features(
    arrays: list[CrisprArray],
    hits: list[SpacerHit],
    fuzzy_hits: list[FuzzySpacerHit] | None = None,
) -> dict[str, object]:
    fuzzy_hits = fuzzy_hits or []
    total_spacers = sum(array.spacer_count for array in arrays)
    score = score_crispr_targeting_evidence(hits)
    pam_scores = [
        hit.pam_compatibility_score
        for hit in hits
        if hit.pam_compatibility_score is not None
    ]
    seed_mismatches = [
        hit.seed_mismatches for hit in hits if hit.seed_mismatches is not None
    ]
    fuzzy_features = _fuzzy_targeting_features(fuzzy_hits)
    features = {
        "crispr_array_count": len(arrays),
        "crispr_spacer_count": total_spacers,
        "crispr_mean_spacers_per_array": (
            round(total_spacers / len(arrays), 6) if arrays else 0.0
        ),
        "spacer_hit_count": len(hits),
        "unique_matching_spacers": len({hit.spacer_id for hit in hits}),
        "targeted_phage_contig_count": len({hit.phage_contig_id for hit in hits}),
        "best_spacer_identity": max((hit.identity for hit in hits), default=0.0),
        "best_spacer_coverage": max((hit.coverage for hit in hits), default=0.0),
        "pam_evaluated_hit_count": sum(hit.pam_match is not None for hit in hits),
        "pam_supported_hit_count": sum(hit.pam_match is True for hit in hits),
        "best_pam_compatibility_score": max(pam_scores) if pam_scores else 0.0,
        "seed_evaluated_hit_count": len(seed_mismatches),
        "best_seed_mismatches": min(seed_mismatches) if seed_mismatches else -1,
        "crispr_targeting_score": score.score,
    }
    features.update(fuzzy_features)
    return features


def _fuzzy_targeting_features(fuzzy_hits: list[FuzzySpacerHit]) -> dict[str, object]:
    if not fuzzy_hits:
        return {
            "fuzzy_spacer_candidate_count": 0,
            "fuzzy_unique_spacers": 0,
            "fuzzy_targeted_phage_contig_count": 0,
            "best_fuzzy_spacer_identity": 0.0,
            "best_fuzzy_spacer_mismatches": -1,
            "best_fuzzy_seed_edge_mismatches": -1,
            "best_fuzzy_distal_mismatches": -1,
            "fuzzy_seed_perfect_hit_count": 0,
            "fuzzy_near_perfect_hit_count": 0,
            "fuzzy_high_confidence_hit_count": 0,
            "graded_crispr_interference_score": 0.0,
        }
    unique_spacers = {hit.spacer_id for hit in fuzzy_hits}
    best_hit = sorted(
        fuzzy_hits,
        key=lambda hit: (hit.seed_edge_mismatches, hit.mismatches, -hit.identity),
    )[0]
    high_confidence_hits = [
        hit
        for hit in fuzzy_hits
        if hit.identity >= 0.90 and hit.seed_edge_mismatches == 0 and hit.mismatches <= 3
    ]
    score = _graded_interference_score(fuzzy_hits, unique_spacers)
    return {
        "fuzzy_spacer_candidate_count": len(fuzzy_hits),
        "fuzzy_unique_spacers": len(unique_spacers),
        "fuzzy_targeted_phage_contig_count": len({hit.phage_contig_id for hit in fuzzy_hits}),
        "best_fuzzy_spacer_identity": max(hit.identity for hit in fuzzy_hits),
        "best_fuzzy_spacer_mismatches": min(hit.mismatches for hit in fuzzy_hits),
        "best_fuzzy_seed_edge_mismatches": best_hit.seed_edge_mismatches,
        "best_fuzzy_distal_mismatches": best_hit.distal_mismatches,
        "fuzzy_seed_perfect_hit_count": sum(hit.seed_edge_mismatches == 0 for hit in fuzzy_hits),
        "fuzzy_near_perfect_hit_count": sum(hit.mismatches <= 2 for hit in fuzzy_hits),
        "fuzzy_high_confidence_hit_count": len(high_confidence_hits),
        "graded_crispr_interference_score": score,
    }


def _graded_interference_score(
    fuzzy_hits: list[FuzzySpacerHit],
    unique_spacers: set[str],
) -> float:
    best_identity = max(hit.identity for hit in fuzzy_hits)
    best_seed = min(hit.seed_edge_mismatches for hit in fuzzy_hits)
    best_mismatches = min(hit.mismatches for hit in fuzzy_hits)
    seed_perfect = any(hit.seed_edge_mismatches == 0 for hit in fuzzy_hits)
    multi_spacer_bonus = min(20.0, log2(len(unique_spacers) + 1) * 10.0)
    identity_component = best_identity * 35.0
    seed_component = 25.0 if seed_perfect else max(0.0, 15.0 - (best_seed * 5.0))
    mismatch_component = max(0.0, 20.0 - (best_mismatches * 4.0))
    score = identity_component + seed_component + mismatch_component + multi_spacer_bonus
    return round(min(100.0, max(0.0, score)), 2)
