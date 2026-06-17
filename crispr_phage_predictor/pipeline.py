from __future__ import annotations

from dataclasses import replace
from dataclasses import dataclass
from typing import Literal
from typing import Callable

import pandas as pd

from crispr_phage_predictor.crispr import (
    DEFAULT_MAX_REPEAT_LENGTH,
    DEFAULT_MIN_REPEAT_LENGTH,
    CrisprArray,
    detect_crispr_arrays,
)
from crispr_phage_predictor.cas_prediction import ArrayCasPrediction
from crispr_phage_predictor.cas_prediction import CURATED_PAM_RULES_BY_SUBTYPE
from crispr_phage_predictor.external.blast import find_spacer_hits_with_blast
from crispr_phage_predictor.external.minced import detect_arrays_with_minced
from crispr_phage_predictor.io import FastaRecord
from crispr_phage_predictor.matching import find_spacer_hits
from crispr_phage_predictor.matching import reverse_complement
from crispr_phage_predictor.matching import summarize_seed_mismatches
from crispr_phage_predictor.matching import SpacerHit
from crispr_phage_predictor.pam import evaluate_pam_rule
from crispr_phage_predictor.scoring import score_crispr_targeting_evidence
from crispr_phage_predictor.scoring import score_experimental_pam_weighted_evidence

DetectionMethod = Literal["internal", "minced"]
MatchingMethod = Literal["internal", "blast"]

CRISPR_ARRAY_COLUMNS = [
    "array_id",
    "bacterium",
    "contig",
    "start",
    "end",
    "repeat_length",
    "repeat_count",
    "spacer_count",
    "mean_spacer_length",
    "repeat_consensus",
]

SPACER_COLUMNS = [
    "array_id",
    "spacer_id",
    "bacterium",
    "contig",
    "spacer_index",
    "spacer_length",
    "spacer_sequence",
]

SPACER_HIT_COLUMNS = [
    "bacterium",
    "phage",
    "array_id",
    "spacer_id",
    "phage_contig",
    "start",
    "end",
    "strand",
    "identity_percent",
    "mismatches",
    "alignment_length",
    "spacer_length",
    "coverage_percent",
    "evalue",
    "bitscore",
    "spacer_sequence",
    "aligned_spacer_sequence",
    "aligned_protospacer_sequence",
    "protospacer_sequence",
    "protospacer_5p_flank",
    "protospacer_3p_flank",
    "genomic_upstream_flank",
    "genomic_downstream_flank",
    "predicted_cas_subtype",
    "cas_subtype_confidence",
    "cas_subtype_prediction_method",
    "pam_rule",
    "pam_rule_source",
    "pam_sequence",
    "pam_match",
    "pam_support_level",
    "pam_compatibility_score",
    "pam_offset_from_protospacer",
    "seed_region",
    "seed_length",
    "seed_mismatches",
    "seed_mismatch_positions",
]

PAM_SUBTYPE_SUPPORT_COLUMNS = [
    "array_id",
    "bacterium",
    "spacer_hit_count",
    "repeat_predicted_cas_subtype",
    "repeat_prediction_confidence",
    "top_pam_supported_subtype",
    "top_pam_support_count",
    "repeat_predicted_subtype_pam_support_count",
    "pam_supported_subtypes",
    "pam_subtype_support_counts",
    "repeat_pam_subtype_agreement",
]

EVIDENCE_MATRIX_COLUMNS = [
    "bacterium",
    "phage",
    "crispr_targeting_score",
    "experimental_pam_weighted_score",
    "hypothetical_resistance_score",
    "spacer_hits",
    "unique_matching_spacers",
    "best_identity_percent",
    "best_coverage_percent",
    "pam_supported_hits",
    "pam_evaluated_hits",
    "pam_support_level",
    "best_pam_compatibility_score",
    "mean_pam_compatibility_score",
    "seed_evaluated_hits",
    "best_seed_mismatches",
    "current_evidence_level",
    "evidence_summary",
    "interpretation",
]


@dataclass(frozen=True)
class InitialRunSummary:
    bacterial_file_count: int
    bacterial_sequence_count: int
    phage_file_count: int
    phage_sequence_count: int

    def stage_table(self) -> pd.DataFrame:
        return pd.DataFrame(
            [
                {
                    "stage": "FASTA parsing",
                    "status": "ready",
                    "notes": "Multiple bacterial and phage FASTA files are accepted.",
                },
                {
                    "stage": "CRISPR array detection",
                    "status": "ready: exact-repeat MVP",
                    "notes": "Finds candidate direct-repeat arrays with plausible spacer lengths.",
                },
                {
                    "stage": "Spacer-phage matching",
                    "status": "ready: exact-match MVP",
                    "notes": "Crosses all extracted spacers against all uploaded phage genomes.",
                },
                {
                    "stage": "Cas type classifier",
                    "status": "ready: repeat/array ExtraTrees model when artifact is present",
                    "notes": "Predicts likely subtype from FASTA-derived repeat and array features.",
                },
                {
                    "stage": "PAM analysis",
                    "status": "ready: subtype-aware curated rule subset",
                    "notes": "Evaluates protospacer flanks when a supported subtype rule is available.",
                },
                {
                    "stage": "CRISPR targeting evidence scoring",
                    "status": "ready",
                    "notes": "Produces a bacteria-by-phage spacer-targeting evidence matrix.",
                },
            ]
        )


def build_initial_run_summary(
    bacteria_records: list[FastaRecord],
    phage_records: list[FastaRecord],
) -> InitialRunSummary:
    return InitialRunSummary(
        bacterial_file_count=len({record.source_file for record in bacteria_records}),
        bacterial_sequence_count=len(bacteria_records),
        phage_file_count=len({record.source_file for record in phage_records}),
        phage_sequence_count=len(phage_records),
    )


def detect_arrays_for_records(
    records: list[FastaRecord],
    method: DetectionMethod = "internal",
    progress_callback: Callable[[int, int, FastaRecord], None] | None = None,
) -> list[CrisprArray]:
    if method == "minced":
        return detect_arrays_with_minced(records, progress_callback=progress_callback)
    if method != "internal":
        raise ValueError(f"Unknown CRISPR detection method: {method}")

    arrays: list[CrisprArray] = []
    total = len(records)
    repeat_scan_steps = DEFAULT_MAX_REPEAT_LENGTH - DEFAULT_MIN_REPEAT_LENGTH + 1
    total_scan_steps = max(total * repeat_scan_steps, 1)
    for index, record in enumerate(records, start=1):
        genome_id = record.source_file

        def update_internal_scan_progress(repeat_step: int, _repeat_total: int) -> None:
            if progress_callback:
                completed_steps = ((index - 1) * repeat_scan_steps) + repeat_step
                progress_callback(completed_steps, total_scan_steps, record)

        arrays.extend(
            detect_crispr_arrays(
                sequence=record.sequence,
                genome_id=genome_id,
                contig_id=record.record_id,
                scan_progress_callback=update_internal_scan_progress,
            )
        )
    return arrays


def find_spacer_hits_for_records(
    crispr_arrays: list[CrisprArray],
    phage_records: list[FastaRecord],
    method: MatchingMethod = "internal",
    blast_min_identity: float = 0.9,
    blast_min_coverage: float = 0.95,
    blast_require_full_query: bool = True,
) -> list[SpacerHit]:
    if method == "blast":
        return find_spacer_hits_with_blast(
            crispr_arrays,
            phage_records,
            min_identity=blast_min_identity,
            min_coverage=blast_min_coverage,
            require_full_query=blast_require_full_query,
        )
    if method != "internal":
        raise ValueError(f"Unknown spacer matching method: {method}")
    return find_spacer_hits(crispr_arrays, phage_records)


def annotate_spacer_hits_with_pam(
    hits: list[SpacerHit],
    pam_rules_by_array: dict[str, str] | None = None,
    pam_rules_by_bacterium: dict[str, str] | None = None,
    cas_predictions_by_array: dict[str, ArrayCasPrediction] | None = None,
    default_pam_rule: str | None = None,
    seed_length: int = 8,
) -> list[SpacerHit]:
    annotated_hits = []
    for hit in hits:
        cas_prediction = (cas_predictions_by_array or {}).get(hit.array_id)
        pam_rule = _select_pam_rule(
            hit=hit,
            cas_prediction=cas_prediction,
            pam_rules_by_array=pam_rules_by_array or {},
            pam_rules_by_bacterium=pam_rules_by_bacterium or {},
            default_pam_rule=default_pam_rule,
        )
        evaluation = evaluate_pam_rule(
            protospacer_5p_flank=hit.protospacer_5p_flank,
            protospacer_3p_flank=hit.protospacer_3p_flank,
            pam_rule=pam_rule,
            genomic_upstream_flank=hit.genomic_upstream_flank,
            genomic_downstream_flank=hit.genomic_downstream_flank,
        )
        seed_summary = summarize_seed_mismatches(
            spacer_sequence=hit.spacer_sequence,
            protospacer_sequence=_strand_oriented_protospacer(hit),
            pam_rule=evaluation.pam_rule,
            seed_length=seed_length,
        )
        annotated_hits.append(
            replace(
                hit,
                pam_rule=evaluation.pam_rule,
                pam_rule_source=_pam_rule_source(
                    hit=hit,
                    cas_prediction=cas_prediction,
                    pam_rules_by_array=pam_rules_by_array or {},
                    pam_rules_by_bacterium=pam_rules_by_bacterium or {},
                    default_pam_rule=default_pam_rule,
                ),
                predicted_cas_subtype=cas_prediction.cas_subtype if cas_prediction else "",
                cas_subtype_confidence=(
                    cas_prediction.cas_subtype_confidence if cas_prediction else None
                ),
                cas_subtype_prediction_method=(
                    cas_prediction.prediction_method if cas_prediction else ""
                ),
                pam_sequence=evaluation.pam_sequence,
                pam_match=evaluation.pam_match,
                pam_support_level=evaluation.pam_support_level,
                pam_compatibility_score=evaluation.compatibility_score,
                pam_offset_from_protospacer=evaluation.pam_offset_from_protospacer,
                seed_region=seed_summary.seed_region if seed_summary else "",
                seed_length=seed_summary.seed_length if seed_summary else None,
                seed_mismatches=seed_summary.seed_mismatches if seed_summary else None,
                seed_mismatch_positions=(
                    seed_summary.seed_mismatch_positions if seed_summary else ""
                ),
            )
        )
    return annotated_hits


def _strand_oriented_protospacer(hit: SpacerHit) -> str:
    if hit.aligned_protospacer_sequence:
        return hit.aligned_protospacer_sequence
    if hit.strand == "-":
        return reverse_complement(hit.protospacer_sequence)
    return hit.protospacer_sequence


def _select_pam_rule(
    hit: SpacerHit,
    cas_prediction: ArrayCasPrediction | None,
    pam_rules_by_array: dict[str, str],
    pam_rules_by_bacterium: dict[str, str],
    default_pam_rule: str | None,
) -> str | None:
    if hit.array_id in pam_rules_by_array:
        return pam_rules_by_array[hit.array_id]
    if hit.bacterium_id in pam_rules_by_bacterium:
        return pam_rules_by_bacterium[hit.bacterium_id]
    if cas_prediction and cas_prediction.pam_rule:
        return cas_prediction.pam_rule
    return default_pam_rule


def _pam_rule_source(
    hit: SpacerHit,
    cas_prediction: ArrayCasPrediction | None,
    pam_rules_by_array: dict[str, str],
    pam_rules_by_bacterium: dict[str, str],
    default_pam_rule: str | None,
) -> str:
    if hit.array_id in pam_rules_by_array:
        return "array_override"
    if hit.bacterium_id in pam_rules_by_bacterium:
        return "bacterium_override"
    if cas_prediction and cas_prediction.pam_rule:
        return cas_prediction.pam_rule_source
    if default_pam_rule:
        return "manual_default"
    return cas_prediction.pam_rule_source if cas_prediction else "not_available"


def summarize_crispr_arrays(arrays: list[CrisprArray]) -> pd.DataFrame:
    rows = [
        {
            "array_id": array.array_id,
            "bacterium": array.genome_id,
            "contig": array.contig_id,
            "start": array.start,
            "end": array.end,
            "repeat_length": array.repeat_length,
            "repeat_count": array.repeat_count,
            "spacer_count": array.spacer_count,
            "mean_spacer_length": round(array.mean_spacer_length, 2),
            "repeat_consensus": array.repeat_consensus,
        }
        for array in arrays
    ]
    return pd.DataFrame(rows, columns=CRISPR_ARRAY_COLUMNS)


def summarize_spacers(arrays: list[CrisprArray]) -> pd.DataFrame:
    rows = []
    for array in arrays:
        for index, spacer in enumerate(array.spacers, start=1):
            rows.append(
                {
                    "array_id": array.array_id,
                    "spacer_id": f"{array.array_id}|spacer_{index}",
                    "bacterium": array.genome_id,
                    "contig": array.contig_id,
                    "spacer_index": index,
                    "spacer_length": len(spacer),
                    "spacer_sequence": spacer,
                }
            )
    return pd.DataFrame(rows, columns=SPACER_COLUMNS)


def summarize_spacer_hits(hits: list[SpacerHit]) -> pd.DataFrame:
    rows = [
        {
            "bacterium": hit.bacterium_id,
            "phage": hit.phage_id,
            "array_id": hit.array_id,
            "spacer_id": hit.spacer_id,
            "phage_contig": hit.phage_contig_id,
            "start": hit.start,
            "end": hit.end,
            "strand": hit.strand,
            "identity_percent": round(hit.identity * 100, 2),
            "mismatches": hit.mismatches,
            "alignment_length": hit.alignment_length,
            "spacer_length": hit.spacer_length,
            "coverage_percent": round(hit.coverage * 100, 2),
            "evalue": hit.evalue,
            "bitscore": hit.bitscore,
            "spacer_sequence": hit.spacer_sequence,
            "aligned_spacer_sequence": hit.aligned_spacer_sequence,
            "aligned_protospacer_sequence": hit.aligned_protospacer_sequence,
            "protospacer_sequence": hit.protospacer_sequence,
            "protospacer_5p_flank": hit.protospacer_5p_flank,
            "protospacer_3p_flank": hit.protospacer_3p_flank,
            "genomic_upstream_flank": hit.genomic_upstream_flank,
            "genomic_downstream_flank": hit.genomic_downstream_flank,
            "predicted_cas_subtype": hit.predicted_cas_subtype,
            "cas_subtype_confidence": hit.cas_subtype_confidence,
            "cas_subtype_prediction_method": hit.cas_subtype_prediction_method,
            "pam_rule": hit.pam_rule,
            "pam_rule_source": hit.pam_rule_source,
            "pam_sequence": hit.pam_sequence,
            "pam_match": hit.pam_match,
            "pam_support_level": hit.pam_support_level,
            "pam_compatibility_score": hit.pam_compatibility_score,
            "pam_offset_from_protospacer": hit.pam_offset_from_protospacer,
            "seed_region": hit.seed_region,
            "seed_length": hit.seed_length,
            "seed_mismatches": hit.seed_mismatches,
            "seed_mismatch_positions": hit.seed_mismatch_positions,
        }
        for hit in hits
    ]
    return pd.DataFrame(rows, columns=SPACER_HIT_COLUMNS)


def summarize_pam_subtype_support(
    hits: list[SpacerHit],
    pam_rules_by_subtype: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Explore subtype support from observed protospacer flanks.

    This is diagnostic only. It does not change repeat-based Cas subtype
    predictions or targeting scores.
    """
    rules = pam_rules_by_subtype or CURATED_PAM_RULES_BY_SUBTYPE
    rows = []
    for array_id in sorted({hit.array_id for hit in hits}):
        array_hits = [hit for hit in hits if hit.array_id == array_id]
        predicted_subtype = _first_nonempty(hit.predicted_cas_subtype for hit in array_hits)
        predicted_confidence = _first_present(
            hit.cas_subtype_confidence for hit in array_hits if hit.cas_subtype_confidence is not None
        )
        support_counts = _pam_subtype_support_counts(array_hits, rules)
        supported_subtypes = [
            subtype
            for subtype, count in sorted(
                support_counts.items(),
                key=lambda item: (-item[1], item[0]),
            )
            if count > 0
        ]
        top_subtype = supported_subtypes[0] if supported_subtypes else ""
        top_count = support_counts.get(top_subtype, 0) if top_subtype else 0
        predicted_count = support_counts.get(predicted_subtype, 0) if predicted_subtype else 0
        rows.append(
            {
                "array_id": array_id,
                "bacterium": _first_nonempty(hit.bacterium_id for hit in array_hits),
                "spacer_hit_count": len(array_hits),
                "repeat_predicted_cas_subtype": predicted_subtype,
                "repeat_prediction_confidence": predicted_confidence,
                "top_pam_supported_subtype": top_subtype,
                "top_pam_support_count": top_count,
                "repeat_predicted_subtype_pam_support_count": predicted_count,
                "pam_supported_subtypes": ";".join(supported_subtypes),
                "pam_subtype_support_counts": ";".join(
                    f"{subtype}:{support_counts[subtype]}"
                    for subtype in sorted(support_counts)
                    if support_counts[subtype] > 0
                ),
                "repeat_pam_subtype_agreement": _repeat_pam_agreement(
                    predicted_subtype=predicted_subtype,
                    top_pam_subtype=top_subtype,
                    predicted_support_count=predicted_count,
                    top_support_count=top_count,
                ),
            }
        )
    return pd.DataFrame(rows, columns=PAM_SUBTYPE_SUPPORT_COLUMNS)


def build_crispr_targeting_evidence_matrix(
    bacteria_records: list[FastaRecord],
    phage_records: list[FastaRecord],
    hits: list[SpacerHit],
) -> pd.DataFrame:
    bacteria_ids = sorted({record.source_file for record in bacteria_records})
    phage_ids = sorted({record.source_file for record in phage_records})
    rows = []
    for bacterium_id in bacteria_ids:
        for phage_id in phage_ids:
            pair_hits = [
                hit for hit in hits if hit.bacterium_id == bacterium_id and hit.phage_id == phage_id
            ]
            unique_spacers = {hit.spacer_id for hit in pair_hits}
            best_identity = max((hit.identity for hit in pair_hits), default=0.0)
            best_coverage = max((hit.coverage for hit in pair_hits), default=0.0)
            pam_supported_hits = sum(hit.pam_match is True for hit in pair_hits)
            pam_evaluated_hits = sum(hit.pam_match is not None for hit in pair_hits)
            pam_support_level = _pair_pam_support_level(pair_hits)
            pam_scores = [
                hit.pam_compatibility_score
                for hit in pair_hits
                if hit.pam_compatibility_score is not None
            ]
            seed_evaluated_hits = sum(hit.seed_mismatches is not None for hit in pair_hits)
            best_seed_mismatches = min(
                (
                    hit.seed_mismatches
                    for hit in pair_hits
                    if hit.seed_mismatches is not None
                ),
                default=None,
            )
            targeting_score = score_crispr_targeting_evidence(pair_hits)
            rows.append(
                {
                    "bacterium": bacterium_id,
                    "phage": phage_id,
                    "crispr_targeting_score": targeting_score.score,
                    "experimental_pam_weighted_score": score_experimental_pam_weighted_evidence(
                        pair_hits
                    ),
                    "hypothetical_resistance_score": targeting_score.score,
                    "spacer_hits": len(pair_hits),
                    "unique_matching_spacers": len(unique_spacers),
                    "best_identity_percent": round(best_identity * 100, 2),
                    "best_coverage_percent": round(best_coverage * 100, 2),
                    "pam_supported_hits": pam_supported_hits,
                    "pam_evaluated_hits": pam_evaluated_hits,
                    "pam_support_level": pam_support_level,
                    "best_pam_compatibility_score": (
                        round(max(pam_scores), 6) if pam_scores else None
                    ),
                    "mean_pam_compatibility_score": (
                        round(sum(pam_scores) / len(pam_scores), 6) if pam_scores else None
                    ),
                    "seed_evaluated_hits": seed_evaluated_hits,
                    "best_seed_mismatches": best_seed_mismatches,
                    "current_evidence_level": targeting_score.evidence_level,
                    "evidence_summary": targeting_score.evidence_summary,
                    "interpretation": targeting_score.interpretation,
                }
            )
    return pd.DataFrame(rows, columns=EVIDENCE_MATRIX_COLUMNS)


def build_resistance_evidence_matrix(
    bacteria_records: list[FastaRecord],
    phage_records: list[FastaRecord],
    hits: list[SpacerHit],
) -> pd.DataFrame:
    """Backward-compatible alias for older callers and saved-run tooling."""
    return build_crispr_targeting_evidence_matrix(bacteria_records, phage_records, hits)


def build_exact_match_heatmap(evidence_matrix: pd.DataFrame) -> pd.DataFrame:
    if evidence_matrix.empty:
        return pd.DataFrame()
    return evidence_matrix.pivot(
        index="bacterium",
        columns="phage",
        values="unique_matching_spacers",
    ).fillna(0).astype(int)


def _pair_pam_support_level(pair_hits: list[SpacerHit]) -> str:
    if not pair_hits:
        return "not_applicable"
    if any(hit.pam_match is True for hit in pair_hits):
        return "compatible"
    if any(hit.pam_match is False for hit in pair_hits):
        return "not_supported"
    if any(hit.pam_support_level == "insufficient_flank" for hit in pair_hits):
        return "insufficient_flank"
    if any(hit.pam_support_level == "invalid_rule" for hit in pair_hits):
        return "invalid_rule"
    return "not_evaluated"


def _pam_subtype_support_counts(
    hits: list[SpacerHit],
    pam_rules_by_subtype: dict[str, str],
) -> dict[str, int]:
    counts = {subtype: 0 for subtype in pam_rules_by_subtype}
    for hit in hits:
        for subtype, rule in pam_rules_by_subtype.items():
            evaluation = evaluate_pam_rule(
                protospacer_5p_flank=hit.protospacer_5p_flank,
                protospacer_3p_flank=hit.protospacer_3p_flank,
                genomic_upstream_flank=hit.genomic_upstream_flank,
                genomic_downstream_flank=hit.genomic_downstream_flank,
                pam_rule=rule,
            )
            if evaluation.pam_match is True:
                counts[subtype] += 1
    return counts


def _repeat_pam_agreement(
    predicted_subtype: str,
    top_pam_subtype: str,
    predicted_support_count: int,
    top_support_count: int,
) -> str:
    if not top_pam_subtype:
        return "no_pam_subtype_support"
    if not predicted_subtype:
        return "no_repeat_subtype_prediction"
    if predicted_support_count <= 0:
        return "repeat_prediction_not_pam_supported"
    if predicted_subtype == top_pam_subtype:
        return "agrees_top_subtype"
    if predicted_support_count == top_support_count:
        return "agrees_tied_subtype"
    return "conflicts_with_top_subtype"


def _first_nonempty(values) -> str:
    for value in values:
        text = str(value or "").strip()
        if text:
            return text
    return ""


def _first_present(values):
    for value in values:
        return value
    return None
