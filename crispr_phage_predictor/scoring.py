from __future__ import annotations

from dataclasses import dataclass
from math import log2

from crispr_phage_predictor.matching import SpacerHit


@dataclass(frozen=True)
class TargetingScore:
    score: float
    evidence_level: str
    interpretation: str
    evidence_summary: str


def score_crispr_targeting_evidence(pair_hits: list[SpacerHit]) -> TargetingScore:
    """Score candidate CRISPR-mediated targeting evidence on a 0-100 scale.

    The score is a transparent spacer-targeting evidence summary, not a
    biological resistance probability. It intentionally rewards independent
    spacer evidence, high identity/coverage, PAM/PFS compatibility, low seed
    mismatch burden, and confident subtype prediction when available.
    """
    if not pair_hits:
        return TargetingScore(
            score=0.0,
            evidence_level="no spacer-match evidence",
            interpretation="No spacer-protospacer match detected by the selected pipeline.",
            evidence_summary="No spacer hits.",
        )

    unique_spacers = {hit.spacer_id for hit in pair_hits}
    best_identity = max(hit.identity for hit in pair_hits)
    best_coverage = max(hit.coverage for hit in pair_hits)
    pam_supported_hits = sum(hit.pam_match is True for hit in pair_hits)
    pam_evaluated_hits = sum(hit.pam_match is not None for hit in pair_hits)
    best_seed_mismatches = min(
        (hit.seed_mismatches for hit in pair_hits if hit.seed_mismatches is not None),
        default=None,
    )
    best_subtype_confidence = max(
        (
            hit.cas_subtype_confidence
            for hit in pair_hits
            if hit.cas_subtype_confidence is not None
        ),
        default=None,
    )

    spacer_component = min(35.0, log2(len(unique_spacers) + 1) * 18.0)
    identity_component = max(0.0, min(25.0, best_identity * best_coverage * 25.0))
    pam_component = _pam_component(pair_hits, pam_supported_hits, pam_evaluated_hits)
    seed_component = _seed_component(best_seed_mismatches, pam_supported_hits)
    subtype_component = (
        max(0.0, min(10.0, best_subtype_confidence * 10.0))
        if best_subtype_confidence is not None
        else 0.0
    )
    raw_score = (
        spacer_component
        + identity_component
        + pam_component
        + seed_component
        + subtype_component
    )
    score_cap = _score_cap(pam_supported_hits, pam_evaluated_hits)
    score = round(min(score_cap, max(0.0, raw_score)), 2)
    evidence_level = _score_level(score)
    return TargetingScore(
        score=score,
        evidence_level=evidence_level,
        interpretation=_interpret_score(
            evidence_level=evidence_level,
            pam_supported_hits=pam_supported_hits,
            pam_evaluated_hits=pam_evaluated_hits,
            best_seed_mismatches=best_seed_mismatches,
        ),
        evidence_summary=(
            f"{len(unique_spacers)} unique spacer(s); best identity "
            f"{best_identity * 100:.2f}%; best coverage {best_coverage * 100:.2f}%; "
            f"{pam_supported_hits}/{pam_evaluated_hits} PAM/PFS-evaluated hit(s) supported; "
            f"best seed mismatches: {_format_optional_int(best_seed_mismatches)}."
        ),
    )


ResistanceScore = TargetingScore


def score_resistance_likelihood(pair_hits: list[SpacerHit]) -> TargetingScore:
    """Backward-compatible alias for older callers.

    New code should use :func:`score_crispr_targeting_evidence` to avoid
    implying that spacer/PAM evidence alone proves biological resistance.
    """
    return score_crispr_targeting_evidence(pair_hits)


def score_experimental_pam_weighted_evidence(pair_hits: list[SpacerHit]) -> float:
    """Experimental score using probabilistic PAM/PFS compatibility.

    This is a diagnostic comparator. It is intentionally not used as the
    production targeting score until benchmark calibration supports it.
    """
    if not pair_hits:
        return 0.0

    unique_spacers = {hit.spacer_id for hit in pair_hits}
    best_identity = max(hit.identity for hit in pair_hits)
    best_coverage = max(hit.coverage for hit in pair_hits)
    pam_scores = [
        hit.pam_compatibility_score
        for hit in pair_hits
        if hit.pam_compatibility_score is not None
    ]
    best_pam_score = max(pam_scores) if pam_scores else None
    mean_pam_score = sum(pam_scores) / len(pam_scores) if pam_scores else None
    best_seed_mismatches = min(
        (hit.seed_mismatches for hit in pair_hits if hit.seed_mismatches is not None),
        default=None,
    )
    best_subtype_confidence = max(
        (
            hit.cas_subtype_confidence
            for hit in pair_hits
            if hit.cas_subtype_confidence is not None
        ),
        default=None,
    )

    spacer_component = min(35.0, log2(len(unique_spacers) + 1) * 18.0)
    identity_component = max(0.0, min(25.0, best_identity * best_coverage * 25.0))
    pam_component = _experimental_pam_component(best_pam_score, mean_pam_score)
    seed_component = _seed_component(
        best_seed_mismatches=best_seed_mismatches,
        pam_supported_hits=1 if (best_pam_score or 0.0) >= 1.0 else 0,
    )
    subtype_component = (
        max(0.0, min(10.0, best_subtype_confidence * 10.0))
        if best_subtype_confidence is not None
        else 0.0
    )
    raw_score = (
        spacer_component
        + identity_component
        + pam_component
        + seed_component
        + subtype_component
    )
    return round(min(_experimental_score_cap(best_pam_score), max(0.0, raw_score)), 2)


def _pam_component(
    pair_hits: list[SpacerHit],
    pam_supported_hits: int,
    pam_evaluated_hits: int,
) -> float:
    if pam_supported_hits:
        return 20.0
    if pam_evaluated_hits:
        return -25.0
    if any(hit.pam_support_level in {"insufficient_flank", "invalid_rule"} for hit in pair_hits):
        return 0.0
    return 0.0


def _experimental_pam_component(
    best_pam_score: float | None,
    mean_pam_score: float | None,
) -> float:
    if best_pam_score is None:
        return 0.0
    mean_score = mean_pam_score if mean_pam_score is not None else best_pam_score
    return round((20.0 * best_pam_score) + (10.0 * mean_score) - 15.0, 6)


def _score_cap(pam_supported_hits: int, pam_evaluated_hits: int) -> float:
    if pam_supported_hits:
        return 100.0
    if pam_evaluated_hits:
        return 39.0
    return 100.0


def _experimental_score_cap(best_pam_score: float | None) -> float:
    if best_pam_score is None:
        return 100.0
    if best_pam_score >= 1.0:
        return 100.0
    if best_pam_score >= 0.75:
        return 70.0
    if best_pam_score >= 0.5:
        return 55.0
    return 39.0


def _seed_component(best_seed_mismatches: int | None, pam_supported_hits: int) -> float:
    if not pam_supported_hits:
        return 0.0
    if best_seed_mismatches is None:
        return 0.0
    if best_seed_mismatches == 0:
        return 10.0
    if best_seed_mismatches <= 2:
        return 5.0
    return -5.0


def _score_level(score: float) -> str:
    if score >= 75:
        return "strong candidate CRISPR targeting evidence"
    if score >= 50:
        return "moderate candidate CRISPR targeting evidence"
    if score > 0:
        return "weak candidate CRISPR targeting evidence"
    return "no spacer-match evidence"


def _interpret_score(
    evidence_level: str,
    pam_supported_hits: int,
    pam_evaluated_hits: int,
    best_seed_mismatches: int | None,
) -> str:
    if evidence_level == "no spacer-match evidence":
        return "No spacer-protospacer match detected by the selected pipeline."
    caveat = (
        "This is a hypothetical CRISPR targeting score, not confirmed resistance; "
        "Cas function, expression, phage escape, anti-CRISPR genes, and assembly quality remain unresolved."
    )
    if pam_supported_hits:
        return (
            f"{evidence_level} with PAM/PFS support in at least one hit. "
            f"Best seed mismatches: {_format_optional_int(best_seed_mismatches)}. {caveat}"
        )
    if pam_evaluated_hits:
        return (
            f"{evidence_level}, but evaluated flanks did not support the supplied or predicted PAM/PFS rule. "
            f"Best seed mismatches: {_format_optional_int(best_seed_mismatches)}. {caveat}"
        )
    return (
        f"{evidence_level}; PAM/PFS compatibility was not evaluated or no curated rule was available. "
        f"Best seed mismatches: {_format_optional_int(best_seed_mismatches)}. {caveat}"
    )


def _format_optional_int(value: int | None) -> str:
    return "not evaluated" if value is None else str(value)
