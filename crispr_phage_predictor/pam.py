from __future__ import annotations

from dataclasses import dataclass


IUPAC_BASES = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "U": {"T"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}


@dataclass(frozen=True)
class PamEvaluation:
    pam_rule: str
    pam_sequence: str
    pam_match: bool | None
    pam_support_level: str
    compatibility_score: float | None = None
    pam_offset_from_protospacer: int | None = None


def evaluate_pam_rule(
    protospacer_5p_flank: str,
    protospacer_3p_flank: str,
    pam_rule: str | None,
    genomic_upstream_flank: str = "",
    genomic_downstream_flank: str = "",
) -> PamEvaluation:
    """Evaluate a simple PAM/PFS rule against protospacer flanks.

    Supported rule format is ``5prime:MOTIF`` or ``3prime:MOTIF`` for
    strand-oriented protospacer flanks. ``genomic_5prime:MOTIF`` and
    ``genomic_3prime:MOTIF`` evaluate the uploaded phage genome coordinate
    flanks directly. Motifs may use IUPAC ambiguity codes, for example
    ``5prime:AWG`` or ``3prime:NGG``.
    """
    if not pam_rule:
        return PamEvaluation("", "", None, "not_evaluated", None, None)

    normalized_rule = pam_rule.strip().upper()
    try:
        side, motif = normalized_rule.split(":", maxsplit=1)
    except ValueError:
        return PamEvaluation(pam_rule, "", None, "invalid_rule", None, None)

    side = _normalize_side(side)
    motif = motif.strip().replace(" ", "")
    if not side or not motif or any(base not in IUPAC_BASES for base in motif):
        return PamEvaluation(pam_rule, "", None, "invalid_rule", None, None)

    flank = _select_flank(
        side=side,
        protospacer_5p_flank=protospacer_5p_flank,
        protospacer_3p_flank=protospacer_3p_flank,
        genomic_upstream_flank=genomic_upstream_flank,
        genomic_downstream_flank=genomic_downstream_flank,
    )
    if len(flank) < len(motif):
        return PamEvaluation(pam_rule, flank, None, "insufficient_flank", None, None)

    adjacent_sequence = (
        flank[-len(motif) :]
        if side in {"5prime", "genomic_5prime"}
        else flank[: len(motif)]
    )
    if _sequence_matches_motif(adjacent_sequence, motif):
        return PamEvaluation(
            pam_rule,
            adjacent_sequence,
            True,
            "compatible",
            _pam_compatibility_score(adjacent_sequence, motif),
            0,
        )

    nearby_sequence, nearby_offset = _nearest_matching_window(
        flank=flank,
        motif=motif,
        side=side,
    )
    if nearby_sequence:
        return PamEvaluation(
            pam_rule,
            nearby_sequence,
            True,
            "compatible_nearby",
            _pam_compatibility_score(nearby_sequence, motif),
            nearby_offset,
        )

    return PamEvaluation(
        pam_rule,
        adjacent_sequence,
        False,
        "not_supported",
        _pam_compatibility_score(adjacent_sequence, motif),
        0,
    )


def _normalize_side(side: str) -> str:
    if side in {"5", "5P", "5PRIME", "UPSTREAM"}:
        return "5prime"
    if side in {"3", "3P", "3PRIME", "DOWNSTREAM"}:
        return "3prime"
    if side in {"GENOMIC_5", "GENOMIC5", "GENOMIC_5P", "GENOMIC_5PRIME", "GENOMIC_UPSTREAM"}:
        return "genomic_5prime"
    if side in {"GENOMIC_3", "GENOMIC3", "GENOMIC_3P", "GENOMIC_3PRIME", "GENOMIC_DOWNSTREAM"}:
        return "genomic_3prime"
    return ""


def _select_flank(
    side: str,
    protospacer_5p_flank: str,
    protospacer_3p_flank: str,
    genomic_upstream_flank: str,
    genomic_downstream_flank: str,
) -> str:
    if side == "5prime":
        return protospacer_5p_flank.upper()
    if side == "3prime":
        return protospacer_3p_flank.upper()
    if side == "genomic_5prime":
        return genomic_upstream_flank.upper()
    if side == "genomic_3prime":
        return genomic_downstream_flank.upper()
    return ""


def _sequence_matches_motif(sequence: str, motif: str) -> bool:
    if len(sequence) != len(motif):
        return False
    for observed, expected in zip(sequence.upper(), motif.upper()):
        allowed = IUPAC_BASES.get(expected)
        if allowed is None or observed not in allowed:
            return False
    return True


def _nearest_matching_window(flank: str, motif: str, side: str) -> tuple[str, int | None]:
    motif_length = len(motif)
    if len(flank) < motif_length:
        return "", None

    if side in {"5prime", "genomic_5prime"}:
        starts = range(len(flank) - motif_length, -1, -1)
        for start in starts:
            candidate = flank[start : start + motif_length]
            if _sequence_matches_motif(candidate, motif):
                offset = len(flank) - (start + motif_length)
                return candidate, offset
        return "", None

    for start in range(0, len(flank) - motif_length + 1):
        candidate = flank[start : start + motif_length]
        if _sequence_matches_motif(candidate, motif):
            return candidate, start
    return "", None


def _pam_compatibility_score(sequence: str, motif: str) -> float:
    if len(sequence) != len(motif) or not motif:
        return 0.0
    position_scores = []
    for observed, expected in zip(sequence.upper(), motif.upper()):
        allowed = IUPAC_BASES.get(expected)
        if not allowed:
            position_scores.append(0.0)
        elif observed in allowed:
            position_scores.append(1.0)
        elif expected == "N":
            position_scores.append(0.25)
        else:
            position_scores.append(0.0)
    return round(sum(position_scores) / len(position_scores), 6)
