from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ResistanceScore:
    bacterium_id: str
    phage_id: str
    score: float
    interpretation: str
    evidence_summary: str


def score_resistance_likelihood() -> list[ResistanceScore]:
    raise NotImplementedError(
        "Resistance likelihood scoring will be implemented after spacer matching and PAM analysis."
    )
