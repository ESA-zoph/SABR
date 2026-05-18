from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from crispr_phage_predictor.crispr import CrisprArray, detect_crispr_arrays
from crispr_phage_predictor.io import FastaRecord
from crispr_phage_predictor.matching import SpacerHit


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
                    "status": "planned",
                    "notes": "Will predict type/subtype from repeat, array, and cas-gene features.",
                },
                {
                    "stage": "PAM analysis",
                    "status": "planned",
                    "notes": "Will evaluate protospacer flanks using predicted system type.",
                },
                {
                    "stage": "Resistance likelihood scoring",
                    "status": "planned",
                    "notes": "Will produce a bacteria-by-phage evidence matrix.",
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


def detect_arrays_for_records(records: list[FastaRecord]) -> list[CrisprArray]:
    arrays: list[CrisprArray] = []
    for record in records:
        genome_id = record.source_file
        arrays.extend(
            detect_crispr_arrays(
                sequence=record.sequence,
                genome_id=genome_id,
                contig_id=record.record_id,
            )
        )
    return arrays


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
    return pd.DataFrame(rows)


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
    return pd.DataFrame(rows)


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
            "spacer_sequence": hit.spacer_sequence,
            "protospacer_sequence": hit.protospacer_sequence,
        }
        for hit in hits
    ]
    return pd.DataFrame(rows)


def build_resistance_evidence_matrix(
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
            rows.append(
                {
                    "bacterium": bacterium_id,
                    "phage": phage_id,
                    "exact_spacer_hits": len(pair_hits),
                    "unique_matching_spacers": len(unique_spacers),
                    "best_identity_percent": 100.0 if pair_hits else 0.0,
                    "current_evidence_level": _evidence_level(len(unique_spacers)),
                    "interpretation": _interpret_exact_match_evidence(len(unique_spacers)),
                }
            )
    return pd.DataFrame(rows)


def build_exact_match_heatmap(evidence_matrix: pd.DataFrame) -> pd.DataFrame:
    if evidence_matrix.empty:
        return pd.DataFrame()
    return evidence_matrix.pivot(
        index="bacterium",
        columns="phage",
        values="unique_matching_spacers",
    ).fillna(0).astype(int)


def _evidence_level(unique_spacer_count: int) -> str:
    if unique_spacer_count >= 3:
        return "strong exact-match evidence"
    if unique_spacer_count == 2:
        return "moderate exact-match evidence"
    if unique_spacer_count == 1:
        return "single-spacer exact-match evidence"
    return "no exact-match evidence"


def _interpret_exact_match_evidence(unique_spacer_count: int) -> str:
    if unique_spacer_count:
        return (
            "Candidate CRISPR targeting evidence. PAM, Cas type, seed mismatches, "
            "and cas gene functionality are not evaluated yet."
        )
    return "No exact spacer-protospacer match detected by the current MVP pipeline."
