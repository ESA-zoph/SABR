from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path

import pandas as pd

from crispr_phage_predictor.crispr import CrisprArray
from crispr_phage_predictor.io import FastaRecord
from crispr_phage_predictor.matching import SpacerHit
from crispr_phage_predictor.ml.model_artifact import model_artifact_metadata
from crispr_phage_predictor.pipeline import (
    summarize_crispr_arrays,
    summarize_pam_subtype_support,
    summarize_spacer_hits,
    summarize_spacers,
)


DEFAULT_OUTPUT_ROOT = Path("outputs") / "runs"


@dataclass(frozen=True)
class AnalysisRunMetadata:
    run_id: str
    created_at: str
    detection_method: str
    matching_method: str
    detection_backend_detail: str | None
    blast_min_identity: float | None
    blast_min_coverage: float | None
    blast_require_full_query: bool | None
    pam_mode: str | None
    pam_rule: str | None
    cas_prediction_count: int | None
    cas_model_artifact: dict | None
    seed_length: int | None
    detection_elapsed_seconds: float | None
    matching_elapsed_seconds: float | None
    total_elapsed_seconds: float | None
    bacterial_file_count: int
    bacterial_sequence_count: int
    phage_file_count: int
    phage_sequence_count: int
    candidate_array_count: int
    extracted_spacer_count: int
    spacer_hit_count: int
    bacterial_duplicate_record_count: int
    phage_duplicate_record_count: int


def save_analysis_run(
    bacteria_records: list[FastaRecord],
    phage_records: list[FastaRecord],
    crispr_arrays: list[CrisprArray],
    spacer_hits: list[SpacerHit],
    evidence_matrix: pd.DataFrame,
    heatmap: pd.DataFrame,
    detection_method: str,
    matching_method: str,
    detection_backend_detail: str | None = None,
    blast_min_identity: float | None = None,
    blast_min_coverage: float | None = None,
    blast_require_full_query: bool | None = None,
    pam_mode: str | None = None,
    pam_rule: str | None = None,
    cas_prediction_count: int | None = None,
    cas_model_artifact: dict | None = None,
    seed_length: int | None = None,
    bacterial_duplicate_record_count: int = 0,
    phage_duplicate_record_count: int = 0,
    detection_elapsed_seconds: float | None = None,
    matching_elapsed_seconds: float | None = None,
    total_elapsed_seconds: float | None = None,
    output_root: Path = DEFAULT_OUTPUT_ROOT,
) -> Path:
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = output_root / run_id
    run_dir.mkdir(parents=True, exist_ok=False)

    metadata = AnalysisRunMetadata(
        run_id=run_id,
        created_at=datetime.now().isoformat(timespec="seconds"),
        detection_method=detection_method,
        matching_method=matching_method,
        detection_backend_detail=detection_backend_detail,
        blast_min_identity=blast_min_identity,
        blast_min_coverage=blast_min_coverage,
        blast_require_full_query=blast_require_full_query,
        pam_mode=pam_mode,
        pam_rule=pam_rule,
        cas_prediction_count=cas_prediction_count,
        cas_model_artifact=(
            cas_model_artifact
            if cas_model_artifact is not None
            else model_artifact_metadata()
        ),
        seed_length=seed_length,
        detection_elapsed_seconds=_round_optional(detection_elapsed_seconds),
        matching_elapsed_seconds=_round_optional(matching_elapsed_seconds),
        total_elapsed_seconds=_round_optional(total_elapsed_seconds),
        bacterial_file_count=len({record.source_file for record in bacteria_records}),
        bacterial_sequence_count=len(bacteria_records),
        phage_file_count=len({record.source_file for record in phage_records}),
        phage_sequence_count=len(phage_records),
        candidate_array_count=len(crispr_arrays),
        extracted_spacer_count=sum(array.spacer_count for array in crispr_arrays),
        spacer_hit_count=len(spacer_hits),
        bacterial_duplicate_record_count=bacterial_duplicate_record_count,
        phage_duplicate_record_count=phage_duplicate_record_count,
    )

    _write_json(run_dir / "run_metadata.json", asdict(metadata))
    _write_records_summary(run_dir / "bacterial_records.csv", bacteria_records)
    _write_records_summary(run_dir / "phage_records.csv", phage_records)
    _write_table(run_dir / "crispr_arrays.csv", summarize_crispr_arrays(crispr_arrays))
    _write_table(run_dir / "spacers.csv", summarize_spacers(crispr_arrays))
    _write_table(run_dir / "spacer_hits.csv", summarize_spacer_hits(spacer_hits))
    _write_table(run_dir / "pam_subtype_support.csv", summarize_pam_subtype_support(spacer_hits))
    _write_table(run_dir / "evidence_matrix.csv", evidence_matrix)
    _write_table(run_dir / "heatmap.csv", heatmap.reset_index())
    return run_dir


def _write_json(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, indent=2), encoding="utf-8")


def _round_optional(value: float | None) -> float | None:
    if value is None:
        return None
    return round(value, 3)


def _write_records_summary(path: Path, records: list[FastaRecord]) -> None:
    rows = [
        {
            "source_file": record.source_file,
            "record_id": record.record_id,
            "description": record.description,
            "length_bp": record.length,
            "gc_percent": round(record.gc_fraction * 100, 2),
            "accession": record.accession,
            "sequence_hash": record.sequence_hash,
        }
        for record in records
    ]
    _write_table(path, pd.DataFrame(rows))


def _write_table(path: Path, table: pd.DataFrame) -> None:
    table.to_csv(path, index=False)
