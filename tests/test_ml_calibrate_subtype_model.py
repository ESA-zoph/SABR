from __future__ import annotations

import pandas as pd

from crispr_phage_predictor.ml.calibrate_subtype_model import (
    calibrate_extra_trees_subtype_model,
)


def _row(genome_id: str, repeat: str, subtype: str) -> dict[str, object]:
    return {
        "source": "test",
        "genome_id": genome_id,
        "organism": genome_id,
        "taxonomy": "",
        "assembly_level": "",
        "contig_id": genome_id,
        "array_start": "",
        "array_end": "",
        "repeat_sequence": repeat,
        "repeat_length": len(repeat),
        "spacer_count": 3,
        "mean_spacer_length": 32.0,
        "cas_type": "Type I",
        "cas_subtype": subtype,
        "label_source": "test",
        "label_confidence": "high",
        "pam_rule": "",
    }


def test_calibrate_extra_trees_subtype_model_returns_expected_tables(tmp_path):
    path = tmp_path / "training.csv"
    pd.DataFrame(
        [
            _row("g1", "AAAAAAAAAAAAAAAAAAAAAAAAAAAA", "I-E"),
            _row("g2", "AAAAAAAATAAAAAAAAAAAAAAAAAAA", "I-E"),
            _row("g3", "AAAAAAAAGAAAAAAAAAAAAAAAAAAA", "I-E"),
            _row("g4", "AAAAAAAACAAAAAAAAAAAAAAAAAAA", "I-E"),
            _row("g5", "CCCCCCCCCCCCCCCCCCCCCCCCCCCC", "I-F"),
            _row("g6", "CCCCCCCCGCCCCCCCCCCCCCCCCCCC", "I-F"),
            _row("g7", "CCCCCCCCTCCCCCCCCCCCCCCCCCCC", "I-F"),
            _row("g8", "CCCCCCCCACCCCCCCCCCCCCCCCCCC", "I-F"),
        ]
    ).to_csv(path, index=False)

    outputs = calibrate_extra_trees_subtype_model(
        path,
        split_strategy="row_random",
        min_class_count=1,
        n_bins=5,
    )

    assert set(outputs) == {
        "summary",
        "predictions",
        "confidence_bins",
        "subtype_confidence",
        "accuracy_by_threshold",
    }
    assert len(outputs["confidence_bins"]) == 5
    assert outputs["summary"].iloc[0]["test_rows"] == len(outputs["predictions"])
