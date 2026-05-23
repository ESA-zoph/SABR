from __future__ import annotations

import pandas as pd

from crispr_phage_predictor.ml.evaluate_external_dataset import evaluate_external_dataset


def _row(genome_id: str, repeat: str, subtype: str) -> dict[str, object]:
    return {
        "source": "test",
        "genome_id": genome_id,
        "organism": "",
        "taxonomy": "",
        "assembly_level": "",
        "contig_id": genome_id,
        "array_start": "",
        "array_end": "",
        "repeat_sequence": repeat,
        "repeat_length": len(repeat),
        "spacer_count": 3,
        "mean_spacer_length": 32.0,
        "cas_type": "Type I" if subtype.startswith("I-") else "Type II",
        "cas_subtype": subtype,
        "label_source": "test",
        "label_confidence": "high",
        "pam_rule": "",
    }


def test_evaluate_external_dataset_excludes_test_labels_absent_from_training(tmp_path):
    train_path = tmp_path / "train.csv"
    test_path = tmp_path / "test.csv"
    pd.DataFrame(
        [
            _row("train-1", "AAAAAAAAAAAAAAAAAAAAAAAAAAAA", "I-E"),
            _row("train-2", "AAAAAAAATAAAAAAAAAAAAAAAAAAA", "I-E"),
            _row("train-3", "CCCCCCCCCCCCCCCCCCCCCCCCCCCC", "I-F"),
            _row("train-4", "CCCCCCCCGCCCCCCCCCCCCCCCCCCC", "I-F"),
        ]
    ).to_csv(train_path, index=False)
    pd.DataFrame(
        [
            _row("test-1", "AAAAAAAAGAAAAAAAAAAAAAAAAAAA", "I-E"),
            _row("test-2", "CCCCCCCCTCCCCCCCCCCCCCCCCCCC", "I-F"),
            _row("test-3", "TTTTTTTTTTTTTTTTTTTTTTTTTTTT", "II-A"),
        ]
    ).to_csv(test_path, index=False)

    result = evaluate_external_dataset(train_path, test_path, min_class_count=1, n_estimators=20)

    assert result["train_rows"] == 4
    assert result["raw_test_rows"] == 3
    assert result["evaluated_test_rows"] == 2
    assert result["excluded_test_rows"] == 1
    assert result["excluded_test_subtypes"] == ["II-A"]
    assert 0.0 <= result["accuracy"] <= 1.0
