from __future__ import annotations

import pandas as pd

from crispr_phage_predictor.ml.build_crisprcasdb_augmented_dataset import (
    build_crisprcasdb_augmented_dataset,
)


def _row(genome_id: str, repeat: str, subtype: str, source: str = "current") -> dict[str, object]:
    return {
        "source": source,
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
        "cas_type": "Type I" if subtype.startswith("I-") else "Type III",
        "cas_subtype": subtype,
        "label_source": source,
        "label_confidence": "high" if source == "current" else "computational_nearby_cas_cluster",
        "pam_rule": "",
    }


def test_build_crisprcasdb_augmented_dataset_keeps_only_novel_nonconflicting_candidates(tmp_path):
    current_path = tmp_path / "current.csv"
    candidate_path = tmp_path / "candidate.csv"
    shared_repeat = "GTTTCAATGCTGCTTCGCCTGCAATGGGTTTAGTAT"
    novel_repeat = "CCAGCCGCCTTCAGGCGGCTGTGTGTTGAAAC"
    conflict_repeat = "GTTTTAGAGCTATGCTGTTTTGAATGGTCCCAAAAC"
    pd.DataFrame([_row("current-1", shared_repeat, "I-F")]).to_csv(current_path, index=False)
    pd.DataFrame(
        [
            _row("candidate-1", shared_repeat, "I-F", source="candidate"),
            _row("candidate-2", novel_repeat, "I-E", source="candidate"),
            _row("candidate-3", conflict_repeat, "I-E", source="candidate"),
            _row("candidate-4", conflict_repeat, "III-A", source="candidate"),
        ]
    ).to_csv(candidate_path, index=False)

    augmented, additions = build_crisprcasdb_augmented_dataset(current_path, candidate_path)

    assert len(additions) == 1
    assert additions.iloc[0]["repeat_sequence"] == novel_repeat
    assert additions.iloc[0]["cas_subtype"] == "I-E"
    assert len(augmented) == 2
    assert list(augmented["cas_subtype"]) == ["I-F", "I-E"]


def test_build_crisprcasdb_augmented_dataset_caps_per_subtype(tmp_path):
    current_path = tmp_path / "current.csv"
    candidate_path = tmp_path / "candidate.csv"
    pd.DataFrame(
        [_row("current-1", "GTTTCAATGCTGCTTCGCCTGCAATGGGTTTAGTAT", "I-F")]
    ).to_csv(current_path, index=False)
    candidates = [
        _row("candidate-1", "CCAGCCGCCTTCAGGCGGCTGTGTGTTGAAAC", "I-E", source="candidate"),
        _row("candidate-2", "CCAGCCGCCTTCAGGCGGCTGTGTGTTGAAAT", "I-E", source="candidate"),
        _row("candidate-3", "CCAGCCGCCTTCAGGCGGCTGTGTGTTGAAG", "I-E", source="candidate"),
    ]
    pd.DataFrame(candidates).to_csv(candidate_path, index=False)

    _, additions = build_crisprcasdb_augmented_dataset(
        current_path,
        candidate_path,
        max_per_subtype=2,
    )

    assert len(additions) == 2
