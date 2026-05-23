from __future__ import annotations

import pandas as pd

from crispr_phage_predictor.ml.audit_crisprcasdb_candidates import (
    audit_crisprcasdb_candidates,
)


def _row(genome_id: str, repeat: str, subtype: str, source: str = "test") -> dict[str, object]:
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
        "label_confidence": "high",
        "pam_rule": "",
    }


def test_audit_crisprcasdb_candidates_reports_overlap_conflicts_and_balanced_subset(tmp_path):
    current_path = tmp_path / "current.csv"
    candidate_path = tmp_path / "candidate.csv"
    shared_repeat = "GTTTCAATGCTGCTTCGCCTGCAATGGGTTTAGTAT"
    novel_repeat = "CCAGCCGCCTTCAGGCGGCTGTGTGTTGAAAC"
    conflict_repeat = "GTTTTAGAGCTATGCTGTTTTGAATGGTCCCAAAAC"
    pd.DataFrame(
        [
            _row("current-1", shared_repeat, "I-F"),
            _row("current-2", "ATTTTAAAGCTATGCTGTTTTGAATGGTCCCAAAAC", "III-A"),
        ]
    ).to_csv(current_path, index=False)
    pd.DataFrame(
        [
            _row("candidate-1", shared_repeat, "I-F", source="candidate"),
            _row("candidate-2", novel_repeat, "I-E", source="candidate"),
            _row("candidate-3", conflict_repeat, "I-E", source="candidate"),
            _row("candidate-4", conflict_repeat, "III-A", source="candidate"),
        ]
    ).to_csv(candidate_path, index=False)

    outputs = audit_crisprcasdb_candidates(current_path, candidate_path, max_per_subtype=1)
    summary = dict(zip(outputs["summary"]["metric"], outputs["summary"]["value"], strict=True))

    assert summary["current_rows"] == 2
    assert summary["candidate_rows"] == 4
    assert summary["overlapping_repeat_hashes"] == 1
    assert summary["candidate_conflicting_repeat_hashes"] == 1
    assert summary["novel_nonconflicting_candidate_rows"] == 1
    assert summary["balanced_candidate_rows"] == 1
    conflicts = outputs["candidate_repeat_conflicts"]
    assert len(conflicts) == 1
    assert conflicts.iloc[0]["subtypes"] == "I-E;III-A"
