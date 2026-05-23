from __future__ import annotations

import pandas as pd

from crispr_phage_predictor.ml.collect_curated_minced_training import (
    collect_curated_training_table,
)


def test_collect_curated_training_table_with_internal_detector(tmp_path):
    repeat = "ATGCGTACGTTAGCTAGCTAGGCTA"
    spacer_a = "AACCGGTTAACCGGTTAACC"
    spacer_b = "TTGGAACCTTGGAACCTTGG"
    fasta_path = tmp_path / "genome.fasta"
    fasta_path.write_text(
        f">contig_1\n{repeat}{spacer_a}{repeat}{spacer_b}{repeat}\n",
        encoding="utf-8",
    )
    manifest_path = tmp_path / "manifest.tsv"
    pd.DataFrame(
        [
            {
                "fasta_path": str(fasta_path),
                "genome_id": "GCF_TEST",
                "organism": "Example bacterium",
                "taxonomy": "Bacteria",
                "assembly_level": "complete genome",
                "cas_type": "Type I",
                "cas_subtype": "I-F",
                "label_source": "curated_literature",
                "label_confidence": "curated",
                "label_scope": "genome",
                "pam_rule": "CC",
            }
        ]
    ).to_csv(manifest_path, sep="\t", index=False)

    table = collect_curated_training_table(manifest_path, detector="internal")

    assert len(table) == 1
    row = table.iloc[0]
    assert row["genome_id"] == "GCF_TEST"
    assert row["contig_id"] == "contig_1"
    assert row["repeat_sequence"] == repeat
    assert row["spacer_count"] == 2
    assert row["cas_subtype"] == "I-F"
    assert row["label_source"] == "curated_literature"
