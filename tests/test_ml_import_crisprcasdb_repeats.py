from __future__ import annotations

import zipfile

from crispr_phage_predictor.ml.import_crisprcasdb_repeats import (
    import_crisprcasdb_direct_repeats,
)


def test_import_crisprcasdb_direct_repeats_builds_unlabeled_inventory(tmp_path):
    zip_path = tmp_path / "dr_34.zip"
    with zipfile.ZipFile(zip_path, "w") as archive:
        archive.writestr(
            "direct_repeat_seqName.fsa",
            "\n".join(
                [
                    ">CP000001.1+CP000002.1",
                    "GTTTCAATGCTGCTTCGCCTGCAATGGGTTTAGTAT",
                    ">uuid-with-ambiguous-base",
                    "YGGTTTATCCCCGCTGGCGCGGGGAA",
                ]
            ),
        )

    table = import_crisprcasdb_direct_repeats(zip_path)

    assert list(table["record_id"]) == ["CP000001.1+CP000002.1", "uuid-with-ambiguous-base"]
    assert table.iloc[0]["source"] == "crisprcasdb_direct_repeat_fasta"
    assert table.iloc[0]["release"] == "34"
    assert table.iloc[0]["first_accession"] == "CP000001.1"
    assert table.iloc[0]["accession_count"] == 2
    assert table.iloc[0]["valid_iupac_dna"]
    assert table.iloc[0]["usable_for_sabr_repeat_features"]
    assert table.iloc[1]["valid_iupac_dna"]
    assert not table.iloc[1]["usable_for_sabr_repeat_features"]


def test_import_crisprcasdb_direct_repeats_can_keep_only_feature_usable_rows(tmp_path):
    zip_path = tmp_path / "dr_34.zip"
    with zipfile.ZipFile(zip_path, "w") as archive:
        archive.writestr(
            "direct_repeat_id.fsa",
            "\n".join(
                [
                    ">usable",
                    "GTTTCAATGCTGCTTCGCCTGCAATGGGTTTAGTAT",
                    ">too_short",
                    "ACGT",
                    ">ambiguous",
                    "YGGTTTATCCCCGCTGGCGCGGGGAA",
                ]
            ),
        )

    table = import_crisprcasdb_direct_repeats(
        zip_path,
        member="direct_repeat_id.fsa",
        only_usable=True,
    )

    assert list(table["record_id"]) == ["usable"]
