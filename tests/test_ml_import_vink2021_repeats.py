from __future__ import annotations

import pandas as pd

from crispr_phage_predictor.ml.import_vink2021_repeats import (
    import_vink2021_candidate_repeats,
)


def test_import_vink2021_candidate_repeats_filters_to_matching_proximity(tmp_path):
    input_path = tmp_path / "vink.csv"
    repeat = "CCAGCCGCCTTCAGGCGGCTGTGTGTTGAAAC"
    pd.DataFrame(
        [
            {
                "spacers": "A" * 32,
                "repeats": repeat,
                "accessionnrs": "CP000001.1",
                "subtype": "CAS-TypeIC",
                "subtypesinproximity": "CAS-TypeIC",
                "genus": "Example",
                "family": "Exampleaceae",
                "order": "Exampleales",
                "class": "Exampleia",
                "phylum": "Exampleota",
                "superkingdom": "Bacteria",
                "PAM": "TTCX",
            },
            {
                "spacers": "C" * 32,
                "repeats": repeat,
                "accessionnrs": "CP000001.1",
                "subtype": "CAS-TypeIC",
                "subtypesinproximity": "CAS-TypeIC",
                "genus": "Example",
                "family": "Exampleaceae",
                "order": "Exampleales",
                "class": "Exampleia",
                "phylum": "Exampleota",
                "superkingdom": "Bacteria",
                "PAM": "TTCX",
            },
            {
                "spacers": "G" * 32,
                "repeats": "G" * 32,
                "accessionnrs": "CP000002.1",
                "subtype": "CAS-TypeIE",
                "subtypesinproximity": "CAS-TypeIF",
                "genus": "Mismatch",
                "family": "",
                "order": "",
                "class": "",
                "phylum": "",
                "superkingdom": "Bacteria",
                "PAM": "",
            },
        ]
    ).to_csv(input_path, index=False)

    table = import_vink2021_candidate_repeats(input_path)

    assert len(table) == 1
    row = table.iloc[0]
    assert row["genome_id"] == "CP000001.1"
    assert row["repeat_sequence"] == repeat
    assert row["spacer_count"] == 2
    assert row["cas_type"] == "Type I"
    assert row["cas_subtype"] == "I-C"
    assert row["label_confidence"] == "computational_proximity"


def test_import_vink2021_candidate_repeats_normalizes_multi_letter_type(tmp_path):
    input_path = tmp_path / "vink.csv"
    repeat = "GTTTTAGAGCTATGCTGTTTTGAATGGTCCCAAAAC"
    pd.DataFrame(
        [
            {
                "spacers": "A" * 30,
                "repeats": repeat,
                "accessionnrs": "CP000003.1",
                "subtype": "CAS-TypeIIA",
                "subtypesinproximity": "CAS-TypeIIA",
                "genus": "Streptococcus",
                "family": "",
                "order": "",
                "class": "",
                "phylum": "",
                "superkingdom": "Bacteria",
                "PAM": "NGG",
            },
            {
                "spacers": "C" * 30,
                "repeats": repeat,
                "accessionnrs": "CP000003.1",
                "subtype": "CAS-TypeIIA",
                "subtypesinproximity": "CAS-TypeIIA",
                "genus": "Streptococcus",
                "family": "",
                "order": "",
                "class": "",
                "phylum": "",
                "superkingdom": "Bacteria",
                "PAM": "NGG",
            },
        ]
    ).to_csv(input_path, index=False)

    table = import_vink2021_candidate_repeats(input_path)

    assert table.iloc[0]["cas_type"] == "Type II"
    assert table.iloc[0]["cas_subtype"] == "II-A"
