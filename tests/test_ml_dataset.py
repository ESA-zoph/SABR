import unittest

import pandas as pd

from crispr_phage_predictor.ml.dataset import (
    build_training_table_from_cctyper,
    cas_type_from_subtype,
    empty_repeat_cas_training_table,
    filter_high_confidence_labels,
    validate_repeat_cas_training_table,
)


class RepeatCasDatasetTests(unittest.TestCase):
    def test_validates_expected_training_schema(self):
        table = empty_repeat_cas_training_table()
        table.loc[0] = {
            "source": "cctyper",
            "genome_id": "GCF_000001",
            "organism": "Example bacterium",
            "taxonomy": "Bacteria",
            "assembly_level": "complete genome",
            "contig_id": "contig_1",
            "array_start": 100,
            "array_end": 400,
            "repeat_sequence": "GTTCACTGCCGTACAGGCAGCTTAGAAA",
            "repeat_length": 28,
            "spacer_count": 4,
            "mean_spacer_length": 32.0,
            "cas_type": "Type I",
            "cas_subtype": "I-E",
            "label_source": "nearby_cas_operon",
            "label_confidence": "high",
            "pam_rule": "5prime:AWG",
        }

        validate_repeat_cas_training_table(table)

    def test_rejects_repeat_length_mismatch(self):
        table = empty_repeat_cas_training_table()
        table.loc[0] = {
            "source": "cctyper",
            "genome_id": "GCF_000001",
            "organism": "",
            "taxonomy": "",
            "assembly_level": "",
            "contig_id": "contig_1",
            "array_start": "",
            "array_end": "",
            "repeat_sequence": "ACGT",
            "repeat_length": 5,
            "spacer_count": 1,
            "mean_spacer_length": 32.0,
            "cas_type": "Type I",
            "cas_subtype": "I-E",
            "label_source": "nearby_cas_operon",
            "label_confidence": "high",
            "pam_rule": "",
        }

        with self.assertRaises(ValueError):
            validate_repeat_cas_training_table(table)

    def test_filters_high_confidence_labels(self):
        table = empty_repeat_cas_training_table()
        base_row = {
            "source": "cctyper",
            "genome_id": "GCF_000001",
            "organism": "",
            "taxonomy": "",
            "assembly_level": "",
            "contig_id": "contig_1",
            "array_start": "",
            "array_end": "",
            "repeat_sequence": "ACGT",
            "repeat_length": 4,
            "spacer_count": 1,
            "mean_spacer_length": 32.0,
            "cas_type": "Type I",
            "cas_subtype": "I-E",
            "label_source": "nearby_cas_operon",
            "label_confidence": "high",
            "pam_rule": "",
        }
        table.loc[0] = base_row
        table.loc[1] = {**base_row, "genome_id": "GCF_000002", "label_confidence": "low"}

        filtered = filter_high_confidence_labels(table)

        self.assertEqual(list(filtered["genome_id"]), ["GCF_000001"])

    def test_builds_training_table_from_cctyper_crisprs_near_cas(self):
        cctyper = pd.DataFrame(
            [
                {
                    "Contig": "NC_000001",
                    "CRISPR": "NC_000001_1",
                    "Start": 100,
                    "End": 500,
                    "Consensus_repeat": "GTTCACTGCCGTACAGGCAGCTTAGAAA",
                    "N_repeats": 5,
                    "Repeat_len": 28,
                    "Spacer_len_avg": 32.5,
                    "Trusted": "TRUE",
                    "Subtype": "I-E",
                    "Subtype_probability": 0.95,
                }
            ]
        )

        table = build_training_table_from_cctyper(
            cctyper,
            genome_id="GCF_000001",
            organism="Example bacterium",
            assembly_level="complete genome",
            pam_rules={"I-E": "5prime:AWG"},
        )

        self.assertEqual(len(table), 1)
        self.assertEqual(table.loc[0, "cas_type"], "Type I")
        self.assertEqual(table.loc[0, "cas_subtype"], "I-E")
        self.assertEqual(table.loc[0, "spacer_count"], 4)
        self.assertEqual(table.loc[0, "label_confidence"], "high")
        self.assertEqual(table.loc[0, "pam_rule"], "5prime:AWG")

    def test_skips_ambiguous_cctyper_subtypes(self):
        cctyper = pd.DataFrame(
            [
                {
                    "Contig": "NC_000001",
                    "Start": 100,
                    "End": 500,
                    "Consensus_repeat": "ACGTACGTACGTACGTACGTACGT",
                    "N_repeats": 3,
                    "Subtype": "I-E/I-F",
                }
            ]
        )

        table = build_training_table_from_cctyper(cctyper, genome_id="GCF_000001")

        self.assertTrue(table.empty)

    def test_maps_subtype_to_broad_type(self):
        self.assertEqual(cas_type_from_subtype("II-A"), "Type II")
        self.assertEqual(cas_type_from_subtype("V-A"), "Type V")
        self.assertEqual(cas_type_from_subtype("VI-B1"), "Type VI")


if __name__ == "__main__":
    unittest.main()
