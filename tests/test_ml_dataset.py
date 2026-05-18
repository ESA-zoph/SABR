import unittest

from crispr_phage_predictor.ml.dataset import (
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


if __name__ == "__main__":
    unittest.main()
