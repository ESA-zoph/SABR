import unittest

from crispr_phage_predictor.ml.dataset import empty_repeat_cas_training_table
from crispr_phage_predictor.ml.features import build_repeat_feature_table, feature_columns


class RepeatFeatureTests(unittest.TestCase):
    def test_builds_numeric_repeat_features(self):
        training = empty_repeat_cas_training_table()
        training.loc[0] = {
            "source": "cctyper",
            "genome_id": "GCF_000001",
            "organism": "",
            "taxonomy": "",
            "assembly_level": "",
            "contig_id": "contig_1",
            "array_start": "",
            "array_end": "",
            "repeat_sequence": "ACGTACGT",
            "repeat_length": 8,
            "spacer_count": 3,
            "mean_spacer_length": 32.5,
            "cas_type": "Type I",
            "cas_subtype": "I-E",
            "label_source": "nearby_cas_operon",
            "label_confidence": "high",
            "pam_rule": "",
        }

        features = build_repeat_feature_table(training, kmer_sizes=(2,))

        self.assertEqual(features.loc[0, "repeat_length"], 8)
        self.assertEqual(features.loc[0, "repeat_gc_percent"], 50.0)
        self.assertEqual(features.loc[0, "spacer_count"], 3)
        self.assertEqual(features.loc[0, "mean_spacer_length"], 32.5)
        self.assertAlmostEqual(features.loc[0, "kmer_2_AC"], 2 / 7)
        self.assertAlmostEqual(features.loc[0, "kmer_2_GT"], 2 / 7)
        self.assertEqual(features.loc[0, "cas_subtype"], "I-E")

    def test_feature_columns_exclude_labels_and_ids(self):
        training = empty_repeat_cas_training_table()
        training.loc[0] = {
            "source": "cctyper",
            "genome_id": "GCF_000001",
            "organism": "",
            "taxonomy": "",
            "assembly_level": "",
            "contig_id": "contig_1",
            "array_start": "",
            "array_end": "",
            "repeat_sequence": "AAAA",
            "repeat_length": 4,
            "spacer_count": 1,
            "mean_spacer_length": 30.0,
            "cas_type": "Type I",
            "cas_subtype": "I-E",
            "label_source": "nearby_cas_operon",
            "label_confidence": "high",
            "pam_rule": "",
        }
        features = build_repeat_feature_table(training, kmer_sizes=(2,))

        columns = feature_columns(features)

        self.assertIn("repeat_length", columns)
        self.assertIn("kmer_2_AA", columns)
        self.assertNotIn("genome_id", columns)
        self.assertNotIn("cas_subtype", columns)


if __name__ == "__main__":
    unittest.main()
