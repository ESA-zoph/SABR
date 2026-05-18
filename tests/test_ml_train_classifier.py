import unittest

from crispr_phage_predictor.ml.dataset import empty_repeat_cas_training_table
from crispr_phage_predictor.ml.train_classifier import evaluate_classifier


def _row(genome_id: str, repeat: str, subtype: str) -> dict[str, object]:
    return {
        "source": "test",
        "genome_id": genome_id,
        "organism": "",
        "taxonomy": "",
        "assembly_level": "",
        "contig_id": "contig_1",
        "array_start": "",
        "array_end": "",
        "repeat_sequence": repeat,
        "repeat_length": len(repeat),
        "spacer_count": 4,
        "mean_spacer_length": 32.0,
        "cas_type": "Type I",
        "cas_subtype": subtype,
        "label_source": "nearby_cas_operon",
        "label_confidence": "high",
        "pam_rule": "",
    }


class TrainClassifierTests(unittest.TestCase):
    def test_evaluates_classifier(self):
        table = empty_repeat_cas_training_table()
        rows = [
            _row("g1", "AAAAAAAAAAAAAAAAAAAAAAAAAAAA", "I-E"),
            _row("g2", "AAAAAAAATAAAAAAAAAAAAAAAAAAA", "I-E"),
            _row("g3", "AAAAAAAAGAAAAAAAAAAAAAAAAAAA", "I-E"),
            _row("g4", "CCCCCCCCCCCCCCCCCCCCCCCCCCCC", "I-F"),
            _row("g5", "CCCCCCCCGCCCCCCCCCCCCCCCCCCC", "I-F"),
            _row("g6", "CCCCCCCCTCCCCCCCCCCCCCCCCCCC", "I-F"),
        ]
        for index, row in enumerate(rows):
            table.loc[index] = row

        result = evaluate_classifier(table, test_size=0.34, random_state=1)

        self.assertEqual(result.train_size + result.test_size, 6)
        self.assertGreaterEqual(result.accuracy, 0.0)
        self.assertIn("I-E", result.labels)
        self.assertIn("I-F", result.labels)


if __name__ == "__main__":
    unittest.main()
