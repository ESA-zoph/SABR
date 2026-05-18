import unittest

from crispr_phage_predictor.ml.classifier import NearestRepeatClassifier, RepeatCasSubtypeClassifier
from crispr_phage_predictor.ml.dataset import empty_repeat_cas_training_table


def _training_row(genome_id: str, repeat: str, subtype: str) -> dict[str, object]:
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


class RepeatCasSubtypeClassifierTests(unittest.TestCase):
    def test_trains_and_predicts_subtype(self):
        training = empty_repeat_cas_training_table()
        training.loc[0] = _training_row("g1", "AAAAAAAAAAAAAAAAAAAAAAAAAAAA", "I-E")
        training.loc[1] = _training_row("g2", "AAAAAAAATAAAAAAAAAAAAAAAAAAA", "I-E")
        training.loc[2] = _training_row("g3", "CCCCCCCCCCCCCCCCCCCCCCCCCCCC", "I-F")
        training.loc[3] = _training_row("g4", "CCCCCCCCGCCCCCCCCCCCCCCCCCCC", "I-F")

        classifier = RepeatCasSubtypeClassifier(kmer_sizes=(2,), n_estimators=50)
        classifier.fit(training)
        prediction = classifier.predict_one("AAAAAAAAAAAAAAAAAAAAAAAAAAAA", spacer_count=4)

        self.assertEqual(prediction.cas_subtype, "I-E")
        self.assertGreaterEqual(prediction.confidence, 0.5)
        self.assertIn("I-E", prediction.probabilities)
        self.assertIn("I-F", prediction.probabilities)

    def test_requires_two_subtypes_for_training(self):
        training = empty_repeat_cas_training_table()
        training.loc[0] = _training_row("g1", "AAAAAAAAAAAAAAAAAAAAAAAAAAAA", "I-E")
        training.loc[1] = _training_row("g2", "AAAAAAAATAAAAAAAAAAAAAAAAAAA", "I-E")

        classifier = RepeatCasSubtypeClassifier(kmer_sizes=(2,), n_estimators=10)

        with self.assertRaises(ValueError):
            classifier.fit(training)

    def test_requires_fit_before_prediction(self):
        classifier = RepeatCasSubtypeClassifier(kmer_sizes=(2,), n_estimators=10)

        with self.assertRaises(ValueError):
            classifier.predict_one("AAAAAAAAAAAAAAAAAAAAAAAAAAAA")


class NearestRepeatClassifierTests(unittest.TestCase):
    def test_predicts_subtype_from_nearest_repeat(self):
        training = empty_repeat_cas_training_table()
        training.loc[0] = _training_row("g1", "AAAAAAAAAAAAAAAAAAAAAAAAAAAA", "I-E")
        training.loc[1] = _training_row("g2", "CCCCCCCCCCCCCCCCCCCCCCCCCCCC", "I-F")

        classifier = NearestRepeatClassifier()
        classifier.fit(training)
        prediction = classifier.predict_one("AAAAAAAATAAAAAAAAAAAAAAAAAAA")

        self.assertEqual(prediction.cas_subtype, "I-E")
        self.assertGreater(prediction.best_identity, 0.9)
        self.assertEqual(prediction.confidence, prediction.best_identity)

    def test_requires_fit_before_similarity_prediction(self):
        classifier = NearestRepeatClassifier()

        with self.assertRaises(ValueError):
            classifier.predict_one("AAAAAAAAAAAAAAAAAAAAAAAAAAAA")


if __name__ == "__main__":
    unittest.main()
