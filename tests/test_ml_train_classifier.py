import unittest

from crispr_phage_predictor.ml.dataset import empty_repeat_cas_training_table
from crispr_phage_predictor.ml.train_classifier import evaluate_classifier, evaluate_methods


def _row(
    genome_id: str,
    repeat: str,
    subtype: str,
    organism: str | None = None,
) -> dict[str, object]:
    organism = organism or ("Alpha example" if subtype == "I-E" else "Beta example")
    return {
        "source": "test",
        "genome_id": genome_id,
        "organism": organism,
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

        self.assertEqual(result.method, "random_forest")
        self.assertEqual(result.train_size + result.test_size, 6)
        self.assertGreaterEqual(result.accuracy, 0.0)
        self.assertIn("I-E", result.labels)
        self.assertIn("I-F", result.labels)

    def test_compares_nearest_repeat_and_random_forest(self):
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

        results = evaluate_methods(table, test_size=0.34, random_state=1)

        self.assertEqual(
            [result.method for result in results],
            [
                "nearest_repeat",
                "logistic_regression",
                "linear_svm",
                "gradient_boosting",
                "extra_trees",
                "hybrid_extra_trees",
                "random_forest",
            ],
        )

    def test_supports_group_holdout_by_genome(self):
        table = empty_repeat_cas_training_table()
        rows = [
            _row("g1", "AAAAAAAAAAAAAAAAAAAAAAAAAAAA", "I-E", organism="Alpha one"),
            _row("g2", "AAAAAAAATAAAAAAAAAAAAAAAAAAA", "I-E", organism="Beta one"),
            _row("g3", "AAAAAAAAGAAAAAAAAAAAAAAAAAAA", "I-E", organism="Gamma one"),
            _row("g4", "CCCCCCCCCCCCCCCCCCCCCCCCCCCC", "I-F", organism="Alpha two"),
            _row("g5", "CCCCCCCCGCCCCCCCCCCCCCCCCCCC", "I-F", organism="Beta two"),
            _row("g6", "CCCCCCCCTCCCCCCCCCCCCCCCCCCC", "I-F", organism="Gamma two"),
            _row("g7", "AAAAAAAACAAAAAAAAAAAAAAAAAAA", "I-E", organism="Delta one"),
            _row("g8", "CCCCCCCCACCCCCCCCCCCCCCCCCCC", "I-F", organism="Delta two"),
        ]
        for index, row in enumerate(rows):
            table.loc[index] = row

        results = evaluate_methods(
            table,
            test_size=0.25,
            random_state=2,
            split_strategy="group_holdout",
            group_column="genome_id",
        )

        self.assertTrue(results)
        self.assertTrue(all(result.split_strategy == "group_holdout" for result in results))
        self.assertGreaterEqual(results[0].test_size, 1)

    def test_supports_derived_genus_holdout(self):
        table = empty_repeat_cas_training_table()
        rows = [
            _row("g1", "AAAAAAAAAAAAAAAAAAAAAAAAAAAA", "I-E", organism="Alpha one"),
            _row("g2", "AAAAAAAATAAAAAAAAAAAAAAAAAAA", "I-E", organism="Beta one"),
            _row("g3", "AAAAAAAAGAAAAAAAAAAAAAAAAAAA", "I-E", organism="Gamma one"),
            _row("g4", "CCCCCCCCCCCCCCCCCCCCCCCCCCCC", "I-F", organism="Alpha two"),
            _row("g5", "CCCCCCCCGCCCCCCCCCCCCCCCCCCC", "I-F", organism="Beta two"),
            _row("g6", "CCCCCCCCTCCCCCCCCCCCCCCCCCCC", "I-F", organism="Gamma two"),
            _row("g7", "AAAAAAAACAAAAAAAAAAAAAAAAAAA", "I-E", organism="Delta one"),
            _row("g8", "CCCCCCCCACCCCCCCCCCCCCCCCCCC", "I-F", organism="Delta two"),
        ]
        for index, row in enumerate(rows):
            table.loc[index] = row

        results = evaluate_methods(
            table,
            test_size=0.25,
            random_state=2,
            split_strategy="group_holdout",
            group_column="genus",
            methods=("nearest_repeat", "random_forest"),
        )

        self.assertEqual([result.method for result in results], ["nearest_repeat", "random_forest"])

    def test_filters_rare_classes_before_evaluation(self):
        table = empty_repeat_cas_training_table()
        rows = [
            _row("g1", "AAAAAAAAAAAAAAAAAAAAAAAAAAAA", "I-E"),
            _row("g2", "AAAAAAAATAAAAAAAAAAAAAAAAAAA", "I-E"),
            _row("g3", "AAAAAAAAGAAAAAAAAAAAAAAAAAAA", "I-E"),
            _row("g4", "CCCCCCCCCCCCCCCCCCCCCCCCCCCC", "I-F"),
            _row("g5", "CCCCCCCCGCCCCCCCCCCCCCCCCCCC", "I-F"),
            _row("g6", "CCCCCCCCTCCCCCCCCCCCCCCCCCCC", "I-F"),
            _row("g7", "GGGGGGGGGGGGGGGGGGGGGGGGGGGG", "V-A"),
        ]
        for index, row in enumerate(rows):
            table.loc[index] = row

        results = evaluate_methods(
            table,
            test_size=0.34,
            random_state=1,
            min_class_count=2,
            methods=("nearest_repeat", "random_forest"),
        )

        self.assertNotIn("V-A", results[0].labels)

    def test_supports_hybrid_extra_trees(self):
        table = empty_repeat_cas_training_table()
        rows = [
            _row("g1", "AAAAAAAAAAAAAAAAAAAAAAAAAAAA", "I-E", organism="Alpha one"),
            _row("g2", "AAAAAAAATAAAAAAAAAAAAAAAAAAA", "I-E", organism="Beta one"),
            _row("g3", "AAAAAAAAGAAAAAAAAAAAAAAAAAAA", "I-E", organism="Gamma one"),
            _row("g4", "CCCCCCCCCCCCCCCCCCCCCCCCCCCC", "I-F", organism="Alpha two"),
            _row("g5", "CCCCCCCCGCCCCCCCCCCCCCCCCCCC", "I-F", organism="Beta two"),
            _row("g6", "CCCCCCCCTCCCCCCCCCCCCCCCCCCC", "I-F", organism="Gamma two"),
            _row("g7", "AAAAAAAACAAAAAAAAAAAAAAAAAAA", "I-E", organism="Delta one"),
            _row("g8", "CCCCCCCCACCCCCCCCCCCCCCCCCCC", "I-F", organism="Delta two"),
        ]
        for index, row in enumerate(rows):
            table.loc[index] = row

        results = evaluate_methods(
            table,
            test_size=0.25,
            random_state=2,
            split_strategy="group_holdout",
            group_column="genus",
            methods=("hybrid_extra_trees",),
        )

        self.assertEqual([result.method for result in results], ["hybrid_extra_trees"])

    def test_supports_hierarchical_extra_trees(self):
        table = empty_repeat_cas_training_table()
        rows = [
            _row("g1", "AAAAAAAAAAAAAAAAAAAAAAAAAAAA", "I-E", organism="Alpha one"),
            _row("g2", "AAAAAAAATAAAAAAAAAAAAAAAAAAA", "I-E", organism="Beta one"),
            _row("g3", "AAAAAAAAGAAAAAAAAAAAAAAAAAAA", "I-E", organism="Gamma one"),
            _row("g4", "CCCCCCCCCCCCCCCCCCCCCCCCCCCC", "I-F", organism="Alpha two"),
            _row("g5", "CCCCCCCCGCCCCCCCCCCCCCCCCCCC", "I-F", organism="Beta two"),
            _row("g6", "CCCCCCCCTCCCCCCCCCCCCCCCCCCC", "I-F", organism="Gamma two"),
            _row("g7", "TTTTTTTTTTTTTTTTTTTTTTTTTTTT", "II-A", organism="Delta one"),
            _row("g8", "TTTTTTTTGTTTTTTTTTTTTTTTTTTT", "II-A", organism="Epsilon one"),
        ]
        for index, row in enumerate(rows):
            if row["cas_subtype"].startswith("II"):
                row["cas_type"] = "Type II"
            table.loc[index] = row

        results = evaluate_methods(
            table,
            test_size=0.25,
            random_state=2,
            split_strategy="group_holdout",
            group_column="genome_id",
            methods=("hierarchical_extra_trees",),
        )

        self.assertEqual([result.method for result in results], ["hierarchical_extra_trees"])
        self.assertGreaterEqual(results[0].accuracy, 0.0)


if __name__ == "__main__":
    unittest.main()
