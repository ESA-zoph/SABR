import tempfile
import unittest
from pathlib import Path

from crispr_phage_predictor.cas_prediction import predict_array_cas_subtypes
from crispr_phage_predictor.crispr import CrisprArray
from crispr_phage_predictor.ml.dataset import empty_repeat_cas_training_table
from crispr_phage_predictor.ml.model_artifact import (
    save_artifact,
    train_extra_trees_artifact,
)


def _training_row(repeat: str, subtype: str) -> dict[str, object]:
    return {
        "source": "test",
        "genome_id": subtype,
        "organism": "",
        "taxonomy": "",
        "assembly_level": "",
        "contig_id": "contig_1",
        "array_start": "",
        "array_end": "",
        "repeat_sequence": repeat,
        "repeat_length": len(repeat),
        "spacer_count": 3,
        "mean_spacer_length": 32.0,
        "cas_type": "Type I",
        "cas_subtype": subtype,
        "label_source": "test",
        "label_confidence": "high",
        "pam_rule": "",
    }


class CasPredictionTests(unittest.TestCase):
    def test_predicts_array_subtype_and_curated_pam_rule(self):
        table = empty_repeat_cas_training_table()
        table.loc[0] = _training_row("GTTCACTGCCGTACAGGCAGCTTAGAAA", "I-F")
        table.loc[1] = _training_row("AAAAAAAAAAAAAAAAAAAAAAAAAAAA", "I-E")
        array = CrisprArray(
            array_id="array_1",
            genome_id="bacterium.fna",
            contig_id="contig_1",
            start=1,
            end=100,
            repeat_consensus="GTTCACTGCCGTACAGGCAGCTTAGAAA",
            spacers=["ACGT" * 8],
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            training_path = Path(temp_dir) / "training.csv"
            table.to_csv(training_path, index=False)
            predictions = predict_array_cas_subtypes([array], training_table_path=training_path)

        prediction = predictions["array_1"]
        self.assertEqual(prediction.cas_subtype, "I-F")
        self.assertEqual(prediction.pam_rule, "genomic_3prime:GG")
        self.assertEqual(prediction.pam_rule_source, "curated_subtype_catalog")

    def test_prefers_extra_trees_artifact_when_available(self):
        table = empty_repeat_cas_training_table()
        table.loc[0] = _training_row("GTTCACTGCCGTACAGGCAGCTTAGAAA", "I-F")
        table.loc[1] = _training_row("GTTCACTGCCGTACAGGCAGCTTAGAAT", "I-F")
        table.loc[2] = _training_row("AAAAAAAAAAAAAAAAAAAAAAAAAAAA", "I-E")
        table.loc[3] = _training_row("AAAAAAAATAAAAAAAAAAAAAAAAAAA", "I-E")
        array = CrisprArray(
            array_id="array_1",
            genome_id="bacterium.fna",
            contig_id="contig_1",
            start=1,
            end=100,
            repeat_consensus="GTTCACTGCCGTACAGGCAGCTTAGAAA",
            spacers=["ACGT" * 8],
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            training_path = Path(temp_dir) / "training.csv"
            artifact_path = Path(temp_dir) / "model.joblib"
            table.to_csv(training_path, index=False)
            artifact = train_extra_trees_artifact(
                training_table_path=training_path,
                min_class_count=1,
                n_estimators=20,
            )
            save_artifact(artifact, artifact_path)
            predictions = predict_array_cas_subtypes(
                [array],
                training_table_path=training_path,
                model_artifact_path=artifact_path,
                min_confidence_for_pam=0.0,
            )

        prediction = predictions["array_1"]
        self.assertEqual(prediction.prediction_method, "extra_trees")
        self.assertEqual(prediction.cas_subtype, "I-F")


if __name__ == "__main__":
    unittest.main()
