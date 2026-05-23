import tempfile
import unittest
import json
from pathlib import Path

import joblib
import pandas as pd

from crispr_phage_predictor.crispr import CrisprArray
from crispr_phage_predictor.io import FastaRecord
from crispr_phage_predictor.matching import SpacerHit
from crispr_phage_predictor.ml.model_artifact import model_artifact_metadata
from crispr_phage_predictor.output import save_analysis_run


class OutputTests(unittest.TestCase):
    def test_saves_analysis_run_tables(self):
        bacteria = [FastaRecord("b.fna", "b1", "b1", "ATGC")]
        phages = [FastaRecord("p.fna", "p1", "p1", "ATGC")]
        arrays = [
            CrisprArray(
                array_id="b.fna|b1|array_1",
                genome_id="b.fna",
                contig_id="b1",
                start=1,
                end=100,
                repeat_consensus="GTTCACTGCCGTACAGGCAGCTTAGAAA",
                spacers=["ACGTACGTACGTACGTACGTACGTACGTACGT"],
            )
        ]
        hits = [
            SpacerHit(
                bacterium_id="b.fna",
                array_id="b.fna|b1|array_1",
                phage_id="p.fna",
                spacer_id="b.fna|b1|array_1|spacer_1",
                phage_contig_id="p1",
                start=1,
                end=32,
                strand="+",
                identity=1.0,
                mismatches=0,
                alignment_length=32,
                spacer_length=32,
                coverage=1.0,
                evalue=None,
                bitscore=None,
                spacer_sequence="ACGTACGTACGTACGTACGTACGTACGTACGT",
                protospacer_sequence="ACGTACGTACGTACGTACGTACGTACGTACGT",
            )
        ]
        evidence = pd.DataFrame(
            [
                {
                    "bacterium": "b.fna",
                    "phage": "p.fna",
                    "spacer_hits": 1,
                    "unique_matching_spacers": 1,
                    "best_identity_percent": 100.0,
                }
            ]
        )
        heatmap = pd.DataFrame({"p.fna": [1]}, index=["b.fna"])

        with tempfile.TemporaryDirectory() as temp_dir:
            run_dir = save_analysis_run(
                bacteria_records=bacteria,
                phage_records=phages,
                crispr_arrays=arrays,
                spacer_hits=hits,
                evidence_matrix=evidence,
                heatmap=heatmap,
                detection_method="internal",
                matching_method="blast",
                output_root=Path(temp_dir),
            )

            self.assertTrue((run_dir / "run_metadata.json").exists())
            self.assertTrue((run_dir / "evidence_matrix.csv").exists())
            self.assertTrue((run_dir / "spacer_hits.csv").exists())
            metadata = json.loads((run_dir / "run_metadata.json").read_text())
            self.assertIn("cas_model_artifact", metadata)
            self.assertIn("artifact_sha256", metadata["cas_model_artifact"])

    def test_reads_model_artifact_metadata(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            artifact_path = Path(temp_dir) / "model.joblib"
            joblib.dump(
                {
                    "model": object(),
                    "feature_names": ["repeat_length"],
                    "classes": ["I-F", "II-A"],
                    "metadata": {
                        "method": "extra_trees",
                        "training_table": "training.csv",
                        "training_rows": 123,
                        "min_class_count": 20,
                        "random_state": 42,
                        "n_estimators": 400,
                    },
                },
                artifact_path,
            )

            metadata = model_artifact_metadata(artifact_path)

            self.assertTrue(metadata["artifact_exists"])
            self.assertEqual(len(metadata["artifact_sha256"]), 64)
            self.assertEqual(metadata["method"], "extra_trees")
            self.assertEqual(metadata["training_rows"], 123)
            self.assertEqual(metadata["classes"], ["I-F", "II-A"])
            self.assertEqual(metadata["load_error"], "")

    def test_missing_model_artifact_metadata_is_graceful(self):
        metadata = model_artifact_metadata(Path("does-not-exist.joblib"))

        self.assertFalse(metadata["artifact_exists"])
        self.assertEqual(metadata["artifact_sha256"], "")
        self.assertEqual(metadata["classes"], [])


if __name__ == "__main__":
    unittest.main()
