import tempfile
import unittest
from pathlib import Path

import pandas as pd

from crispr_phage_predictor.ml.run_cctyper_batch import discover_cctyper_jobs, write_manifest


class RunCCTyperBatchTests(unittest.TestCase):
    def test_discovers_resistant_and_susceptible_fastas(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            (root / "resistant").mkdir()
            (root / "susceptible").mkdir()
            (root / "resistant" / "PA14 genome.fasta").write_text(">seq\nACGT\n")
            (root / "susceptible" / "DGCC7710.fna").write_text(">seq\nACGT\n")
            (root / "susceptible" / "notes.txt").write_text("ignore me\n")

            jobs = discover_cctyper_jobs(root, root / "runs")

            self.assertEqual(len(jobs), 2)
            self.assertEqual([job.phenotype_label for job in jobs], ["resistant", "susceptible"])
            self.assertEqual(jobs[0].genome_id, "PA14_genome")
            self.assertEqual(jobs[1].genome_id, "DGCC7710")
            self.assertEqual(jobs[0].output_dir, root / "runs" / "resistant" / "PA14_genome_cctyper")

    def test_writes_collect_cctyper_manifest_columns(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            (root / "resistant").mkdir()
            (root / "resistant" / "example.fasta").write_text(">seq\nACGT\n")

            jobs = discover_cctyper_jobs(root, root / "runs")
            manifest_path = root / "manifest.csv"
            write_manifest(jobs, manifest_path)

            manifest = pd.read_csv(manifest_path)

            self.assertIn("cctyper_output_dir", manifest.columns)
            self.assertIn("genome_id", manifest.columns)
            self.assertEqual(manifest.loc[0, "source_group"], "resistant")
            self.assertEqual(manifest.loc[0, "phenotype_label"], "resistant")


if __name__ == "__main__":
    unittest.main()
