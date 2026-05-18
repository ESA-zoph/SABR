import tempfile
import unittest
from pathlib import Path

import pandas as pd

from crispr_phage_predictor.ml.collect_cctyper_training import collect_cctyper_training_table


class CollectCCTyperTrainingTests(unittest.TestCase):
    def test_collects_cctyper_output_directories(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            output_dir = root / "GCF_000001_cctyper"
            output_dir.mkdir()
            pd.DataFrame(
                [
                    {
                        "Contig": "NC_000001",
                        "Start": 10,
                        "End": 200,
                        "Consensus_repeat": "GTTCACTGCCGTACAGGCAGCTTAGAAA",
                        "N_repeats": 4,
                        "Repeat_len": 28,
                        "Spacer_len_avg": 32,
                        "Trusted": "TRUE",
                        "Subtype": "I-E",
                        "Subtype_probability": 0.95,
                    }
                ]
            ).to_csv(output_dir / "crisprs_near_cas.tab", sep="\t", index=False)

            manifest_path = root / "manifest.csv"
            pd.DataFrame(
                [
                    {
                        "cctyper_output_dir": str(output_dir),
                        "genome_id": "GCF_000001",
                        "organism": "Example bacterium",
                        "assembly_level": "complete genome",
                    }
                ]
            ).to_csv(manifest_path, index=False)

            output_path = root / "repeats_cas_types.csv"
            table = collect_cctyper_training_table(manifest_path, output_path)

            self.assertTrue(output_path.exists())
            self.assertEqual(len(table), 1)
            self.assertEqual(table.loc[0, "genome_id"], "GCF_000001")
            self.assertEqual(table.loc[0, "cas_subtype"], "I-E")


if __name__ == "__main__":
    unittest.main()
