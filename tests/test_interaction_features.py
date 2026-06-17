import tempfile
import unittest
from pathlib import Path

import pandas as pd

from crispr_phage_predictor.accession_linkage import empty_accession_linkage_table
from crispr_phage_predictor.interaction_features import build_hybrid_interaction_feature_table


class InteractionFeatureTests(unittest.TestCase):
    def test_builds_features_for_hybrid_ready_pair(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            host_path = root / "host.fasta"
            phage_path = root / "phage.fasta"
            host_path.write_text(">host\nACGTACGT\n", encoding="utf-8")
            phage_path.write_text(">phage\nAAAACCCC\n", encoding="utf-8")

            interactions = pd.DataFrame(
                [
                    {
                        "interaction_id": "pair_1",
                        "source_key": "source",
                        "bacterium": "Bacterium",
                        "strain": "strain_1",
                        "phage": "Phage",
                        "eop_class": "high",
                        "susceptibility_label": "susceptible",
                        "eop_value": "1.0",
                    }
                ]
            )
            linkage = empty_accession_linkage_table()
            linkage.loc[0] = {
                "linkage_id": "host",
                "entity_type": "bacterium",
                "source_key": "source",
                "display_name": "Bacterium",
                "strain_or_isolate": "*",
                "accession": "NC_1",
                "accession_database": "RefSeq",
                "assembly_level": "complete_genome",
                "sequence_status": "local",
                "linkage_status": "reference_proxy",
                "confidence": "medium",
                "local_path": str(host_path),
                "notes": "",
            }
            linkage.loc[1] = {
                "linkage_id": "phage",
                "entity_type": "phage",
                "source_key": "source",
                "display_name": "Phage",
                "strain_or_isolate": "Phage",
                "accession": "NC_2",
                "accession_database": "GenBank",
                "assembly_level": "complete_genome",
                "sequence_status": "local",
                "linkage_status": "exact",
                "confidence": "high",
                "local_path": str(phage_path),
                "notes": "",
            }
            coverage = pd.DataFrame(
                [
                    {
                        "interaction_id": "pair_1",
                        "pair_hybrid_ready": True,
                    }
                ]
            )

            features = build_hybrid_interaction_feature_table(
                interactions,
                linkage,
                coverage,
                k=2,
            )

        self.assertEqual(len(features), 1)
        self.assertEqual(features.loc[0, "binary_susceptibility"], "susceptible")
        self.assertTrue(bool(features.loc[0, "uses_reference_proxy_host"]))
        self.assertEqual(features.loc[0, "host_total_bp"], 8)
        self.assertEqual(features.loc[0, "phage_total_bp"], 8)
        self.assertIn("host_kmer_2_AC_fraction", features.columns)
        self.assertIn("phage_kmer_2_AA_fraction", features.columns)


if __name__ == "__main__":
    unittest.main()
