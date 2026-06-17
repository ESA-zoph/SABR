import tempfile
import unittest
from pathlib import Path

import pandas as pd

from crispr_phage_predictor.crispr import CrisprArray
from crispr_phage_predictor.interaction_targeting_features import (
    add_targeting_features,
    find_fuzzy_spacer_hits,
)
from crispr_phage_predictor.io import FastaRecord


class InteractionTargetingFeatureTests(unittest.TestCase):
    def test_adds_targeting_features_from_linked_fastas(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            spacer = "ACGTACGTACGTACGTACGTACGTACGTAC"
            repeat = "GTTCACTGCCGTACAGGCAGCTTAGAAA"
            host_sequence = repeat + spacer + repeat + spacer + repeat
            host_path = root / "host.fasta"
            phage_path = root / "phage.fasta"
            host_path.write_text(f">host\n{host_sequence}\n", encoding="utf-8")
            phage_path.write_text(f">phage\nTTTT{spacer}GGGG\n", encoding="utf-8")
            table = pd.DataFrame(
                [
                    {
                        "interaction_id": "pair_1",
                        "host_local_path": str(host_path),
                        "phage_local_path": str(phage_path),
                    }
                ]
            )

            augmented = add_targeting_features(table)

        self.assertEqual(augmented.loc[0, "crispr_array_count"], 1)
        self.assertEqual(augmented.loc[0, "crispr_spacer_count"], 2)
        self.assertGreaterEqual(augmented.loc[0, "spacer_hit_count"], 1)
        self.assertGreater(augmented.loc[0, "crispr_targeting_score"], 0)
        self.assertGreaterEqual(augmented.loc[0, "fuzzy_spacer_candidate_count"], 1)
        self.assertGreater(augmented.loc[0, "graded_crispr_interference_score"], 0)

    def test_fuzzy_matching_allows_distal_mismatch_with_seed_intact(self):
        spacer = "ACGTACGTACGTACGTACGTACGTACGTAC"
        protospacer = "ACGTACGTACGTACGTACGTACGTACGTTC"
        arrays = [
            CrisprArray(
                array_id="array_1",
                genome_id="host",
                contig_id="host_contig",
                start=1,
                end=100,
                repeat_consensus="GTTCACTGCCGTACAGGCAGCTTAGAAA",
                spacers=[spacer],
            )
        ]
        records = [
            FastaRecord(
                source_file="phage.fasta",
                record_id="phage",
                description="phage",
                sequence=f"TTTT{protospacer}GGGG",
            )
        ]

        hits = find_fuzzy_spacer_hits(arrays, records, max_mismatches=2)

        self.assertGreaterEqual(len(hits), 1)
        self.assertEqual(hits[0].mismatches, 1)
        self.assertEqual(hits[0].seed_edge_mismatches, 0)
        self.assertGreater(hits[0].identity, 0.95)


if __name__ == "__main__":
    unittest.main()
