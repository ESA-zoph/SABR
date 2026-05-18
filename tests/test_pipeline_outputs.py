import unittest

from crispr_phage_predictor.io import FastaRecord
from crispr_phage_predictor.matching import SpacerHit
from crispr_phage_predictor.pipeline import (
    build_exact_match_heatmap,
    build_resistance_evidence_matrix,
)


class PipelineOutputTests(unittest.TestCase):
    def test_builds_exact_match_heatmap(self):
        bacteria = [
            FastaRecord("b1.fna", "c1", "c1", "ATGC"),
            FastaRecord("b2.fna", "c1", "c1", "ATGC"),
        ]
        phages = [
            FastaRecord("p1.fna", "p1", "p1", "ATGC"),
            FastaRecord("p2.fna", "p2", "p2", "ATGC"),
        ]
        hits = [
            SpacerHit(
                bacterium_id="b1.fna",
                array_id="array_1",
                phage_id="p1.fna",
                spacer_id="spacer_1",
                phage_contig_id="p1",
                start=1,
                end=4,
                strand="+",
                identity=1.0,
                mismatches=0,
                spacer_sequence="ATGC",
                protospacer_sequence="ATGC",
            )
        ]

        matrix = build_resistance_evidence_matrix(bacteria, phages, hits)
        heatmap = build_exact_match_heatmap(matrix)

        self.assertEqual(heatmap.loc["b1.fna", "p1.fna"], 1)
        self.assertEqual(heatmap.loc["b1.fna", "p2.fna"], 0)
        self.assertEqual(heatmap.loc["b2.fna", "p1.fna"], 0)
        self.assertEqual(heatmap.loc["b2.fna", "p2.fna"], 0)


if __name__ == "__main__":
    unittest.main()
