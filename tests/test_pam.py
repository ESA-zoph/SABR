import unittest

from crispr_phage_predictor.pam import evaluate_pam_rule


class PamEvaluationTests(unittest.TestCase):
    def test_matches_3prime_iupac_pam(self):
        result = evaluate_pam_rule(
            protospacer_5p_flank="TTTT",
            protospacer_3p_flank="AGGCT",
            pam_rule="3prime:NGG",
        )

        self.assertEqual(result.pam_sequence, "AGG")
        self.assertTrue(result.pam_match)
        self.assertEqual(result.pam_support_level, "compatible")
        self.assertEqual(result.compatibility_score, 1.0)

    def test_matches_5prime_iupac_pam_adjacent_to_protospacer(self):
        result = evaluate_pam_rule(
            protospacer_5p_flank="TTAAG",
            protospacer_3p_flank="CCCC",
            pam_rule="5prime:AWG",
        )

        self.assertEqual(result.pam_sequence, "AAG")
        self.assertTrue(result.pam_match)
        self.assertEqual(result.compatibility_score, 1.0)

    def test_matches_genomic_3prime_pam(self):
        result = evaluate_pam_rule(
            protospacer_5p_flank="CATGCTTGTGATCGAGGACC",
            protospacer_3p_flank="GCACCGGCACCACTGGCAGC",
            genomic_upstream_flank="GCTGCCAGTGGTGCCGGTGC",
            genomic_downstream_flank="GGTCCTCGATCACAAGCATG",
            pam_rule="genomic_3prime:GG",
        )

        self.assertEqual(result.pam_sequence, "GG")
        self.assertTrue(result.pam_match)

    def test_reports_unsupported_pam(self):
        result = evaluate_pam_rule(
            protospacer_5p_flank="TTTT",
            protospacer_3p_flank="ATG",
            pam_rule="3prime:NGG",
        )

        self.assertEqual(result.pam_sequence, "ATG")
        self.assertFalse(result.pam_match)
        self.assertEqual(result.pam_support_level, "not_supported")
        self.assertEqual(result.compatibility_score, round(2 / 3, 6))

    def test_scores_fully_incompatible_pam_as_zero(self):
        result = evaluate_pam_rule(
            protospacer_5p_flank="TTTT",
            protospacer_3p_flank="AAA",
            pam_rule="3prime:NGG",
        )

        self.assertEqual(result.pam_sequence, "AAA")
        self.assertFalse(result.pam_match)
        self.assertEqual(result.compatibility_score, round(1 / 3, 6))

    def test_reports_insufficient_flank(self):
        result = evaluate_pam_rule(
            protospacer_5p_flank="A",
            protospacer_3p_flank="",
            pam_rule="5prime:AWG",
        )

        self.assertEqual(result.pam_sequence, "A")
        self.assertIsNone(result.pam_match)
        self.assertEqual(result.pam_support_level, "insufficient_flank")
        self.assertIsNone(result.compatibility_score)

    def test_reports_invalid_rule(self):
        result = evaluate_pam_rule(
            protospacer_5p_flank="AAA",
            protospacer_3p_flank="CCC",
            pam_rule="middle:NGG",
        )

        self.assertEqual(result.pam_sequence, "")
        self.assertIsNone(result.pam_match)
        self.assertEqual(result.pam_support_level, "invalid_rule")


if __name__ == "__main__":
    unittest.main()
