import unittest

from crispr_phage_predictor.matching import SpacerHit
from crispr_phage_predictor.scoring import score_crispr_targeting_evidence
from crispr_phage_predictor.scoring import score_experimental_pam_weighted_evidence
from crispr_phage_predictor.scoring import score_resistance_likelihood


def _hit(**overrides) -> SpacerHit:
    values = {
        "bacterium_id": "b.fna",
        "array_id": "array_1",
        "phage_id": "p.fna",
        "spacer_id": "spacer_1",
        "phage_contig_id": "p1",
        "start": 1,
        "end": 32,
        "strand": "+",
        "identity": 1.0,
        "mismatches": 0,
        "alignment_length": 32,
        "spacer_length": 32,
        "coverage": 1.0,
        "evalue": None,
        "bitscore": None,
        "spacer_sequence": "A" * 32,
        "protospacer_sequence": "A" * 32,
    }
    values.update(overrides)
    return SpacerHit(**values)


class TargetingScoringTests(unittest.TestCase):
    def test_scores_no_hits_as_zero(self):
        score = score_crispr_targeting_evidence([])

        self.assertEqual(score.score, 0.0)
        self.assertEqual(score.evidence_level, "no spacer-match evidence")

    def test_pam_supported_seed_clean_hit_scores_strongly(self):
        score = score_crispr_targeting_evidence(
            [
                _hit(
                    pam_match=True,
                    pam_support_level="compatible",
                    seed_mismatches=0,
                    cas_subtype_confidence=1.0,
                )
            ]
        )

        self.assertGreaterEqual(score.score, 70.0)
        self.assertIn("PAM/PFS support", score.interpretation)

    def test_pam_unsupported_hit_is_penalized(self):
        supported = score_crispr_targeting_evidence(
            [
                _hit(
                    pam_match=True,
                    pam_support_level="compatible",
                    seed_mismatches=0,
                    cas_subtype_confidence=1.0,
                )
            ]
        )
        unsupported = score_crispr_targeting_evidence(
            [
                _hit(
                    pam_match=False,
                    pam_support_level="not_supported",
                    seed_mismatches=0,
                    cas_subtype_confidence=1.0,
                )
            ]
        )

        self.assertLess(unsupported.score, supported.score)
        self.assertLess(unsupported.score, 50.0)
        self.assertIn("did not support", unsupported.interpretation)

    def test_pam_unsupported_multi_spacer_evidence_is_capped_low(self):
        score = score_crispr_targeting_evidence(
            [
                _hit(
                    spacer_id=f"spacer_{index}",
                    pam_match=False,
                    pam_support_level="not_supported",
                    seed_mismatches=0,
                    cas_subtype_confidence=1.0,
                )
                for index in range(1, 4)
            ]
        )

        self.assertLessEqual(score.score, 39.0)
        self.assertEqual(score.evidence_level, "weak candidate CRISPR targeting evidence")
        self.assertIn("did not support", score.interpretation)

    def test_experimental_pam_weighted_score_uses_partial_pam_compatibility(self):
        perfect = score_experimental_pam_weighted_evidence(
            [
                _hit(
                    pam_match=True,
                    pam_support_level="compatible",
                    pam_compatibility_score=1.0,
                    seed_mismatches=0,
                    cas_subtype_confidence=1.0,
                )
            ]
        )
        partial = score_experimental_pam_weighted_evidence(
            [
                _hit(
                    pam_match=False,
                    pam_support_level="not_supported",
                    pam_compatibility_score=0.5,
                    seed_mismatches=0,
                    cas_subtype_confidence=1.0,
                )
            ]
        )
        weak = score_experimental_pam_weighted_evidence(
            [
                _hit(
                    pam_match=False,
                    pam_support_level="not_supported",
                    pam_compatibility_score=1 / 3,
                    seed_mismatches=0,
                    cas_subtype_confidence=1.0,
                )
            ]
        )

        self.assertGreater(perfect, partial)
        self.assertGreater(partial, weak)
        self.assertLessEqual(weak, 39.0)

    def test_old_resistance_name_remains_alias(self):
        self.assertEqual(
            score_resistance_likelihood([]).score,
            score_crispr_targeting_evidence([]).score,
        )


if __name__ == "__main__":
    unittest.main()
