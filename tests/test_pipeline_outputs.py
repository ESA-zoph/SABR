import unittest

from crispr_phage_predictor.io import FastaRecord
from crispr_phage_predictor.cas_prediction import ArrayCasPrediction
from crispr_phage_predictor.matching import SpacerHit
from crispr_phage_predictor.pipeline import (
    annotate_spacer_hits_with_pam,
    build_crispr_targeting_evidence_matrix,
    build_exact_match_heatmap,
    build_resistance_evidence_matrix,
    summarize_pam_subtype_support,
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
                alignment_length=4,
                spacer_length=4,
                coverage=1.0,
                evalue=None,
                bitscore=None,
                spacer_sequence="ATGC",
                protospacer_sequence="ATGC",
            )
        ]

        matrix = build_resistance_evidence_matrix(bacteria, phages, hits)
        heatmap = build_exact_match_heatmap(matrix)

        self.assertIn("crispr_targeting_score", matrix.columns)
        self.assertIn("experimental_pam_weighted_score", matrix.columns)
        self.assertIn("hypothetical_resistance_score", matrix.columns)
        self.assertEqual(heatmap.loc["b1.fna", "p1.fna"], 1)
        self.assertEqual(heatmap.loc["b1.fna", "p2.fna"], 0)
        self.assertEqual(heatmap.loc["b2.fna", "p1.fna"], 0)
        self.assertEqual(heatmap.loc["b2.fna", "p2.fna"], 0)

    def test_rolls_pam_support_into_evidence_matrix(self):
        bacteria = [FastaRecord("b1.fna", "c1", "c1", "ATGC")]
        phages = [FastaRecord("p1.fna", "p1", "p1", "ATGC")]
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
                alignment_length=4,
                spacer_length=4,
                coverage=1.0,
                evalue=None,
                bitscore=None,
                spacer_sequence="ATGC",
                protospacer_sequence="ATGC",
                protospacer_3p_flank="AGGTT",
            )
        ]

        annotated_hits = annotate_spacer_hits_with_pam(hits, default_pam_rule="3prime:NGG")
        matrix = build_crispr_targeting_evidence_matrix(bacteria, phages, annotated_hits)

        self.assertEqual(annotated_hits[0].pam_sequence, "AGG")
        self.assertTrue(annotated_hits[0].pam_match)
        self.assertEqual(annotated_hits[0].pam_compatibility_score, 1.0)
        self.assertEqual(matrix.loc[0, "pam_supported_hits"], 1)
        self.assertEqual(matrix.loc[0, "pam_support_level"], "compatible")
        self.assertEqual(matrix.loc[0, "best_pam_compatibility_score"], 1.0)
        self.assertEqual(matrix.loc[0, "mean_pam_compatibility_score"], 1.0)
        self.assertEqual(matrix.loc[0, "seed_evaluated_hits"], 1)
        self.assertEqual(matrix.loc[0, "best_seed_mismatches"], 0)

    def test_uses_cas_prediction_pam_rule_without_manual_rule(self):
        hits = [
            SpacerHit(
                bacterium_id="b1.fna",
                array_id="array_1",
                phage_id="p1.fna",
                spacer_id="spacer_1",
                phage_contig_id="p1",
                start=10,
                end=13,
                strand="+",
                identity=1.0,
                mismatches=0,
                alignment_length=4,
                spacer_length=4,
                coverage=1.0,
                evalue=None,
                bitscore=None,
                spacer_sequence="ATGC",
                protospacer_sequence="ATGC",
                protospacer_5p_flank="AACCC",
            )
        ]
        cas_predictions = {
            "array_1": ArrayCasPrediction(
                array_id="array_1",
                cas_subtype="I-F",
                cas_subtype_confidence=1.0,
                prediction_method="nearest_repeat",
                pam_rule="genomic_3prime:GG",
                pam_rule_source="curated_subtype_catalog",
            )
        }

        annotated_hits = annotate_spacer_hits_with_pam(
            hits,
            cas_predictions_by_array=cas_predictions,
        )

        self.assertEqual(annotated_hits[0].predicted_cas_subtype, "I-F")
        self.assertEqual(annotated_hits[0].pam_rule, "genomic_3prime:GG")
        self.assertEqual(annotated_hits[0].pam_sequence, "")
        self.assertIsNone(annotated_hits[0].pam_match)
        self.assertEqual(annotated_hits[0].pam_support_level, "insufficient_flank")

    def test_uses_cas_prediction_downstream_type_if_pam_rule(self):
        hits = [
            SpacerHit(
                bacterium_id="b1.fna",
                array_id="array_1",
                phage_id="p1.fna",
                spacer_id="spacer_1",
                phage_contig_id="p1",
                start=10,
                end=13,
                strand="+",
                identity=1.0,
                mismatches=0,
                alignment_length=4,
                spacer_length=4,
                coverage=1.0,
                evalue=None,
                bitscore=None,
                spacer_sequence="ATGC",
                protospacer_sequence="ATGC",
                protospacer_3p_flank="TTAAA",
                genomic_downstream_flank="GGAAA",
            )
        ]
        cas_predictions = {
            "array_1": ArrayCasPrediction(
                array_id="array_1",
                cas_subtype="I-F",
                cas_subtype_confidence=1.0,
                prediction_method="nearest_repeat",
                pam_rule="genomic_3prime:GG",
                pam_rule_source="curated_subtype_catalog",
            )
        }

        annotated_hits = annotate_spacer_hits_with_pam(
            hits,
            cas_predictions_by_array=cas_predictions,
        )

        self.assertEqual(annotated_hits[0].pam_sequence, "GG")
        self.assertTrue(annotated_hits[0].pam_match)
        self.assertEqual(annotated_hits[0].seed_region, "3prime:1-4")
        self.assertEqual(annotated_hits[0].seed_mismatches, 0)

    def test_summarizes_exploratory_pam_subtype_support(self):
        hits = [
            SpacerHit(
                bacterium_id="b1.fna",
                array_id="array_1",
                phage_id="p1.fna",
                spacer_id="spacer_1",
                phage_contig_id="p1",
                start=10,
                end=13,
                strand="+",
                identity=1.0,
                mismatches=0,
                alignment_length=4,
                spacer_length=4,
                coverage=1.0,
                evalue=None,
                bitscore=None,
                spacer_sequence="ATGC",
                protospacer_sequence="ATGC",
                genomic_downstream_flank="GGAAA",
                predicted_cas_subtype="I-F",
                cas_subtype_confidence=1.0,
            )
        ]

        table = summarize_pam_subtype_support(hits)

        self.assertEqual(table.loc[0, "top_pam_supported_subtype"], "I-F")
        self.assertEqual(table.loc[0, "top_pam_support_count"], 1)
        self.assertEqual(table.loc[0, "repeat_pam_subtype_agreement"], "agrees_top_subtype")

    def test_flags_repeat_pam_subtype_conflict_diagnostically(self):
        hits = [
            SpacerHit(
                bacterium_id="b1.fna",
                array_id="array_1",
                phage_id="p1.fna",
                spacer_id="spacer_1",
                phage_contig_id="p1",
                start=10,
                end=13,
                strand="+",
                identity=1.0,
                mismatches=0,
                alignment_length=4,
                spacer_length=4,
                coverage=1.0,
                evalue=None,
                bitscore=None,
                spacer_sequence="ATGC",
                protospacer_sequence="ATGC",
                genomic_downstream_flank="GGAAA",
                predicted_cas_subtype="II-A",
                cas_subtype_confidence=1.0,
            )
        ]

        table = summarize_pam_subtype_support(hits)

        self.assertEqual(table.loc[0, "top_pam_supported_subtype"], "I-F")
        self.assertEqual(
            table.loc[0, "repeat_pam_subtype_agreement"],
            "repeat_prediction_not_pam_supported",
        )


if __name__ == "__main__":
    unittest.main()
