import unittest

from crispr_phage_predictor.crispr import CrisprArray
from crispr_phage_predictor.io import FastaRecord
from crispr_phage_predictor.matching import (
    extract_protospacer_context,
    find_spacer_hits,
    reverse_complement,
    summarize_seed_mismatches,
)


class SpacerMatchingTests(unittest.TestCase):
    def test_finds_forward_exact_spacer_hit(self):
        array = CrisprArray(
            array_id="bacterium.fna|contig_1|array_1",
            genome_id="bacterium.fna",
            contig_id="contig_1",
            start=1,
            end=100,
            repeat_consensus="GTTCACTGCCGTACAGGCAGCTTAGAAA",
            spacers=["ACGTACGTACGTACGTACGTACGTACGTACGT"],
        )
        phage = FastaRecord(
            source_file="phage.fna",
            record_id="phage_contig_1",
            description="phage_contig_1",
            sequence="TTTT" + array.spacers[0] + "GGGG",
        )

        hits = find_spacer_hits([array], [phage])

        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0].strand, "+")
        self.assertEqual(hits[0].start, 5)
        self.assertEqual(hits[0].identity, 1.0)
        self.assertEqual(hits[0].protospacer_5p_flank, "TTTT")
        self.assertEqual(hits[0].protospacer_3p_flank, "GGGG")

    def test_finds_reverse_complement_exact_spacer_hit(self):
        spacer = "ACGTACGTACGTACGTACGTACGTACGTAACC"
        array = CrisprArray(
            array_id="bacterium.fna|contig_1|array_1",
            genome_id="bacterium.fna",
            contig_id="contig_1",
            start=1,
            end=100,
            repeat_consensus="GTTCACTGCCGTACAGGCAGCTTAGAAA",
            spacers=[spacer],
        )
        phage = FastaRecord(
            source_file="phage.fna",
            record_id="phage_contig_1",
            description="phage_contig_1",
            sequence="TTTT" + reverse_complement(spacer) + "GGGG",
        )

        hits = find_spacer_hits([array], [phage])

        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0].strand, "-")
        self.assertEqual(hits[0].protospacer_5p_flank, "CCCC")
        self.assertEqual(hits[0].protospacer_3p_flank, "AAAA")

    def test_extracts_reverse_strand_flanks_in_protospacer_orientation(self):
        context = extract_protospacer_context(
            phage_sequence="AAAACCCCGGGGTTTT",
            start=5,
            end=8,
            strand="-",
            flank_length=4,
        )

        self.assertEqual(context.genomic_upstream_flank, "AAAA")
        self.assertEqual(context.genomic_downstream_flank, "GGGG")
        self.assertEqual(context.protospacer_5p_flank, "CCCC")
        self.assertEqual(context.protospacer_3p_flank, "TTTT")

    def test_summarizes_5prime_pam_proximal_seed_mismatches(self):
        summary = summarize_seed_mismatches(
            spacer_sequence="AACCGGTT",
            protospacer_sequence="ATCCGGTA",
            pam_rule="5prime:CCN",
            seed_length=4,
        )

        self.assertIsNotNone(summary)
        assert summary is not None
        self.assertEqual(summary.seed_region, "5prime:1-4")
        self.assertEqual(summary.seed_mismatches, 1)
        self.assertEqual(summary.seed_mismatch_positions, "2")

    def test_summarizes_3prime_pam_proximal_seed_mismatches(self):
        summary = summarize_seed_mismatches(
            spacer_sequence="AACCGGTT",
            protospacer_sequence="ATCCGGTA",
            pam_rule="3prime:NGG",
            seed_length=4,
        )

        self.assertIsNotNone(summary)
        assert summary is not None
        self.assertEqual(summary.seed_region, "3prime:5-8")
        self.assertEqual(summary.seed_mismatches, 1)
        self.assertEqual(summary.seed_mismatch_positions, "4")


if __name__ == "__main__":
    unittest.main()
