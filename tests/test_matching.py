import unittest

from crispr_phage_predictor.crispr import CrisprArray
from crispr_phage_predictor.io import FastaRecord
from crispr_phage_predictor.matching import find_spacer_hits, reverse_complement


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


if __name__ == "__main__":
    unittest.main()
