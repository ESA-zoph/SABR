import unittest

from crispr_phage_predictor.crispr import detect_crispr_arrays


class CrisprDetectionTests(unittest.TestCase):
    def test_detects_exact_repeat_array(self):
        repeat = "GTTCACTGCCGTACAGGCAGCTTAGAAA"
        spacers = [
            "ACGTACGTACGTACGTACGTACGTACGTACGT",
            "TTGCAAGCTTAGGCTAATCGGATCCGATTAAC",
            "CGATCGATTTAAACCCGGGCGATCGATCGATC",
        ]
        sequence = "AAAACCCC" + repeat + repeat.join(spacers) + repeat + "GGGGTTTT"

        arrays = detect_crispr_arrays(sequence, genome_id="bacterium_a", contig_id="contig_1")

        self.assertEqual(len(arrays), 1)
        self.assertEqual(arrays[0].repeat_consensus, repeat)
        self.assertEqual(arrays[0].spacers, spacers)
        self.assertEqual(arrays[0].repeat_count, 4)
        self.assertEqual(arrays[0].spacer_count, 3)

    def test_ignores_two_repeat_candidate(self):
        repeat = "GTTCACTGCCGTACAGGCAGCTTAGAAA"
        spacer = "ACGTACGTACGTACGTACGTACGTACGTACGT"
        sequence = repeat + spacer + repeat

        arrays = detect_crispr_arrays(sequence, genome_id="bacterium_a", contig_id="contig_1")

        self.assertEqual(arrays, [])


if __name__ == "__main__":
    unittest.main()
