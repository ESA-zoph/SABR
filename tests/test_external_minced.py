import unittest

from crispr_phage_predictor.external.minced import parse_minced_output
from crispr_phage_predictor.io import FastaRecord


class MincedParserTests(unittest.TestCase):
    def test_parses_minced_text_output(self):
        record = FastaRecord(
            source_file="bacterium.fna",
            record_id="contig_1",
            description="contig_1",
            sequence="A" * 500,
        )
        output = """
Sequence 'record_1' (500 bp)
CRISPR 1   Range: 10 - 180
POSITION      REPEAT                         SPACER
10            GTTCACTGCCGTACAGGCAGCTTAGAAA  ACGTACGTACGTACGTACGTACGTACGTACGT
71            GTTCACTGCCGTACAGGCAGCTTAGAAA  TTGCAAGCTTAGGCTAATCGGATCCGATTAAC
132           GTTCACTGCCGTACAGGCAGCTTAGAAA
Repeats: 3
"""

        arrays = parse_minced_output(output, {"record_1": record})

        self.assertEqual(len(arrays), 1)
        self.assertEqual(arrays[0].genome_id, "bacterium.fna")
        self.assertEqual(arrays[0].contig_id, "contig_1")
        self.assertEqual(arrays[0].start, 10)
        self.assertEqual(arrays[0].end, 180)
        self.assertEqual(arrays[0].repeat_consensus, "GTTCACTGCCGTACAGGCAGCTTAGAAA")
        self.assertEqual(len(arrays[0].spacers), 2)


if __name__ == "__main__":
    unittest.main()

