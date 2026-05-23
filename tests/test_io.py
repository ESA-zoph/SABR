import unittest

from crispr_phage_predictor.io import (
    FastaRecord,
    deduplicate_records,
    extract_accession,
    sequence_hash,
    summarize_accession_conflicts,
    summarize_duplicate_records,
)


class FastaDeduplicationTests(unittest.TestCase):
    def test_extracts_ncbi_accessions(self):
        self.assertEqual(extract_accession("JX495041.1", ""), "JX495041.1")
        self.assertEqual(extract_accession("gi|1|ref|NC_008463.1|", ""), "NC_008463.1")
        self.assertEqual(
            extract_accession("seq1", "AF085222.2 Streptococcus phage DT1"),
            "AF085222.2",
        )

    def test_sequence_hash_normalizes_case_and_whitespace(self):
        self.assertEqual(sequence_hash("acgt\nacgt"), sequence_hash("ACGTACGT"))

    def test_deduplicates_records_by_sequence_hash(self):
        records = [
            FastaRecord("first.fna", "r1", "r1", "ACGTACGT"),
            FastaRecord("duplicate.fna", "r2", "r2", "acgtacgt"),
            FastaRecord("second.fna", "r3", "r3", "TTTTCCCC"),
        ]

        unique = deduplicate_records(records)

        self.assertEqual([record.source_file for record in unique], ["first.fna", "second.fna"])

    def test_summarizes_duplicate_records(self):
        records = [
            FastaRecord("first.fna", "r1", "r1", "ACGTACGT"),
            FastaRecord("duplicate.fna", "r2", "r2", "ACGTACGT"),
        ]

        summary = summarize_duplicate_records(records)

        self.assertEqual(len(summary), 1)
        self.assertEqual(summary.loc[0, "kept_source_file"], "first.fna")
        self.assertEqual(summary.loc[0, "duplicate_count"], 1)
        self.assertIn("duplicate.fna:r2", summary.loc[0, "duplicate_records"])

    def test_summarizes_accession_conflicts(self):
        records = [
            FastaRecord("first.fna", "NC_000001.1", "NC_000001.1 record", "ACGT"),
            FastaRecord("second.fna", "NC_000001.1", "NC_000001.1 record", "TGCA"),
        ]

        summary = summarize_accession_conflicts(records)

        self.assertEqual(len(summary), 1)
        self.assertEqual(summary.loc[0, "accession"], "NC_000001.1")
        self.assertEqual(summary.loc[0, "distinct_sequence_hashes"], 2)


if __name__ == "__main__":
    unittest.main()
