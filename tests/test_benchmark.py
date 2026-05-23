import tempfile
import unittest
from pathlib import Path

from crispr_phage_predictor.benchmark import (
    empty_benchmark_label_table,
    evaluate_benchmark_run,
    load_benchmark_label_table,
    summarize_benchmark_evaluation,
    validate_benchmark_label_table,
)


def _valid_row(pair_id: str = "pair_1") -> dict[str, str]:
    return {
        "pair_id": pair_id,
        "label_version": "v0.1",
        "label_status": "validated",
        "benchmark_split": "validation",
        "bacterium": "Pseudomonas aeruginosa",
        "strain": "PA14",
        "bacterial_accession": "NC_008463.1",
        "phage": "JBD18",
        "phage_accession": "JX495041.1",
        "local_bacterium_file": "Paeruginosa_PA14.fasta",
        "local_phage_file": "PhageJBD18.fasta",
        "phenotype_label": "resistant",
        "crispr_resistance_label": "crispr_resistant",
        "crispr_evidence_level": "experimental",
        "pam_evidence_level": "validated_pam",
        "anti_crispr_status": "absent",
        "host_range_status": "host",
        "expected_sabr_behavior": "high_score_expected",
        "source_keys": "Cady2012_PA14_CRISPR",
        "curation_confidence": "high",
        "notes": "",
    }


class BenchmarkSchemaTests(unittest.TestCase):
    def test_validates_benchmark_table(self):
        table = empty_benchmark_label_table()
        table.loc[0] = _valid_row()

        validate_benchmark_label_table(table)

    def test_rejects_invalid_enum_value(self):
        table = empty_benchmark_label_table()
        row = _valid_row()
        row["phenotype_label"] = "maybe"
        table.loc[0] = row

        with self.assertRaisesRegex(ValueError, "invalid phenotype_label"):
            validate_benchmark_label_table(table)

    def test_rejects_duplicate_pair_ids(self):
        table = empty_benchmark_label_table()
        table.loc[0] = _valid_row("pair_1")
        table.loc[1] = _valid_row("pair_1")

        with self.assertRaisesRegex(ValueError, "duplicate pair_id"):
            validate_benchmark_label_table(table)

    def test_loads_tsv_table(self):
        table = empty_benchmark_label_table()
        table.loc[0] = _valid_row()

        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "benchmark.tsv"
            table.to_csv(path, sep="\t", index=False)
            loaded = load_benchmark_label_table(path)

        self.assertEqual(len(loaded), 1)
        self.assertEqual(loaded.loc[0, "pair_id"], "pair_1")

    def test_evaluates_benchmark_run_by_local_files(self):
        table = empty_benchmark_label_table()
        table.loc[0] = _valid_row()

        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            benchmark_path = root / "benchmark.tsv"
            table.to_csv(benchmark_path, sep="\t", index=False)
            run_dir = root / "run"
            run_dir.mkdir()
            (run_dir / "evidence_matrix.csv").write_text(
                "bacterium,phage,hypothetical_resistance_score,spacer_hits,"
                "unique_matching_spacers,best_identity_percent,best_coverage_percent,"
                "pam_support_level,current_evidence_level,interpretation\n"
                "Paeruginosa_PA14.fasta,PhageJBD18.fasta,82,5,5,100,100,"
                "compatible,strong candidate CRISPR targeting evidence,ok\n",
                encoding="utf-8",
            )
            (run_dir / "bacterial_records.csv").write_text(
                "source_file,accession\nPaeruginosa_PA14.fasta,NC_008463.1\n",
                encoding="utf-8",
            )
            (run_dir / "phage_records.csv").write_text(
                "source_file,accession\nPhageJBD18.fasta,JX495041.1\n",
                encoding="utf-8",
            )

            evaluation = evaluate_benchmark_run(run_dir, benchmark_path)

        self.assertEqual(len(evaluation), 1)
        self.assertEqual(evaluation.loc[0, "run_match_status"], "matched")
        self.assertEqual(evaluation.loc[0, "score_expectation_result"], "meets_expectation")
        self.assertEqual(evaluation.loc[0, "run_phage_accession"], "JX495041.1")

    def test_summarizes_benchmark_evaluation(self):
        table = empty_benchmark_label_table()
        table.loc[0] = _valid_row()

        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            benchmark_path = root / "benchmark.tsv"
            table.to_csv(benchmark_path, sep="\t", index=False)
            run_dir = root / "run"
            run_dir.mkdir()
            (run_dir / "evidence_matrix.csv").write_text(
                "bacterium,phage,unique_matching_spacers\n"
                "Paeruginosa_PA14.fasta,PhageJBD18.fasta,5\n",
                encoding="utf-8",
            )
            evaluation = evaluate_benchmark_run(run_dir, benchmark_path)
            summary = summarize_benchmark_evaluation(evaluation)

        self.assertEqual(summary.loc[0, "rows"], 1)
        self.assertEqual(summary.loc[0, "matched_rows"], 1)


if __name__ == "__main__":
    unittest.main()
