import tempfile
import unittest
from pathlib import Path

from crispr_phage_predictor.interactions import (
    empty_interaction_table,
    eop_class_from_value,
    load_interaction_table,
    normalize_interaction_table,
    parse_eop,
    susceptibility_from_eop_class,
    validate_interaction_table,
)


def _valid_row(interaction_id: str = "pa14_jbd18_eop") -> dict[str, str]:
    return {
        "interaction_id": interaction_id,
        "source_key": "Cady2012_PA14_CRISPR",
        "source_type": "paper",
        "pmid": "",
        "doi": "",
        "assay_type": "eop",
        "bacterium": "Pseudomonas aeruginosa",
        "strain": "PA14",
        "bacterial_accession": "NC_008463.1",
        "phage": "JBD18",
        "phage_accession": "JX495041.1",
        "reference_host": "",
        "raw_eop": "<1e-6",
        "eop_relation": "<",
        "eop_value": "1e-6",
        "eop_class": "none",
        "susceptibility_label": "resistant",
        "plaque_result": "no_plaques",
        "anti_crispr_status": "absent",
        "anti_crispr_genes": "",
        "crispr_interference_evidence": "experimental",
        "other_defense_evidence": "",
        "experimental_conditions": "",
        "curation_status": "curated",
        "curation_confidence": "high",
        "notes": "",
    }


class InteractionSchemaTests(unittest.TestCase):
    def test_validates_interaction_table(self):
        table = empty_interaction_table()
        table.loc[0] = _valid_row()

        validate_interaction_table(table)

    def test_rejects_duplicate_interaction_ids(self):
        table = empty_interaction_table()
        table.loc[0] = _valid_row("pair_1")
        table.loc[1] = _valid_row("pair_1")

        with self.assertRaisesRegex(ValueError, "duplicate interaction_id"):
            validate_interaction_table(table)

    def test_rejects_invalid_susceptibility_label(self):
        table = empty_interaction_table()
        row = _valid_row()
        row["susceptibility_label"] = "maybe"
        table.loc[0] = row

        with self.assertRaisesRegex(ValueError, "invalid susceptibility_label"):
            validate_interaction_table(table)

    def test_loads_tsv_table(self):
        table = empty_interaction_table()
        table.loc[0] = _valid_row()

        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "interactions.tsv"
            table.to_csv(path, sep="\t", index=False)
            loaded = load_interaction_table(path)

        self.assertEqual(len(loaded), 1)
        self.assertEqual(loaded.loc[0, "interaction_id"], "pa14_jbd18_eop")

    def test_parses_common_eop_values(self):
        self.assertEqual(parse_eop("<1e-6"), ("<", 1e-6))
        self.assertEqual(parse_eop("0.25"), ("=", 0.25))
        self.assertEqual(parse_eop("0.04±0.02"), ("=", 0.04))
        self.assertEqual(parse_eop("not reported"), ("not_reported", None))

    def test_maps_eop_to_class_and_label(self):
        self.assertEqual(eop_class_from_value(0.8), "high")
        self.assertEqual(eop_class_from_value(0.2), "medium")
        self.assertEqual(eop_class_from_value(0.01), "low")
        self.assertEqual(eop_class_from_value(1e-5), "trace")
        self.assertEqual(eop_class_from_value(0), "none")
        self.assertEqual(eop_class_from_value(0.001, "<"), "none")
        self.assertEqual(susceptibility_from_eop_class("low"), "reduced_susceptibility")
        self.assertEqual(susceptibility_from_eop_class("trace"), "resistant")

    def test_normalizes_missing_eop_fields(self):
        table = empty_interaction_table()
        row = _valid_row()
        row["raw_eop"] = "0.02"
        row["eop_relation"] = ""
        row["eop_value"] = ""
        row["eop_class"] = ""
        row["susceptibility_label"] = "unknown"
        table.loc[0] = row

        normalized = normalize_interaction_table(table)

        self.assertEqual(normalized.loc[0, "eop_relation"], "=")
        self.assertEqual(normalized.loc[0, "eop_value"], "0.02")
        self.assertEqual(normalized.loc[0, "eop_class"], "low")
        self.assertEqual(
            normalized.loc[0, "susceptibility_label"], "reduced_susceptibility"
        )


if __name__ == "__main__":
    unittest.main()
