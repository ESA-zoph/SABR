import unittest

import pandas as pd

from crispr_phage_predictor.accession_linkage import (
    accession_coverage,
    empty_accession_linkage_table,
    validate_accession_linkage_table,
)


def _valid_row(linkage_id: str = "pa14") -> dict[str, str]:
    return {
        "linkage_id": linkage_id,
        "entity_type": "bacterium",
        "source_key": "Cady2012_PA14_CRISPR",
        "display_name": "Pseudomonas aeruginosa",
        "strain_or_isolate": "PA14",
        "accession": "NC_008463.1",
        "accession_database": "RefSeq",
        "assembly_level": "complete_genome",
        "sequence_status": "available",
        "linkage_status": "exact",
        "confidence": "high",
        "local_path": "",
        "notes": "",
    }


class AccessionLinkageTests(unittest.TestCase):
    def test_validates_linkage_table(self):
        table = empty_accession_linkage_table()
        table.loc[0] = _valid_row()

        validate_accession_linkage_table(table)

    def test_rejects_exact_linkage_without_accession(self):
        table = empty_accession_linkage_table()
        row = _valid_row()
        row["accession"] = ""
        table.loc[0] = row

        with self.assertRaisesRegex(ValueError, "requires accession"):
            validate_accession_linkage_table(table)

    def test_rejects_duplicate_linkage_id(self):
        table = empty_accession_linkage_table()
        table.loc[0] = _valid_row("dup")
        table.loc[1] = _valid_row("dup")

        with self.assertRaisesRegex(ValueError, "duplicate linkage_id"):
            validate_accession_linkage_table(table)

    def test_reports_pair_coverage_from_existing_accession_fields(self):
        interactions = pd.DataFrame(
            [
                {
                    "interaction_id": "pair_1",
                    "source_key": "source",
                    "bacterium": "Bacterium",
                    "strain": "strain",
                    "bacterial_accession": "NC_1",
                    "phage": "Phage",
                    "phage_accession": "NC_2",
                }
            ]
        )
        linkage = empty_accession_linkage_table()

        coverage = accession_coverage(interactions, linkage)

        self.assertTrue(bool(coverage.loc[0, "pair_genome_ready"]))
        self.assertTrue(bool(coverage.loc[0, "pair_hybrid_ready"]))
        self.assertEqual(coverage.loc[0, "dataset_tier"], "tier1_exact_pair")

    def test_reports_pair_coverage_from_linkage_manifest(self):
        interactions = pd.DataFrame(
            [
                {
                    "interaction_id": "pair_1",
                    "source_key": "source",
                    "bacterium": "Bacterium",
                    "strain": "strain",
                    "bacterial_accession": "",
                    "phage": "Phage",
                    "phage_accession": "",
                }
            ]
        )
        linkage = empty_accession_linkage_table()
        linkage.loc[0] = {
            **_valid_row("b"),
            "source_key": "source",
            "display_name": "Bacterium",
            "strain_or_isolate": "strain",
            "accession": "NC_1",
        }
        linkage.loc[1] = {
            **_valid_row("p"),
            "entity_type": "phage",
            "source_key": "source",
            "display_name": "Phage",
            "strain_or_isolate": "Phage",
            "accession": "NC_2",
        }

        coverage = accession_coverage(interactions, linkage)

        self.assertTrue(bool(coverage.loc[0, "pair_genome_ready"]))
        self.assertTrue(bool(coverage.loc[0, "pair_hybrid_ready"]))

    def test_reports_hybrid_pair_coverage_from_reference_proxy(self):
        interactions = pd.DataFrame(
            [
                {
                    "interaction_id": "pair_1",
                    "source_key": "source",
                    "bacterium": "Bacterium",
                    "strain": "panel_strain",
                    "bacterial_accession": "",
                    "phage": "Phage",
                    "phage_accession": "",
                }
            ]
        )
        linkage = empty_accession_linkage_table()
        linkage.loc[0] = {
            **_valid_row("b"),
            "source_key": "source",
            "display_name": "Bacterium",
            "strain_or_isolate": "*",
            "accession": "NC_1",
            "linkage_status": "reference_proxy",
        }
        linkage.loc[1] = {
            **_valid_row("p"),
            "entity_type": "phage",
            "source_key": "source",
            "display_name": "Phage",
            "strain_or_isolate": "Phage",
            "accession": "NC_2",
        }

        coverage = accession_coverage(interactions, linkage)

        self.assertFalse(bool(coverage.loc[0, "pair_genome_ready"]))
        self.assertTrue(bool(coverage.loc[0, "pair_hybrid_ready"]))
        self.assertEqual(coverage.loc[0, "dataset_tier"], "tier2_proxy_host_exact_phage")


if __name__ == "__main__":
    unittest.main()
