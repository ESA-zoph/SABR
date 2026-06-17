import tempfile
import unittest
from pathlib import Path

import pandas as pd

from crispr_phage_predictor.phage_annotation_features import (
    add_phage_annotation_features,
    phage_annotation_features_from_genbank,
)


GENBANK_TEXT = """LOCUS       TESTPHAGE              100 bp    DNA     linear   PHG 01-JAN-2000
DEFINITION  synthetic phage.
ACCESSION   TEST000001
VERSION     TEST000001.1
FEATURES             Location/Qualifiers
     source          1..100
                     /organism="synthetic phage"
     CDS             1..30
                     /product="integrase"
     CDS             31..60
                     /product="tail fiber protein"
     CDS             61..90
                     /product="holin"
ORIGIN
        1 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac
       51 gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt
//
"""


class PhageAnnotationFeatureTests(unittest.TestCase):
    def test_parses_keyword_features_from_genbank(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "TEST000001.gb"
            path.write_text(GENBANK_TEXT, encoding="utf-8")

            features = phage_annotation_features_from_genbank(path)

        self.assertEqual(features["phage_cds_count"], 3)
        self.assertEqual(features["phage_integrase_count"], 1)
        self.assertEqual(features["phage_tail_fiber_count"], 1)
        self.assertEqual(features["phage_holin_count"], 1)
        self.assertEqual(features["phage_temperate_marker_count"], 1)
        self.assertGreaterEqual(features["phage_structural_marker_count"], 1)

    def test_appends_annotation_features_by_accession(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            genbank_dir = root / "gb"
            genbank_dir.mkdir()
            (genbank_dir / "TEST000001.1.gb").write_text(GENBANK_TEXT, encoding="utf-8")
            table = pd.DataFrame([{"phage_accession": "TEST000001.1"}])

            augmented = add_phage_annotation_features(table, genbank_dir=genbank_dir)

        self.assertEqual(augmented.loc[0, "phage_cds_count"], 3)
        self.assertEqual(augmented.loc[0, "phage_integrase_count"], 1)


if __name__ == "__main__":
    unittest.main()
