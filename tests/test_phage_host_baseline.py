import unittest

import pandas as pd

from crispr_phage_predictor.phage_host_baseline import (
    evaluate_baseline_models,
    numeric_feature_columns,
)


def _feature_table() -> pd.DataFrame:
    rows = []
    for index in range(12):
        susceptible = index % 2 == 0
        rows.append(
            {
                "interaction_id": f"pair_{index}",
                "source_key": "source_a" if index < 6 else "source_b",
                "bacterium": "Bacterium",
                "strain": f"strain_{index}",
                "phage": f"phage_{index // 3}",
                "eop_class": "high" if susceptible else "none",
                "susceptibility_label": "susceptible" if susceptible else "resistant",
                "binary_susceptibility": "susceptible" if susceptible else "resistant",
                "eop_value": 1.0 if susceptible else 0.0,
                "host_linkage_status": "reference_proxy",
                "host_accession": "NC_1",
                "host_local_path": "host.fasta",
                "phage_linkage_status": "exact",
                "phage_accession": "NC_2",
                "phage_local_path": "phage.fasta",
                "uses_reference_proxy_host": True,
                "host_total_bp": 1000 + index,
                "phage_total_bp": 100 + index,
                "host_gc_percent": 50.0,
                "phage_gc_percent": 40.0 + index,
                "phage_host_gc_delta": -10.0 + index,
            }
        )
    return pd.DataFrame(rows)


class PhageHostBaselineTests(unittest.TestCase):
    def test_numeric_feature_columns_excludes_metadata(self):
        features = numeric_feature_columns(_feature_table())

        self.assertIn("host_total_bp", features)
        self.assertIn("uses_reference_proxy_host", features)
        self.assertNotIn("phage", features)
        self.assertNotIn("binary_susceptibility", features)

    def test_evaluates_baseline_models(self):
        results = evaluate_baseline_models(
            _feature_table(),
            split_strategy="row_random",
            test_size=0.25,
            random_state=1,
        )

        self.assertEqual({result.method for result in results}, {
            "majority_baseline",
            "logistic_regression",
            "random_forest",
            "extra_trees",
        })
        self.assertTrue(all(result.test_size == 3 for result in results))


if __name__ == "__main__":
    unittest.main()
