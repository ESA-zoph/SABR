import unittest
from unittest.mock import patch

from crispr_phage_predictor.ml.check_cctyper_environment import check_cctyper_environment


class CheckCCTyperEnvironmentTests(unittest.TestCase):
    def test_reports_missing_commands_and_database(self):
        with (
            patch("shutil.which", return_value=None),
            patch("importlib.util.find_spec", return_value=None),
            patch.dict("os.environ", {}, clear=True),
        ):
            checks = check_cctyper_environment()

        names = {check.name: check for check in checks}

        self.assertFalse(names["cctyper"].ok)
        self.assertFalse(names["prodigal"].ok)
        self.assertFalse(names["CCTYPER_DB"].ok)


if __name__ == "__main__":
    unittest.main()
