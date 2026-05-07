from __future__ import annotations

import unittest
from pathlib import Path

from app.analysis_api import run_analysis


class TestAnalysisAPI(unittest.TestCase):
    def test_run_analysis_linking_number_for_2sap_fixture(self) -> None:
        path = Path(__file__).resolve().parents[3] / "data" / "CreatorAll" / "All_2SAPs" / "All2SAPsL2M1span2num1.txt"
        result = run_analysis(path, "linking_number")
        self.assertEqual(result["action"], "linking_number")
        self.assertGreaterEqual(result["metrics"]["objects"], 1)
        self.assertEqual(result["metrics"]["reported_objects"], len(result["rows"]))
        self.assertTrue(all("linking_number" in row or "error" in row for row in result["rows"]))
        self.assertIn("linking_number", result["rows"][0])
