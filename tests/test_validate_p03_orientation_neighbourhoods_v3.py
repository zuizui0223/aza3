#!/usr/bin/env python3
from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))

import validate_p03_orientation_neighbourhoods_v3 as target


class P03NeighbourhoodTests(unittest.TestCase):
    def test_registry_is_fail_closed(self) -> None:
        rows = target.validate_registry()
        self.assertEqual(len(rows), 4)
        self.assertEqual(sum(r["status"] == "BACKBONE_CONTRAST_READY" for r in rows), 2)
        self.assertEqual(rows[2]["status"], "TRAIT_GAP_FIRST")
        self.assertEqual(rows[3]["status"], "OPEN_DO_NOT_FORCE")

    def test_jpn36_is_not_auto_counted(self) -> None:
        rows = target.validate_registry()
        self.assertIn("JPN36", rows[3]["why_not_double_count"])


if __name__ == "__main__":
    unittest.main()
