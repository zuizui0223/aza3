#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SUMMARY = ROOT / "data/planning/chapter3_exact_hole_fill_summary_v7.json"


def rows(name: str) -> list[dict[str, str]]:
    with (ROOT / "data/planning" / name).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


class ExactHoleFillV7Tests(unittest.TestCase):
    def test_snapshot_counts(self) -> None:
        s = json.loads(SUMMARY.read_text(encoding="utf-8"))
        self.assertEqual(s["source_snapshot"]["published_or_no_new_label_core_species"], 128)
        self.assertEqual(s["published_core_hole_fill"]["new_tree_samples_total"], 228)
        self.assertEqual(s["source_snapshot"]["unpublished_new_name_species"], 25)
        self.assertEqual(s["source_snapshot"]["provisional_name_species"], 1)

    def test_empty_and_broken_slots(self) -> None:
        p1 = rows("chapter3_holefill_priority1_empty_blocks_v7.csv")
        p2 = rows("chapter3_holefill_priority2_broken_tips_v7.csv")
        self.assertEqual(len(p1), 33)
        self.assertEqual(sum(int(r["new_tree_samples_needed"]) for r in p1), 66)
        self.assertEqual({r["species_binomial"] for r in p2}, {
            "Cirsium pendulum", "Cirsium sieboldii", "Cirsium nipponicum",
            "Cirsium kamtschaticum", "Cirsium buergeri", "Cirsium microspicatum",
        })

    def test_focal_species_are_nested_not_duplicated(self) -> None:
        hama = rows("chapter3_holefill_priority3_hamaazami_block_v7.csv")
        p6 = rows("chapter3_holefill_priority6_add_own_trait_link_v7.csv")
        self.assertTrue({"Cirsium brevicaule", "Cirsium irimtiense"}.issubset({r["species_binomial"] for r in hama}))
        self.assertTrue({"Cirsium dipsacolepis", "Cirsium lineare"}.issubset({r["species_binomial"] for r in p6}))

    def test_taxonomy_extension_and_conflicts_are_separate(self) -> None:
        ext = rows("chapter3_holefill_taxonomy_extension_v7.csv")
        conflicts = rows("chapter3_holefill_taxonomy_conflicts_v7.csv")
        self.assertEqual(len(ext), 26)
        self.assertEqual(sum(int(r["new_samples_if_extension_admitted"]) for r in ext), 51)
        self.assertEqual({r["member_id"] for r in conflicts}, {"JPN_29", "JPN_31", "JPN_33"})


if __name__ == "__main__":
    unittest.main()
