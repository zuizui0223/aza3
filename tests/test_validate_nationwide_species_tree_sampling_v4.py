#!/usr/bin/env python3
import importlib.util, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("v4", ROOT / "analysis/validate_nationwide_species_tree_sampling_v4.py")
mod = importlib.util.module_from_spec(SPEC); SPEC.loader.exec_module(mod)

class NationwideSpeciesTreeV4Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls): cls.d = mod.validate()
    def test_primary_product_is_nationwide_tree(self):
        self.assertEqual(self.d["primary_product"], "JAPAN_WIDE_NUCLEAR_SPECIES_TREE_AND_NETWORK")
    def test_planning_numbers(self):
        self.assertEqual(self.d["taxonomic_census"]["planning_species_baseline"], 120)
        self.assertEqual(self.d["nationwide_sampling"]["nominal_full_bank_at_120_species"], 290)
        self.assertEqual(self.d["nationwide_sampling"]["full_bank_envelope_if_125_species"], 300)
    def test_target_capture_precedes_rad(self):
        self.assertEqual(self.d["target_capture"]["role"], "PRIMARY_CROSS_SPECIES_SCAFFOLD")
        self.assertEqual(self.d["population_rad"]["not_role"], "NOT_THE_NATIONWIDE_SPECIES_TREE")
    def test_nested_rad_numbers(self):
        self.assertEqual(self.d["population_rad"]["target_focal_primary_rad_total"], 168)
        self.assertEqual(self.d["population_rad"]["nominal_total_unique_physical_plants_if_full_v4_is_completed"], 438)
    def test_fail_closed(self):
        self.assertEqual(self.d["current_state"]["physical_samples"], 0)
        self.assertFalse(self.d["current_state"]["nationwide_collection_authorized"])

if __name__ == "__main__": unittest.main()
