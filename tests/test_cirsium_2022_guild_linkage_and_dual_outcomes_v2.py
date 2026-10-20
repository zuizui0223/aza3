#!/usr/bin/env python3
"""Synthetic tests of source linkage and dual-stage inference; NO synthetic field result."""
import pathlib
import sys
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "analysis"))
import audit_cirsium_2022_four_arm_head_outcomes_v1 as v1
import audit_cirsium_2022_guild_linkage_and_dual_outcomes_v2 as v2


def seed(site, block, arm, viable, emerged, flowers=100):
    return {
        "Site_Code": site, "Block_Code": block, "Treatment_Code": arm,
        "Vial_Code": "1", "Flowers": str(flowers),
        "Viable_Seeds": str(viable), "NViable_Seeds": str(flowers-viable),
        "Damaged_Seeds": "0", "Missing_Seeds": "0", "Insects": str(emerged)
    }


def insect(site, block, arm, guild, number):
    return {
        "Site_Code": site, "Block_Code": block, "Treatment_Code": arm,
        "Vial_Code": "1", "Feeding_Guild": guild, "Number": str(number),
        "Order": "Hymenoptera", "Sub_Order": "NA", "Family": "NA",
        "Sp_Code": "NA", "Stage": "adult", "size": "1"
    }


class LinkageTest(unittest.TestCase):
    def test_normalizes_only_casing_known_to_differ(self):
        self.assertEqual(v2.plant_key(seed("ditch", "B4", "Control", 5, 0)),
                         v2.plant_key(insect("Ditch", "b4", "Control", "Herbivore", 1)))

    def test_no_identity_claim_from_matching_aggregates(self):
        ss = [seed("ditch", "B4", "Control", 60, 2)]
        hh = [insect("Ditch", "B4", "Control", "Parasitoid", 2)]
        observed = v2.linkage_summary(ss, hh, strict=False)
        self.assertEqual(observed["n_matching_aggregate_counts"], 1)
        self.assertFalse(observed["confirmed_head_or_individual_insect_identity"])
        self.assertFalse(observed["confirmed_seed_feeder_suppression"])

    def test_zero_in_herbivore_file_is_not_automatic_confirmed_absence(self):
        ss = [seed("ditch", "B4", "Control", 60, 3)]
        r = v2.linkage_summary(ss, [], strict=False)
        self.assertEqual(r["n_groups_with_no_identification_rows"], 1)
        self.assertEqual(r["n_groups_no_identification_but_positive_seed_count"], 1)
        self.assertEqual(r["n_mismatching_aggregate_counts"], 1)

    def test_nonmatching_counts_do_not_join_at_head_level(self):
        ss = [seed("ditch", "B4", "Control", 60, 3)]
        hh = [insect("Ditch", "B4", "Control", "Herbivore", 1)]
        r = v2.linkage_summary(ss, hh, strict=False)
        self.assertEqual(r["n_mismatching_aggregate_counts"], 1)
        self.assertEqual(r["n_identified_lower_than_source_insects"], 1)

    def test_counterpart_guilds_not_merged(self):
        ss = [seed("ditch", "B4", "Control", 60, 3)]
        hh = [insect("Ditch", "B4", "Control", "Herbivore", 2),
              insect("Ditch", "B4", "Control", "Parasitoid", 1)]
        r = v2.linkage_summary(ss, hh, strict=False)
        self.assertEqual(r["identified_guild_count"], {"Herbivore": 2, "Parasitoid": 1})

    def test_invalid_insect_counts_rejected(self):
        data = ("Site_Code,Block_Code,Treatment_Code,Vial_Code,Feeding_Guild,"
                "Order,Sub_Order,Family,Sp_Code,Number,Stage,size\n"
                "Ditch,B4,Control,1,Herbivore,Diptera,NA,NA,NA,-1,adult,1\n").encode()
        with self.assertRaisesRegex(ValueError, "INVALID_INSECT_NUMBER"):
            v2.parse_identification_rows(data, strict=False)

    def test_unsupported_guild_rejected(self):
        data = ("Site_Code,Block_Code,Treatment_Code,Vial_Code,Feeding_Guild,"
                "Order,Sub_Order,Family,Sp_Code,Number,Stage,size\n"
                "Ditch,B4,Control,1,Pollinator,Diptera,NA,NA,NA,1,adult,1\n").encode()
        with self.assertRaisesRegex(ValueError, "UNEXPECTED_GUILD"):
            v2.parse_identification_rows(data, strict=False)

    def test_positive_signed_eaten_filter_is_outcome_selection(self):
        a=seed("ditch","B4","Control",80,2,100)
        b=seed("ditch","B4","Control",80,1,100)
        c=seed("ditch","B4","Pollinator Exclusion",10,0,100)
        a["Missing_Seeds"]="5"
        a["NViable_Seeds"]="15"
        b["Missing_Seeds"]="-2"
        b["NViable_Seeds"]="22"
        c["Missing_Seeds"]="0"
        result=v2.signed_eaten_filter_audit([a,b,c],strict=False)
        self.assertEqual(result["n_eligible_heads"],3)
        self.assertEqual(result["n_eaten_positive_heads"],1)
        self.assertEqual(result["n_excluded_due_to_eaten_not_positive"],2)
        self.assertTrue(result["post_outcome_selection"])
        self.assertFalse(result["valid_unbiased_florivore_inference"])

    def test_filter_never_calls_negative_residual_a_zero_insect(self):
        a=seed("ditch","B4","Control",80,5,100)
        a["Missing_Seeds"]="-2"
        a["NViable_Seeds"]="22"
        q=v2.signed_eaten_filter_audit([a],strict=False)
        self.assertEqual(q["n_eaten_positive_heads"],0)
        self.assertEqual(q["n_excluded_due_to_eaten_not_positive"],1)
        self.assertEqual(q["by_treatment"]["Control"]["insects_per_head_all_eligible"],5)

    def test_insects_are_head_mean_not_floret_fraction(self):
        rows = [seed("ditch", "B4", "Control", 50, 2, 100),
                seed("ditch", "B4", "Control", 5, 4, 10)]
        estimates = v1.group_estimate(rows)
        self.assertAlmostEqual(estimates["emerged_insects_per_head"], 3)
        self.assertAlmostEqual(estimates["viable_fraction"], 55/110)

    def test_missing_not_zero_and_fraction_not_fitness(self):
        ss = [seed("ditch", "B4", "Control", 60, 3)]
        summary = v2.linkage_summary(ss, [], strict=False)
        self.assertTrue(summary["claim_limits"]["no_herbivore_row_is_not_undamaged_head"])
        self.assertFalse(summary["claim_limits"]["identifies_physical_spine_barrier"])
        self.assertFalse(summary["claim_limits"]["identifies_pollination_resource_path"])

    def test_four_complete_blocks_do_not_pseudoreplicate_heads(self):
        rows=[]
        for i,site in enumerate(v1.SITES):
            for j,arm in enumerate(v1.ARMS):
                rows.append(seed(site,"B"+str(i),arm,70-j*12,5-j))
                rows.append(seed(site,"B"+str(i),arm,70-j*12,5-j))
        _, proxy, blocks, complete=v1.prep(rows, strict=False)
        self.assertEqual(len(complete),4)
        r=v2.dual_stage_summary(complete, reps=102, seed=20261020)
        self.assertEqual(r["n_four_arm_blocks"],4)
        self.assertEqual(r["n_source_heads_treated_as_independent_replicates"],0)
        self.assertEqual(r["bootstrap_draws_conditional_on_fixed_sites"],102)
        self.assertLess(r["contrasts"]["Pollinator Exclusion MINUS Control"]["viable_fraction"]["mean"],0)

    def test_same_block_nested_model_comparison_identity(self):
        rows=[]
        for i,site in enumerate(v1.SITES):
            for j,arm in enumerate(v1.ARMS):
                rows.append(seed(site,"B"+str(i),arm,70+j,6-j))
        _, _, _, complete=v1.prep(rows, strict=False)
        r=v2.dual_stage_summary(complete,reps=102)
        for m in ("viable_fraction","emerged_insects_per_head"):
            c=r["contrasts"]
            self.assertAlmostEqual(
                c["Predator Model MINUS Control"][m]["mean"],
                c["Predator Exclusion MINUS Control"][m]["mean"]+
                c["Predator Model MINUS Predator Exclusion"][m]["mean"]
            )


if __name__=="__main__":
    unittest.main()
