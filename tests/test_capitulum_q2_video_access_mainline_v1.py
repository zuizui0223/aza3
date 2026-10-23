"""Protect aza3 Q2 field-intake inference when there are NO real observations.

Source test functions originate in the audited but still-unmerged aza3 PR #44
(capitulm functional access), ported without importing its stacked branches.
All generated episodes here are SYNTHETIC and not included in intake files.
"""
from pathlib import Path
import csv
import json
import sys
import unittest
import copy

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import validate_capitulum_video_effort_denominator_v1 as e
import validate_capitulum_guild_access_event_ledger_v1 as v

EFFORT_CONTRACT=ROOT/"data/contracts/capitulum_video_effort_denominator_v1.json"
EVENT_CONTRACT=ROOT/"data/contracts/capitulum_guild_access_event_contract_v1.json"
EFFORT_CSV=ROOT/"data/intake/capitulum_video_effort_denominator_v1.csv"
EVENT_CSV=ROOT/"data/intake/capitulum_guild_access_event_ledger_v1.csv"

class ActualAccessIntakeReadiness(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.es=json.loads(EFFORT_CONTRACT.read_text(encoding="utf-8"))
        cls.vs=json.loads(EVENT_CONTRACT.read_text(encoding="utf-8"))

    def test_exact_mainline_schema_matches_empty_real_intake(self):
        h,rows=e.read_csv(EFFORT_CSV)
        eh,events=e.read_csv(EVENT_CSV)
        self.assertEqual(h,self.es["columns"])
        self.assertEqual(eh,self.vs["required_fields"])
        self.assertEqual(rows,[])
        self.assertEqual(events,[])
        self.assertEqual(len(h),27)
        self.assertEqual(len(eh),32)

    def test_no_real_video_returns_not_identifiable_not_zero_visitor_rate(self):
        z=e.validate([],[],self.es)
        self.assertEqual(z["status"],"NOT_IDENTIFIABLE_NO_REAL_VIDEO_BOUTS")
        self.assertIsNone(z["observed_raw_arrival_rate_per_minute"])
        self.assertEqual(z["n_effort_bouts"],0)

    def test_synthetic_video_denominator_positive_and_negative_checks(self):
        r=e.tests(self.es)
        self.assertEqual(r["status"],"PASS_SYNTHETIC_DENOMINATOR_FIREWALL")
        self.assertGreaterEqual(r["negative_cases"],13)
        self.assertTrue(r["no_real_organismal_field_data_submitted"])

    def test_synthetic_event_ontology_checks_no_false_pollen_or_egg_success(self):
        r=v.test_fixture(self.vs)
        self.assertEqual(r["status"],"PASS_SYNTHETIC_SEMANTIC_FIREWALL_NOT_FIELD_DATA")
        self.assertEqual(r["n_expected_failure_checks"],8)

    def test_verified_zero_is_distinct_from_missing_camera_coverage(self):
        raw=e.make_effort()
        ok=e.validate([copy.deepcopy(raw)],[],self.es)
        self.assertEqual(ok["n_verified_zero_approach_bouts"],1)
        self.assertEqual(ok["observed_raw_arrival_rate_per_minute"],0)
        bad={**raw,"video_review_mode":"event_triggered_only"}
        with self.assertRaisesRegex(ValueError,"FALSE_VALID_RATE_DENOMINATOR"):
            e.validate([bad],[],self.es)

    def test_unmeasured_other_cirsium_hosts_cannot_be_fabricated_zero(self):
        raw=e.make_effort()
        raw["other_cirsium_buds_count"]="0"
        with self.assertRaisesRegex(ValueError,"UNMEASURED_OTHER_CIRSIUM"):
            e.validate([raw],[],self.es)

    def test_visitor_event_must_be_linked_to_real_camera_bout(self):
        evt=e.make_event()
        with self.assertRaisesRegex(ValueError,"UNMATCHED_EFFORT_BOUT"):
            e.validate([], [evt], self.es)
        good=e.make_effort()
        good["zero_approaches_confirmed"]="0"
        result=e.validate([good],[evt],self.es)
        self.assertEqual(result["eligible_approach_episode_count"],1)
        self.assertEqual(result["observed_raw_arrival_rate_per_minute"],0.5)

    def test_pre_entry_guild_not_backfilled_from_success(self):
        c=self.vs
        valid={name:"" for name in c["required_fields"]}
        valid.update({
            "record_id":"SYNTHETIC_1",
            "individual_id":"TEST_PLANT",
            "population_id":"TEST_POP",
            "capitulum_id":"TEST_HEAD",
            "observation_bout_id":"TEST_BOUT",
            "approach_episode_id":"TEST_E1",
            "video_clip_id":"synthetic/clip",
            "time_from_bout_start_s":"1",
            "phenological_stage":"full_anthesis",
            "reproductive_sex_state":"hermaphroditic",
            "visitor_taxon_label":"undetermined insect",
            "pre_entry_guild":"seed_feeder_candidate",
            "pre_entry_guild_evidence":"unresolved",
            "behavioral_role":"unknown",
            "role_evidence":"not_assessed",
            "entry_route_world":"undetermined",
            "entry_route_head":"undetermined",
            "approached":"1",
            "landed_or_hover_access":"0",
            "reproductive_zone_reached":"0",
            "anther_stigma_contact_observed":"0",
            "pollen_deposition_assay":"not_assessed",
            "oviposition_confirmed":"NA",
            "predator_or_parasitoid_search_observed":"NA",
            "host_attacked_confirmed":"NA",
            "host_killed_before_seed_damage":"NA",
            "seed_damage_timing_assessed":"NA",
            "seeds_saved_inferred":"NA",
            "evidence_uri":"synthetic/clip",
            "scorer_blinded_to_head_traits":"1",
        })
        field=["individual_id","population_id","capitulum_id",
              *c["head_traits_from_existing_field_ledger"]]
        bout=["individual_id","population_id","capitulum_id","observation_bout_id"]
        plant=["individual_id","population_id","phenology_census_id"]
        with self.assertRaisesRegex(ValueError,"POST_OUTCOME_ROLE_CANNOT_DEFINE_ARRIVAL_GUILD"):
            v.validate([valid],c,field,bout,plant)

    def test_no_pretend_head_fitness_in_video_only_data(self):
        self.assertIn("No field intervention authorization",self.vs["authority"])
        self.assertIn("No field",self.vs["authority"])
        self.assertIn("NOT_IDENTIFIABLE_NO_REAL_VIDEO_BOUTS",e.validate([],[],self.es)["status"])
        r=v.validate([],self.vs,
            ["individual_id","population_id","capitulum_id",
              *self.vs["head_traits_from_existing_field_ledger"]],
            ["individual_id","population_id","capitulum_id","observation_bout_id"],
            ["individual_id","population_id","phenology_census_id"])
        self.assertEqual(r["n_approach_episodes"],0)
        self.assertEqual(r["direct_trait_impact_on_seed_fitness"],"NOT_IDENTIFIED")

if __name__=="__main__":
    unittest.main()
