"""Fail-closed semantics for manually source-anchored nonexclusive Cirsium visitors."""
import copy
import json
import pathlib
import sys
import unittest

sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"analysis"))
import audit_cirsium_na_independent_visitor_labels_v1 as v

HERE=pathlib.Path(__file__).resolve().parents[1]
CATALOG=HERE/"data/evidence/cirsium_na_independent_visitor_source_audit_v1.json"


class TestOriginalLabelAndObservationSeparation(unittest.TestCase):
    def setUp(self):
        self.audit=json.loads(CATALOG.read_text(encoding="utf-8"))
        self.traits={z["taxon"]:{"Species":z["taxon"], **z["author_trait_codes"]}
                     for z in self.audit["cases"]}
        self.labels={z["taxon"]:z["author_syndrome"] for z in self.audit["cases"]}

    def test_original_classification_and_independent_visits_not_equivalent(self):
        r=v.audit_four(self.audit,self.traits,self.labels,expect_full=False)
        self.assertEqual(r["n_taxa_spotchecked"],4)
        self.assertEqual(r["n_insect_coded_with_external_hummingbird_visits"],1)
        self.assertEqual(r["n_bird_coded_with_external_insect_visits"],3)
        self.assertFalse(r["direct_effective_pollination_confirmed"])

    def test_wrightii_remains_insect_syndrome_after_actual_bird_visitor(self):
        x=next(z for z in self.audit["cases"] if z["taxon"]=="Cirsium_wrightii")
        self.assertEqual(x["author_syndrome"],"Insect")
        self.assertTrue(x["visitor_bird_observed"])
        self.assertEqual(x["reproductive_service_by_bird"],"UNKNOWN")
        self.assertFalse(v.claim_limits()["bird_visit_means_bird_effective_pollination"])

    def test_rejects_relabeling_insect_as_bird_just_from_visit(self):
        a=copy.deepcopy(self.audit)
        a["cases"][0]["author_syndrome"]="Hummingbird"
        with self.assertRaisesRegex(ValueError,"ORIGINAL_SOURCES_DISAGREE"):
            v.audit_four(a,self.traits,self.labels,expect_full=False)

    def test_invalid_ext_taxon_or_misassigned_morphology_fails(self):
        a=copy.deepcopy(self.audit)
        a["cases"][0]["author_trait_codes"]["HeadPosition"]="1"
        with self.assertRaisesRegex(ValueError,"ORIGINAL_SOURCES_DISAGREE"):
            v.audit_four(a,self.traits,self.labels,expect_full=False)

    def test_rejects_unknown_effectiveness_promoted_to_confirmed(self):
        a=copy.deepcopy(self.audit)
        a["cases"][0]["reproductive_service_by_bird"]="CONFIRMED"
        with self.assertRaisesRegex(ValueError,"UNSUPPORTED_POLLINATION_EFFECTIVENESS"):
            v.audit_four(a,self.traits,self.labels,expect_full=False)

    def test_rejects_invented_direct_raw_visitor_effort(self):
        a=copy.deepcopy(self.audit)
        a["summary"]["independent_bird_pollination_efficiency_experiments_in_curated_records"]=1
        with self.assertRaisesRegex(ValueError,"NO_INDEPENDENT_EFFECTIVENESS_EXPERIMENT"):
            v.audit_four(a,self.traits,self.labels,expect_full=False)

    def test_fna_source_required_for_every_curated_taxon(self):
        a=copy.deepcopy(self.audit)
        a["cases"][1]["external_urls"]=["https://some-blog.example/"]
        with self.assertRaisesRegex(ValueError,"MISSING_PRIMARY_FLORA_REFERENCE"):
            v.audit_four(a,self.traits,self.labels,expect_full=False)

    def test_no_repeated_taxon_as_fake_independent_record(self):
        a=copy.deepcopy(self.audit)
        a["cases"][1]["taxon"]=a["cases"][0]["taxon"]
        with self.assertRaisesRegex(ValueError,"DUPLICATE_TAXON"):
            v.audit_four(a,self.traits,self.labels,expect_full=False)

    def test_original_species_complex_broad_visitors_not_as_two_subspecies(self):
        a=copy.deepcopy(self.audit)
        a["broader_complex_warning"]["may_not_assign_to_both_original_subspecies_as_independent_field_visits"]=False
        with self.assertRaisesRegex(ValueError,"ARIZONICUM_AGGREGATE_IS_NOT_SUBSPECIES_REPLICATION"):
            v.audit_four(a,self.traits,self.labels,expect_full=False)

    def test_curated_source_cannot_claim_measured_selective_advantage(self):
        x=v.claim_limits()
        self.assertFalse(x["actual_same_head_visit_pollen_and_fitness"])
        self.assertFalse(x["head_architecture_as_mechanical_bird_access_barrier_identified"])
        self.assertFalse(x["evolutionary_gain_or_reexpression_identified"])

if __name__=="__main__":
    unittest.main()
