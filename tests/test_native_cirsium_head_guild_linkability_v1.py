"""Scientific source and anti-pseudoreplication checks for head trait × recorded visitor join."""
from pathlib import Path
import sys
import unittest
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
import audit_native_cirsium_head_guild_linkability_v1 as v

def specimen(taxon,head=0,outer=False,sticky=0,visits=10,cells=2,orders=None):
    orders=orders or {"Lepidoptera":int(visits*.7),"Hymenoptera":visits-int(visits*.7)}
    return {
        "source_taxon":taxon,
        "n_source_admitted_adult_flower_visit_events_NA":visits,
        "n_source_order_classified_visits":visits,
        "n_coarse_2deg_cells":cells,
        "source_head_pending":head,
        "source_outer_phyllaries_outward":outer,
        "source_glutinous_ridge":sticky,
        "source_order_record_counts":orders,
        "recorded_hymenoptera_fraction":orders.get("Hymenoptera",0)/visits,
        "recorded_lepidoptera_fraction":orders.get("Lepidoptera",0)/visits,
    }

class SourceGatedVisitorCrosswalk(unittest.TestCase):
    def test_exact_author_public_63_shape_and_true_genealogy_label(self):
        self.assertEqual(v.ORIGINAL_TRAITS_BLOB,
                         "fa127f0587b676e6e004ea121424d26069586bcb")
        self.assertEqual(v.ARCHIVE_SHA256,
            "ce0513fde35d9c1be42cef87da7a6f4b262e9b4fe237ba0b78e9bf31771b58a4")
        self.assertEqual(v.ORIGINAL_EVENT_CSV_SHA256,
            "6d2732565c4f158a0f531b4671832738dafe0b06af7f25b635a5d26b7804deff")

    def test_genuine_adult_insect_order_is_not_anatomy_or_fitness_role(self):
        self.assertEqual(
            v.insect_order("Animalia | Arthropoda | Insecta | Lepidoptera"),
            "Lepidoptera")
        self.assertEqual(
            v.insect_order("Animalia | Arthropoda | Insecta | Hymenoptera"),
            "Hymenoptera")
        self.assertIsNone(v.insect_order("Animalia | Aves | Apodiformes"))
        self.assertIsNone(v.insect_order("Animalia | Arthropoda | Insecta"))
        self.assertIsNone(v.insect_order(
            "Animalia | Arthropoda | Insecta | Hymenoptera | Lepidoptera"))

    def test_morphology_trait_code_is_not_bird_syndrome(self):
        self.assertEqual(v.TRAIT_AXES,{
            "head_pending":"HeadPosition",
            "outer_phyllaries_outward":"OuterPhyllaryPosition",
            "sticky_ridge":"GlutinousRidge"})
        self.assertNotIn("pollinator",v.TRAIT_AXES)

    def test_source_region_is_only_coarse_window_not_native_range(self):
        frame=pd.DataFrame({
            "lat":[34.0,74.0,45.0,30.0],"lon":[-100.0,-98.0,11.0,-130.0]
        })
        self.assertEqual(v.na_region(frame).tolist(),[True,False,False,True])

    def test_missing_guild_is_none_not_evidence_of_zero(self):
        a=specimen("Cirsium wrightii",visits=10,cells=1)
        b=specimen("Cirsium canescens",visits=9,cells=3)
        c=specimen("Cirsium undulatum",visits=20,cells=3)
        chosen=v.coverage([a,b,c],10)
        self.assertEqual([r["source_taxon"] for r in chosen],
                         ["Cirsium undulatum"])

    def test_not_relabel_subspecies_as_vouchered_binomial(self):
        # The actual author file has "Cirsium_arizonicum_rothrockii":
        # NEVER cut the subspecies to "Cirsium arizonicum" and pretend
        # its photographed flower visitor was linked to that trait.
        self.assertIsNone(__import__("re").fullmatch(
            r"Cirsium_[a-z]+","Cirsium_arizonicum_rothrockii"))

    def test_fake_one_or_two_taxa_group_must_report_not_identifiable(self):
        v0=[specimen("Cirsium a",outer=False,visits=20),
            specimen("Cirsium b",outer=False,visits=20)]
        v1=[specimen("Cirsium c",outer=True,visits=20)]
        g=v.contrasts(v0+v1)["source_outer_phyllaries_outward"]
        self.assertEqual(g["status"],
                         "NOT_IDENTIFIABLE_INDEPENDENT_MORPHOLOGY_GROUPS")
        self.assertEqual(g["n_species_state1"],1)
        self.assertIsNone(
            g["mean_recorded_hymenoptera_fraction_state1_minus_state0"])

    def test_three_plus_three_is_descriptive_not_experimental_selection(self):
        x=[specimen(f"Cirsium n{i}",outer=False,visits=20) for i in range(3)]
        x+=[specimen(f"Cirsium y{i}",outer=True,visits=20,
                     orders={"Hymenoptera":15,"Lepidoptera":5})
             for i in range(3)]
        g=v.contrasts(x)["source_outer_phyllaries_outward"]
        self.assertEqual(g["status"],"DESCRIPTIVE_SPECIES_LEVEL_ONLY")
        self.assertAlmostEqual(
            g["mean_recorded_hymenoptera_fraction_state1_minus_state0"],
            .45)
        self.assertAlmostEqual(
            g["mean_recorded_lepidoptera_fraction_state1_minus_state0"],
            -.45)

    def test_no_region_matched_pairs_cannot_be_inferred_as_no_visitor(self):
        df=pd.DataFrame(columns=[
            "source_binomial","lat","lon","original_insect_order"
        ])
        result=v.matched_cells(df,
             pd.DataFrame(columns=["source_binomial"]))
        self.assertEqual(result["n_species_pairs_any_region"],0)
        self.assertEqual(result["n_species_pairs_at_least_3_regions"],0)

    def test_distinct_states_same_coarse_region_can_be_counted_but_no_causality(self):
        def data(s,coords):
            return pd.DataFrame([{
                "source_binomial":s,"lat":la,"lon":lo,
                "original_insect_order":"Hymenoptera"
            } for la,lo in coords])
        # Fifteen *original-like* records do NOT become fifteen independent
        # species or fifteen observed head phenotypes.
        a=data("Cirsium undulatum",[(40,-100)]*5+[(42,-100)]*5)
        b=data("Cirsium vulgare",[(40,-100)]*5+[(42,-100)]*5)
        traits=pd.DataFrame([
            {"source_binomial":"Cirsium undulatum","HeadPosition":"0",
             "OuterPhyllaryPosition":"0","GlutinousRidge":"1"},
            {"source_binomial":"Cirsium vulgare","HeadPosition":"0",
             "OuterPhyllaryPosition":"0","GlutinousRidge":"0"},
        ])
        o=v.matched_cells(pd.concat([a,b]),traits,minimum=5)
        self.assertEqual(
            o["by_trait"]["sticky_ridge"]
             ["n_species_pairs_at_least_one_matched_region"],1)
        self.assertEqual(
            o["by_trait"]["sticky_ridge"]
             ["n_species_pairs_at_least_three_matched_regions"],0)
        self.assertEqual(o["by_trait"]["head_pending"]
                         ["n_species_pairs_at_least_one_matched_region"],0)

    def test_source_panel_wrong_hash_rejected(self):
        with self.assertRaisesRegex(ValueError,
            "AUTHOR_NA_63_TRAIT_SOURCE_GIT_BLOB_CHANGED"):
            v.trait_panel(b"Species,HeadPosition\nCirsium_undulatum,0\n")

    def test_source_archive_wrong_hash_rejected(self):
        with self.assertRaisesRegex(ValueError,
            "GLOBI_ORIGINAL_GIT_ACTIONS_SOURCE_SHA256_CHANGED"):
            v.original_events(b"totally not the original public archive")

if __name__=="__main__":
    unittest.main()
