"""Tests for exact source-backed Cirsium trait vs pollination association.

Synthetic trees only test algorithms, NEVER evidence of adaptive convergence.
"""
import io
import math
import pathlib
import sys
import unittest
import numpy as np
from Bio import Phylo
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"analysis"))
import audit_cirsium_na_head_pollinator_phylo_v1 as v


class TestOriginalSourceSemantics(unittest.TestCase):
    def test_fisher_exact_source_contingency(self):
        t=v.fisher_binary(6,1,10,46)
        self.assertAlmostEqual(t["odds_ratio"],27.6)
        self.assertAlmostEqual(t["p_two_sided"],0.0007009516685550099,places=12)
        self.assertAlmostEqual(t["bird_minus_insect_prevalence"],6/7-10/56)

    def test_no_mechanical_significance_if_trait_constant(self):
        t=v.fisher_binary(0,7,0,56)
        self.assertEqual(t["p_two_sided"],1.0)
        self.assertIsNone(t["odds_ratio"])

    def test_unique_source_taxa_required(self):
        a=[{"Species":"Cirsium_A","OuterPhyllaryPosition":"1"},
           {"Species":"Cirsium_A","OuterPhyllaryPosition":"0"}]
        with self.assertRaisesRegex(ValueError,"DUPLICATE_TAXON"):
            v.validate_unique_taxa(a,"Species")

    def test_prohibited_trait_code_fails(self):
        with self.assertRaisesRegex(ValueError,"UNKNOWN_HEAD_TRAIT_CODE"):
            v.trait_positive({"HeadPosition":"3"},"head_pendulous")

    def test_bh_multiple_test_adjustment(self):
        x=v.bh_adjust([0.01,0.03,0.2,0.0005])
        self.assertEqual(len(x),4)
        self.assertAlmostEqual(x[3],0.002)
        self.assertAlmostEqual(x[0],0.02)
        self.assertAlmostEqual(x[1],0.04)
        self.assertAlmostEqual(x[2],0.2)

    def test_source_matches_tree_exact_tip_names(self):
        tree=Phylo.read(io.StringIO("((A:1,B:1):1,C:2);"),"newick")
        v.verify_tip_set(tree,{"A","B","C"})
        with self.assertRaisesRegex(ValueError,"TREE_TIP_SET_MISMATCH"):
            v.verify_tip_set(tree,{"A","C"})

    def test_binary_parsimony_counts_changes_not_gain_events(self):
        t=Phylo.read(io.StringIO("((A:1,B:1):1,(C:1,D:1):1);"),"newick")
        st={"A":0,"B":1,"C":0,"D":1}
        self.assertEqual(v.parsimony_cost(t,st,root_state=0),2)

    def test_Mk_root_likelihood_finite(self):
        t=Phylo.read(io.StringIO("((A:1,B:1):1,C:2);"),"newick")
        z=v.binary_mk_loglik(t,{"A":0,"B":1,"C":0},a=0.2,b=0.5)
        self.assertTrue(math.isfinite(z))
        self.assertLess(z,0)

    def test_conditional_simulation_sees_exact_target_count(self):
        t=Phylo.read(io.StringIO("((A:1,B:1):1,(C:1,D:1):1);"),"newick")
        q=v.simulate_mk_conditional(t,{"A":0,"B":1,"C":0,"D":1},
                                     {"A":[1,0],"B":[0,1],"C":[1,0],"D":[0,1]},
                                     a=0.2,b=0.2,n_sims=10000,seed=10,min_accepted=100)
        self.assertGreater(q["n_conditioned_draws"],100)
        self.assertEqual(q["observed_bird_tips"],2)
        self.assertTrue(all(0<=d["upper_tail_probability"]<=1 for d in q["phenotypes"].values()))
        self.assertEqual(q["null_condition"],"SIMULATED_BINARY_CTMC_CONDITION_ON_EXACT_OBSERVED_BIRD_TIP_COUNT")

    def test_equal_rate_Mk_fit_is_separate_sensitivity(self):
        t=Phylo.read(io.StringIO("((A:1,B:1):1,(C:1,D:1):1);"),"newick")
        model=v.fit_equal_rate_binary_mk(t,{"A":0,"B":1,"C":0,"D":1})
        self.assertTrue(model["fit_converged"])
        self.assertGreater(model["equal_rate"],0)
        self.assertTrue(math.isfinite(model["maximum_loglik"]))

    def test_nonred_head_association_is_a_separate_small_sample(self):
        rows=[
            {"Species":"A","CorollaColor":"0","OuterPhyllaryPosition":"1"},
            {"Species":"B","CorollaColor":"1","OuterPhyllaryPosition":"1"},
            {"Species":"C","CorollaColor":"2","OuterPhyllaryPosition":"0"},
            {"Species":"D","CorollaColor":"2","OuterPhyllaryPosition":"1"},
            {"Species":"E","CorollaColor":"2","OuterPhyllaryPosition":"0"},
        ]
        poll={"A":"Hummingbird","B":"Hummingbird","C":"Hummingbird",
              "D":"Insect","E":"Insect"}
        r=v.nonred_phenotype_associations(rows,poll,
                                            ["outer_phyllaries_spreading_or_reflexed"])
        self.assertEqual(r["n_bird_nonred"],2)
        self.assertEqual(r["n_insect_nonred"],2)
        self.assertEqual(r["outcomes"]["outer_phyllaries_spreading_or_reflexed"]["bird_positive"],1)

    def test_head_module_mosaic_not_forced_to_be_synchronous(self):
        rows=[{"Species":"A","OuterPhyllaryPosition":"1","InnerPhyllaryColor":"0","HeadPosition":"0"},
              {"Species":"B","OuterPhyllaryPosition":"0","InnerPhyllaryColor":"1","HeadPosition":"1"},
              {"Species":"C","OuterPhyllaryPosition":"1","InnerPhyllaryColor":"1","HeadPosition":"1"},
              {"Species":"D","OuterPhyllaryPosition":"0","InnerPhyllaryColor":"0","HeadPosition":"0"}]
        p={"A":"Hummingbird","B":"Hummingbird","C":"Insect","D":"Insect"}
        z=v.head_state_combinations(rows,p)
        self.assertEqual(z["Hummingbird"]["outer_only"],1)
        self.assertEqual(z["Hummingbird"]["inner_only"],1)
        self.assertEqual(z["Insect"]["both_outer_and_inner"],1)

    def test_result_not_an_adaptive_mechanism(self):
        d=v.claim_limits()
        for key in ("same_head_morphology_and_fitness_measured",
                    "direct_bird_pollination_assayed","adaptation_or_convergent_selection_identified",
                    "independent_origins_of_correlated_trait_changes_identified"):
            self.assertIs(d[key],False)

if __name__=="__main__":
    unittest.main()
