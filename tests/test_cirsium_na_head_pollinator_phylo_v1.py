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

    def test_result_not_an_adaptive_mechanism(self):
        d=v.claim_limits()
        for key in ("same_head_morphology_and_fitness_measured",
                    "direct_bird_pollination_assayed","adaptation_or_convergent_selection_identified",
                    "independent_origins_of_correlated_trait_changes_identified"):
            self.assertIs(d[key],False)

if __name__=="__main__":
    unittest.main()
