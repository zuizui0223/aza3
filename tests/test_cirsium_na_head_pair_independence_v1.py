"""Non-source synthetic unit tests; no manufactured evolutionary observations."""
import io
import math
import pathlib
import sys
import unittest
import numpy as np
from Bio import Phylo
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"analysis"))
import audit_cirsium_na_head_pair_independence_v1 as a


def tree(): return Phylo.read(io.StringIO("((A:1,B:1):2,(C:1,D:1):2);"),"newick")
def states(): return {"A":(0,0),"B":(0,1),"C":(1,0),"D":(1,1)}


class TestModels(unittest.TestCase):
    def test_source_binary_semantics(self):
        with self.assertRaisesRegex(ValueError,"SOURCE_N63_CHANGED"):
            a.parse_source("Species,HeadPosition,OuterPhyllaryPosition\nX,0,0\n", "(X:1);")

    def test_q_generator_and_nested_parameters(self):
        ind=a.q_matrix([.1,.2,.3,.4])
        dep=a.q_matrix([.1,.2,.1,.2,.3,.4,.3,.4],True)
        np.testing.assert_allclose(ind,dep)
        np.testing.assert_allclose(ind.sum(axis=1),np.zeros(4),atol=1e-12)
        self.assertEqual(ind[0,3],0)
        self.assertEqual(ind[2,1],0)

    def test_double_change_not_single_transition(self):
        q=a.q_matrix([.1,.2,.3,.4])
        self.assertEqual(q[0,3],0)
        self.assertEqual(q[1,2],0)

    def test_fit_marginal_loglike_equals_independent_four_likelihood(self):
        t=tree(); x=states()
        head={n:s[0] for n,s in x.items()}
        outer={n:s[1] for n,s in x.items()}
        ll=(a.binary_loglik(t,head,.25,.15,0)+
             a.binary_loglik(t,outer,.3,.4,0))
        four=a.four_loglik(t,x,[.3,.4,.25,.15],False,0)
        self.assertAlmostEqual(ll,four,places=10)

    def test_conditional_simulation_correct_observed_overlap(self):
        t=tree(); x=states()
        head={n:s[0] for n,s in x.items()}
        outer={n:s[1] for n,s in x.items()}
        q=a.conditional_tip_overlap(t,head,outer,.2,.2,n_sims=20000,seed=13,min_accepted=300)
        self.assertEqual(q["observed_overlap"],1)
        self.assertEqual(q["observed_var_positives"],2)
        self.assertEqual(sum(q["overlap_histogram"]),q["accepted_exact_marginal"])
        self.assertEqual(q["condition"],"PRESERVE_TREE_AND_OTHER_TRAIT_TIP_PATTERN_AND_SIMULATED_TRAIT_POSITIVE_COUNT")
        self.assertGreater(q["p_greater_equal"],0)
        self.assertLessEqual(q["p_greater_equal"],1)

    def test_synthetic_evidence_not_a_fitness_measure(self):
        t=tree(); x=states()
        out=a.audit(t,x,15000)
        self.assertFalse(out["interpretation"].get("no_independent_genetic_module_mapping") is False)
        self.assertTrue(out["interpretation"]["no_pollinator_or_seed_predator_benefits_measured"])
        self.assertGreaterEqual(out["four_state_joint_model_comparison"]["lrt_descriptive"],0)
        self.assertTrue(out["four_state_joint_model_comparison"]["chi_square_lrt_reference_p_omitted_due_to_sparse_cells_and_boundary"])

    def test_invalid_tree_root_and_rate_shape_rejected(self):
        with self.assertRaisesRegex(ValueError,"INVALID_RATE_DIMENSIONS"):
            a.q_matrix([.1,.2],False)
        with self.assertRaisesRegex(ValueError,"INVALID_FOUR_ROOT_STATE"):
            a.four_loglik(tree(),states(),[.1]*4,root_state=5)

if __name__=="__main__":unittest.main()
