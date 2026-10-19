"""Negative controls for Cirsium protein-family topology (not functional orthology)."""
from io import StringIO
from pathlib import Path
import shutil
import sys
import unittest

from Bio import Phylo

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis/anto_empirical"))
from cirsium_flavonoid_gene_tree_v5 import (
    CASES, load_sources, load_reference_proteins,
    chosen_source_peptide,source_alignment_input,run_mafft,
    _pairwise_p_distance,distance_matrix,unrooted_family_split,
    bootstrap_family_support,run_panel
)

DATA=ROOT/"data/anto_source_handoff/molecular"


class FamilyTopologyUnitTests(unittest.TestCase):
    def setUp(self):
        self.pos={"D1","D2","D3"}
        self.neg={"A1","A2"}

    def _classify(self,tree):
        return unrooted_family_split(Phylo.read(StringIO(tree),"newick"),
                                     "Q",self.pos,self.neg)

    def test_positive_side_requires_group_separator(self):
        self.assertEqual(self._classify("((((D1,D2),Q),D3),(A1,A2));"),
                         "POSITIVE_FAMILY_SIDE")

    def test_related_enzyme_side_is_not_forced_into_dfr(self):
        self.assertEqual(self._classify("((D1,D2,D3),((A1,Q),A2));"),
                         "RELATED_FAMILY_SIDE")

    def test_middle_of_families_is_not_a_certified_gene(self):
        self.assertEqual(self._classify("((D1,D2,D3),Q,(A1,A2));"),
                         "BETWEEN_FAMILY_ANCHORS")

    def test_scrambled_reference_anchors_are_not_classifiable(self):
        self.assertEqual(self._classify("(((D1,A1),D2),((D3,A2),Q));"),
                         "REFERENCE_FAMILIES_NOT_SEPARATED")

    def test_differing_branch_lengths_do_not_change_split(self):
        self.assertEqual(self._classify(
            "((((D1:1,D2:2):1,Q:15):1,D3:1):1,(A1:1,A2:1):1);"),
            "POSITIVE_FAMILY_SIDE")

    def test_degenerate_one_family_rejected(self):
        t=Phylo.read(StringIO("((((D1,D2),Q),D3),(A1,A2));"),"newick")
        with self.assertRaises(ValueError):
            unrooted_family_split(t,"Q",set(),self.neg)
        with self.assertRaises(ValueError):
            unrooted_family_split(t,"D1",self.pos,self.neg)

    def test_gaps_and_unknown_residues_are_excluded_from_identity(self):
        aa1="A"*70+"-"*20
        aa2="A"*60+"B"*10+"-"*20
        d,n=_pairwise_p_distance(aa1,aa2)
        self.assertEqual(n,70)
        self.assertAlmostEqual(d,10/70)

    def test_short_overlap_cannot_create_false_clade(self):
        with self.assertRaises(ValueError):
            _pairwise_p_distance("A"*40+"-"*70,"A"*110)
        with self.assertRaises(ValueError):
            distance_matrix({"A":"A"*70,"B":"T"*65})


class RealPublishedMaterialTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.src,_=load_sources(
            DATA/"cirsium_japonicum_2020_26_source_sequences_v6.fasta",
            DATA/"cirsium_japonicum_2020_26_sequence_provenance_v6.json")
        cls.refs=load_reference_proteins(
            DATA/"DFR_ANR_ANS_FLS_reference_proteins_v5.fasta",
            DATA/"DFR_ANR_ANS_FLS_reference_provenance_v5.json")

    def test_eleven_are_not_mixed_into_a_single_tree(self):
        self.assertEqual(len(self.src),26)
        self.assertEqual(len(self.refs),9)
        self.assertEqual(len(CASES),3)
        for case,v in CASES.items():
            groups={r["group"] for r in self.refs.values()
                    if r["group"] in {v["positive"],v["related"]}}
            self.assertEqual(groups,{v["positive"],v["related"]})

    def test_primary_and_alt_2020_isoforms_are_not_same_tree_independent_replicates(self):
        a,b=(set(CASES[key]["source"]) for key in
             ("DFR_ANR_primary","DFR_ANR_alt_isoform"))
        self.assertEqual(len(a),3)
        self.assertEqual(len(b),3)
        self.assertEqual(len(a&b),2)
        self.assertNotEqual(a,b)
        for c in ("DFR_ANR_primary","DFR_ANR_alt_isoform"):
            ids=list(CASES[c]["source"])
            self.assertEqual(len({x.rsplit("_i",1)[0] for x in ids}),3)

    def test_same_published_transcripts_are_selected_by_exact_id(self):
        for c in CASES:
            records,metadata=source_alignment_input(c,self.src,self.refs)
            self.assertEqual(len(records),len(metadata))
            self.assertEqual({m["source_id"] for m in metadata.values()
                              if not m["is_reference"]},set(CASES[c]["source"]))
            self.assertTrue(all(m["aa_segment_length"]>=80
                                for m in metadata.values() if not m["is_reference"]))
            self.assertTrue(all(m.get("protein") is None
                                for m in metadata.values()))

    @unittest.skipUnless(shutil.which("mafft"),"External MAFFT required for integration")
    def test_mafft_and_small_phylogeny_pipeline_cannot_claim_flower_reexpression(self):
        r=run_panel(
            DATA/"cirsium_japonicum_2020_26_source_sequences_v6.fasta",
            DATA/"cirsium_japonicum_2020_26_sequence_provenance_v6.json",
            DATA/"DFR_ANR_ANS_FLS_reference_proteins_v5.fasta",
            DATA/"DFR_ANR_ANS_FLS_reference_provenance_v5.json",
            bootstraps=12)
        self.assertFalse(r["phylogenetic_orthology_established"])
        self.assertFalse(r["molecular_phenotype_reexpression_demonstrated"])
        self.assertFalse(r["source_ryukyu_2026_sequences_included"])
        self.assertEqual(len(r["trees"]),3)
        self.assertEqual(sum(len(x["queries"]) for x in r["trees"].values()),7)
        self.assertTrue(all(x["min_observed_shared_sites_across_pairs"]>=60
                            for x in r["trees"].values()))


if __name__=="__main__":
    unittest.main()
