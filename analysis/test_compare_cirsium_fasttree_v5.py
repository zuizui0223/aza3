"""Independent FastTree vs NJ unrooted placement checks; no orthology claims."""
from io import StringIO
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis"))
from compare_cirsium_fasttree_v5 import compare_methods,run_fasttree

class FastTreeAlternativeTopologyTest(unittest.TestCase):
    def test_missing_program_or_input_fails_closed(self):
        with self.assertRaises((OSError,ValueError,RuntimeError)):
            run_fasttree(Path("/nonexistent/cirsium_gene_family_alignment.fasta"))

    @unittest.skipUnless(shutil.which("FastTree") or shutil.which("fasttree"),
                         "FastTree requires external executable")
    def test_small_artificial_protein_tree_validated(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"ref.fasta"
            p.write_text(
                "".join(f">{label}\n{seq}\n" for label,seq in (
                    ("DFR1","ACDEFGHIKLMNPQRSTVWY"*7),
                    ("DFR2","ACDEFGHIKLMNPQRSTVWY"*7),
                    ("ANR1","YCDEFGHIKLMNPQRSTVWA"*7),
                    ("ANR2","YCDEFGHIKLMNPQRSTVWA"*7),
                    ("CJ","ACDEFGHIKLMNPQRSTVWY"*7))),
                encoding="utf-8")
            tree,newick,digest=run_fasttree(p)
            self.assertEqual(len(tree.get_terminals()),5)
            self.assertTrue(newick.endswith(";"))
            self.assertEqual(len(digest),64)

    @unittest.skipUnless(shutil.which("FastTree") or shutil.which("fasttree"),
                         "FastTree requires external executable")
    def test_comparison_requires_real_nj_and_mafft_files(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(FileNotFoundError):
                compare_methods(Path(d)/"unavailable.json",Path(d))
        # A full real-source integration run is checked by the CI workflow
        # after the MAFFT/NJ output files exist; no invented sequence results.

if __name__=="__main__":
    unittest.main()
