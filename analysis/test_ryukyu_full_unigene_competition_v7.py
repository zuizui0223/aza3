"""Offline counterexamples for full 104890-transcriptome competition against 26 genes.

NO real Ryukyu raw reads are loaded by these tests. The bounded network
workflow supplies independent public 8MiB R1 observations when successful.
"""
from hashlib import sha256
from io import BytesIO
from pathlib import Path
import sys
import tempfile
import unittest
import zipfile
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis/anto_empirical"))

from ryukyu_full_unigene_competition_v7 import (
    prepare_competing_reference,parse_competitor_sam,build_index,run
)
from cirsium_source_homology_v3 import read_fasta


def toy_source():
    txt=(b"Published unigene supplementary sheet\n"
         b">TR1|c0_g1_i1 len=120\n"+b"ACGT"*30+b"\n"
         b">TR2|c0_g2_i1 len=120\n"+b"TGCA"*30+b"\n"
         b">TR3|c0_g3_i1 len=120\n"+b"AGCT"*30+b"\n")
    data=BytesIO()
    with zipfile.ZipFile(data,"w",compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("ao0c04001_si_002.txt",txt)
    return data.getvalue(),txt


class WholeUnigeneSourceTests(unittest.TestCase):
    def test_source_fasta_composition_and_all_reference_members(self):
        z,txt=toy_source()
        with tempfile.TemporaryDirectory() as folder:
            f=Path(folder)/"selected.fa"
            out=Path(folder)/"full.fa"
            f.write_text(">TR1|c0_g1_i1\n"+("ACGT"*30)+"\n"
                         ">TR2|c0_g2_i1\n"+("TGCA"*30)+"\n")
            stats=prepare_competing_reference(
                z,f,out,expected_count=3,
                expected_sha=sha256(txt).hexdigest(),
                expected_bytes=len(txt),subset_n=2)
            self.assertEqual(stats["full_publisher_unigene_count"],3)
            self.assertEqual(stats["candidate_sequence_ids_verified_in_full_reference"],2)
            seq=read_fasta(out)
            self.assertEqual(len(seq),5)
            self.assertEqual(seq["TR1|c0_g1_i1"],"ACGT"*30)
            self.assertEqual(seq["TR3|c0_g3_i1"],"AGCT"*30)
            self.assertIn("DECOY_TR1|c0_g1_i1",seq)

    def test_whole_source_hash_drift_fails(self):
        z,txt=toy_source()
        with tempfile.TemporaryDirectory() as folder:
            f=Path(folder)/"selected.fa"
            out=Path(folder)/"full.fa"
            f.write_text(">TR1|c0_g1_i1\n"+("ACGT"*30)+"\n")
            with self.assertRaisesRegex(ValueError,"hash/count"):
                prepare_competing_reference(z,f,out,expected_count=3,
                       expected_sha="f"*64,expected_bytes=len(txt),subset_n=1)

    def test_missing_registered_selected_transcript_fails(self):
        z,txt=toy_source()
        with tempfile.TemporaryDirectory() as folder:
            f=Path(folder)/"selected.fa"
            out=Path(folder)/"full.fa"
            f.write_text(">TR4|c0_g4_i1\n"+("ACGT"*30)+"\n")
            with self.assertRaisesRegex(ValueError,"source candidates"):
                prepare_competing_reference(z,f,out,expected_count=3,
                       expected_sha=sha256(txt).hexdigest(),
                       expected_bytes=len(txt),subset_n=1)

    def test_reference_sequence_swap_does_not_silently_match(self):
        z,txt=toy_source()
        with tempfile.TemporaryDirectory() as folder:
            f=Path(folder)/"selected.fa"
            out=Path(folder)/"full.fa"
            f.write_text(">TR1|c0_g1_i1\n"+("TGCA"*30)+"\n")
            with self.assertRaisesRegex(ValueError,"mismatch"):
                prepare_competing_reference(z,f,out,expected_count=3,
                       expected_sha=sha256(txt).hexdigest(),
                       expected_bytes=len(txt),subset_n=1)

    def test_anonymous_full_reference_reads_are_not_forced_into_colour_genes(self):
        def sam(ref,mapq=60,span=100,nm=4):
            cols=["read_"+ref,"0",ref,"1",str(mapq),f"{span}M",
                  "*","0","0","A"*span,"I"*span,f"NM:i:{nm}"]
            return "\t".join(cols)
        lines=[
            sam("TR1|c0_g1_i1"),
            sam("TR1000|c0_g9_i2"),
            sam("DECOY_TR1|c0_g1_i1"),
            sam("TR1|c0_g1_i1",mapq=0),
            sam("TR1|c0_g1_i1",span=50),
            sam("TR1|c0_g1_i1",nm=40),
        ]
        result=parse_competitor_sam(lines,{"TR1|c0_g1_i1"})
        counts=result["primary_counts"]
        self.assertEqual(counts["strict_original_26_candidate"],1)
        self.assertEqual(counts["strict_other_104864_public_transcript"],1)
        self.assertEqual(counts["strict_composition_shuffled_decoy"],1)
        self.assertEqual(counts["low_mapq_or_coverage_or_identity_original_26_candidate"],3)
        self.assertEqual(result["source_candidate_strict_26_counts"],{"TR1|c0_g1_i1":1})
        self.assertEqual(result["composition_shuffled_decoy_strict_counts"],{"DECOY_TR1|c0_g1_i1":1})

    def test_unmapped_sam_is_not_false_negative_gene_call(self):
        lines=["r1\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*"]
        result=parse_competitor_sam(lines,{"TR1|c0_g1_i1"})
        self.assertEqual(result["primary_counts"],{})
        self.assertTrue(result["observed_not_gene_function"])

    def test_unregistered_source_rna_runs_rejected(self):
        with self.assertRaises(ValueError):
            run(runs=["SRR11111"],cap_bytes=1024*1024)
        with self.assertRaises(ValueError):
            run(runs=["SRR35152725"],cap_bytes=9*1024*1024)

    def test_reference_index_failure_is_not_called_gene_absence(self):
        with tempfile.TemporaryDirectory() as folder:
            fa=Path(folder)/"missing.fa"
            index=Path(folder)/"x.mmi"
            with patch("shutil.which",return_value=None):
                with self.assertRaisesRegex(RuntimeError,"minimap2"):
                    build_index(fa,index)


if __name__=="__main__":
    unittest.main()
