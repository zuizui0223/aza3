"""Offline negative-control tests of bounded Ryukyu public R1 prefix triage."""
from gzip import compress
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase,main
import json
import re
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis/anto_empirical"))
from ryukyu_bounded_read_prefix_v6 import (
    source_manifest,prepare_reference,decode_complete_fastq_prefix,
    parse_minimap2_sam,_cigar_alignment_span,analyze_public_prefixes,
    STRICT_MIN_IDENTITY,
)


def synthetic_fastq(n=1200,read_len=150):
    lines=[]
    for i in range(n):
        seq=("ACTG"*((read_len+3)//4))[:read_len]
        lines.extend([f"@synthetic_{i}",seq,"+","I"*read_len])
    return "\n".join(lines)+"\n"


class SourcePilotSafetyTest(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest=source_manifest()

    def test_manifest_contains_six_actual_public_accessions_with_only_leaf_tissue(self):
        rows=self.manifest["runs"]
        self.assertEqual(len(rows),6)
        self.assertEqual(len({r["biosample_accession"] for r in rows}),6)
        self.assertEqual(len({r["source_specimen_code"] for r in rows}),6)
        self.assertEqual({r["tissue"] for r in rows},{"young leaf"})
        self.assertTrue(all(not r["per_individual_corolla_pigmentation_quantified"]
                            for r in rows))

    def test_reference_has_twenty_six_source_transcripts_and_twenty_six_nulls(self):
        fasta,lookup=prepare_reference()
        self.assertEqual(len(lookup),52)
        self.assertEqual(sum(x["kind"]=="source_transcript" for x in lookup.values()),26)
        self.assertEqual(sum(x["kind"]=="composition_shuffled_decoy" for x in lookup.values()),26)
        self.assertEqual(len(re.findall(r"^>",fasta,flags=re.MULTILINE)),52)
        self.assertEqual((fasta,lookup),prepare_reference())

    def test_fastq_gzip_truncation_emits_only_complete_reads(self):
        gz=compress(synthetic_fastq().encode("ascii"))
        for cut in (gz,gz[:len(gz)-8]):
            fastq,report=decode_complete_fastq_prefix(cut)
            self.assertGreaterEqual(report["read1_prefix_fully_completed_records"],1000)
            self.assertLessEqual(report["read1_prefix_fully_completed_records"],1200)
            self.assertEqual(len(fastq.strip().split("\n"))%4,0)
            self.assertEqual(set(report["read_length_histogram"]),{150})
            self.assertTrue(30<report["mean_quality_phred33_descriptive"]<42)
            self.assertEqual(report["ambiguous_base_fraction_descriptive"],0)
        self.assertFalse(decode_complete_fastq_prefix(gz[:len(gz)-8])[1]["gzip_member_reached_true_eof"])

    def test_corrupt_gzip_magic_rejected(self):
        with self.assertRaisesRegex(ValueError,"gzip"):
            decode_complete_fastq_prefix(b"not a gzip file")

    def test_malformed_fastq_rejected(self):
        s=synthetic_fastq().replace("@synthetic_5","synthetic_5")
        with self.assertRaises(ValueError):
            decode_complete_fastq_prefix(compress(s.encode("ascii")))

    def test_sam_primary_secondary_decoy_and_identity_not_orthology(self):
        _,refs=prepare_reference()
        src=next(x for x in refs if refs[x]["kind"]=="source_transcript")
        dec=next(x for x in refs if refs[x]["kind"]=="composition_shuffled_decoy")
        seq="A"*150
        qual="I"*150
        sam=[
            "@HD\tVN:1.6\tSO:unknown",
            f"r1\t0\t{src}\t10\t60\t150M\t*\t0\t0\t{seq}\t{qual}\tNM:i:2",
            f"r1\t256\t{dec}\t10\t0\t150M\t*\t0\t0\t{seq}\t{qual}\tNM:i:5",
            f"r2\t0\t{dec}\t20\t60\t150M\t*\t0\t0\t{seq}\t{qual}\tNM:i:0",
            f"r3\t4\t*\t0\t0\t*\t*\t0\t0\t{seq}\t{qual}"
        ]
        with TemporaryDirectory() as temp:
            p=Path(temp)/"test.sam"
            p.write_text("\n".join(sam)+"\n",encoding="utf-8")
            report=parse_minimap2_sam(p,refs)
        self.assertEqual(report["n_primary_records"],3)
        self.assertEqual(report["n_primary_unmapped"],1)
        self.assertEqual(report["n_primary_mapped_source"],1)
        self.assertEqual(report["n_primary_mapped_shuffled_decoy"],1)
        self.assertEqual(report["n_secondary_or_supplementary_alignments"],1)
        self.assertEqual(report["n_strict_source_read_alignments"],1)
        self.assertEqual(report["n_strict_shuffled_decoy_alignments"],1)
        self.assertTrue(report["sequence_similarity_not_orthology"])

    def test_invalid_sam_cigar_and_bounds_rejected(self):
        self.assertEqual(_cigar_alignment_span("30S120M"),(120,120))
        self.assertEqual(_cigar_alignment_span("10S100M2I5D"),(107,102))
        with self.assertRaises(ValueError):
            _cigar_alignment_span("50Mbad")
        with self.assertRaises(ValueError):
            analyze_public_prefixes(64*1024*1024)
        with self.assertRaises(ValueError):
            analyze_public_prefixes(1024*1024,only_run="SRR0000000")

    def test_spurious_floral_phenotype_entry_rejected(self):
        doc=json.loads((ROOT/"data/anto_source_handoff/public_sra/ryukyu_six_R1_prefix_source_v6.json").read_text())
        doc["runs"][0]["per_individual_corolla_pigmentation_quantified"]=True
        with TemporaryDirectory() as tmp:
            p=Path(tmp)/"bad.json"
            p.write_text(json.dumps(doc))
            with self.assertRaises(ValueError):
                source_manifest(p)


if __name__=="__main__":
    main()
