"""Offline safety/QC tests for capped live Ryukyu young-leaf RNA R1 prefix pilot."""
from gzip import compress
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import mock
import json
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis/anto_empirical"))
from ryukyu_ena_prefix_homology_v6 import (
    MAX_BYTES,MIN_MAPQ,validate_manifest,read_fasta,
    download_official_r1_prefix,partial_gzip_to_fastq,
    prepare_target_panel,cigar_aligned_span,parse_mapped_sam
)


def sam(name,ref,mapq=60,matched=100,edits=3,flag=0):
    seq="A"*matched
    quals="I"*matched
    return "\t".join((
        name,str(flag),ref,"1",str(mapq),f"{matched}M",
        "*","0","0",seq,quals,f"NM:i:{edits}"
    ))


class CappedRyukyuSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.records,cls.dna,cls.provenance=validate_manifest()

    def test_only_six_real_source_specimens_and_official_ena_R1(self):
        self.assertEqual(len(self.records),6)
        self.assertEqual(len({r["biosample_accession"] for r in self.records}),6)
        self.assertEqual(len({r["source_specimen_code"] for r in self.records}),6)
        self.assertEqual(len(self.dna),26)
        for r in self.records:
            self.assertTrue(r["R1_https_url"].startswith("https://ftp.sra.ebi.ac.uk/"))
            self.assertTrue(r["R1_https_url"].endswith(r["run_accession"]+"_1.fastq.gz"))
            self.assertFalse(r["whole_run_complete"])
            self.assertFalse(r["per_individual_corolla_pigmentation_quantified"])
            self.assertEqual(r["tissue"],"young leaf")

    def test_budget_and_nonofficial_sources_fail_closed(self):
        with TemporaryDirectory() as d:
            file=Path(d)/"raw.fastq.gz"
            record=self.records[0]
            for x in (MAX_BYTES+1,0,-10):
                with self.assertRaises(ValueError):
                    download_official_r1_prefix(record["R1_https_url"],
                                                record["R1_full_file_bytes"],x,file)
            with self.assertRaises(ValueError):
                download_official_r1_prefix(
                    "https://insecure.example/vol1/fastq/invalid.fq.gz",
                    record["R1_full_file_bytes"],1024,file)

    def test_http_200_does_not_download_entire_1gb_run(self):
        record=self.records[0]
        class Unbounded:
            status=200
            headers={"Content-Range":""}
            def __enter__(self):return self
            def __exit__(self,*args):return None
            def read(self,*args):raise AssertionError("Should not read unbounded response")
        with TemporaryDirectory() as d:
            with mock.patch("urllib.request.urlopen",return_value=Unbounded()):
                with self.assertRaisesRegex(ValueError,"did not honor"):
                    download_official_r1_prefix(
                        record["R1_https_url"],record["R1_full_file_bytes"],
                        1024,Path(d)/"bad.gz")

    def test_partial_gzip_ignores_trailing_incomplete_record(self):
        line_set=[]
        for i in range(13):
            line_set.extend([
                f"@test_2026_{i}\n".encode(),("ACGT"*38+"\n").encode(),
                b"+\n",("I"*152+"\n").encode()])
        fastq_bytes=b"".join(line_set)+b"@broken\nACCTG\n+\n"
        raw=compress(fastq_bytes)
        with TemporaryDirectory() as d:
            source=Path(d)/"first_gzip_member_partial.gz"
            source.write_bytes(raw[:-8]) # no gzip stream end; all read records complete
            output=Path(d)/"reads.fq"
            qc=partial_gzip_to_fastq(source,output)
            self.assertEqual(qc["complete_valid_reads"],13)
            self.assertEqual(qc["total_nucleotide_bases"],13*152)
            self.assertFalse(qc["full_fastq_file_validated"])
            self.assertFalse(qc["representative_random_sample"])
            self.assertNotIn("@broken",output.read_text(encoding="utf-8"))

    def test_gzip_with_prefix_truncation_produces_no_fake_mapped_reads(self):
        with TemporaryDirectory() as d:
            gz=Path(d)/"bad.gz"
            gz.write_bytes(b"\x1f\x8b\x08"+b"\x00"*100)
            with self.assertRaises((ValueError,zlib.error,Exception)):
                partial_gzip_to_fastq(gz,Path(d)/"out.fq")

    def test_reference_plus_shuffled_decoys_keep_source_ids(self):
        with TemporaryDirectory() as d:
            file=Path(d)/"ref.fasta"
            result=prepare_target_panel(self.dna,self.provenance,file)
            self.assertEqual(result["n_real_publisher_transcripts"],26)
            self.assertEqual(result["n_sequence_shuffle_decoys"],26)
            fasta=read_fasta(file)
            self.assertEqual(len(fasta),52)
            self.assertEqual({r for r in fasta if not r.startswith("DECOY_")},set(self.dna))
            self.assertFalse(any(fasta["DECOY_"+k]==fasta[k] for k in self.dna))
            for k in self.dna:
                self.assertEqual(sorted(fasta[k].replace("U","T").replace("-","N")),
                                 sorted(fasta["DECOY_"+k]))

    def test_primary_mapq_length_edit_and_decoy_gates(self):
        target="TR39766|c0_g1_i2"
        data=[
            sam("x1",target,mapq=60,matched=110,edits=5),
            sam("x2",target,mapq=60,matched=110,edits=4),
            sam("x3",target,mapq=0,matched=110,edits=2),
            sam("x4",target,mapq=60,matched=65,edits=0),
            sam("x5",target,mapq=60,matched=110,edits=35),
            sam("decoy","DECOY_"+target,mapq=60,matched=120,edits=2),
            sam("secondary",target,mapq=60,matched=120,edits=0,flag=256),
        ]
        score=parse_mapped_sam(data,"SRR35152730")
        self.assertEqual(score["source_transcript_counts"],{target:2})
        self.assertEqual(score["passing_decoy_reads"],1)
        self.assertEqual(score["passing_real_transcript_reads"],2)
        self.assertGreaterEqual(sum(score["low_or_ambiguous_alignment_counts"].values()),3)
        self.assertFalse(score.get("gene_absence_established",False))

    def test_cigar_minimal_span_validation(self):
        self.assertEqual(cigar_aligned_span("15S130M5I10M2S"),140)
        self.assertEqual(cigar_aligned_span("5H10=5X"),15)
        with self.assertRaises(ValueError):cigar_aligned_span("30B4M")


if __name__=="__main__":
    unittest.main()
