"""Synthetic FASTQ fixtures and real SOURCE metadata-only validation.

All synthetic R1 tests explicitly cannot establish real Ryukyu homolog reads.
Only an independently successful bounded HTTP 206 live PR job can produce a
source-backed biological *sequence similarity* report.
"""
import gzip
import io
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis/anto_empirical"))
from ryukyu_R1_prefix_kmer_pilot_v6 import (
    CAP_MAX,K,MANIFEST,TARGETS,PROVENANCE,allowed_seed,
    reverse_complement,canonical,checked_targets,indexes,
    seed_index,dinucleotide_preserving_shuffle,validate_manifest,
    fetch_range,partial_gzip,fastq_records,supporting_groups,
    scan_prefix,run_six
)


def compressed_records(records):
    text=b"".join(
        b"@fixture."+str(i).encode()+b"/1\n"+seq.encode()+b"\n+\n"+b"I"*len(seq)+b"\n"
        for i,seq in enumerate(records))
    return gzip.compress(text,compresslevel=6,mtime=0)


class DummyResponse:
    def __init__(self,content,status,content_range):
        self._content=io.BytesIO(content)
        self.status=status
        self.headers={"Content-Range":content_range}
    def __enter__(self):return self
    def __exit__(self,*args):return False
    def read(self,n):return self._content.read(n)
    def getcode(self):return self.status


class SourceAndSequenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source=checked_targets()
        cls.prepared=indexes()
        cls.runs=validate_manifest()

    def test_six_real_provenance_records_not_flower_samples(self):
        self.assertEqual(len(self.runs),6)
        self.assertEqual({x["taxon_as_experiment_titled"] for x in self.runs},
                         {"Cirsium brevicaule","Cirsium irumtiense"})
        self.assertTrue(all(x["tissue"]=="young leaf" and
                            not x["per_individual_corolla_pigmentation_quantified"]
                            for x in self.runs))
        self.assertEqual(self.prepared["source_transcripts"],26)
        self.assertEqual(self.prepared["source_assembly_groups"],23)

    def test_exact_31mer_reverse_complements_are_equivalent(self):
        seq="ACGTTGCACTTGACGGTACGTCCTAGTACGA"
        self.assertEqual(len(seq),31)
        self.assertTrue(allowed_seed(seq[:K]))
        self.assertEqual(canonical(seq[:K]),canonical(reverse_complement(seq[:K])))

    def test_dinucleotide_shuffle_conserves_base_and_dinucleotide_counts(self):
        import random
        from collections import Counter
        s=self.source["TR39766|c0_g1_i2"]
        s="".join(x for x in s if x in "ACGT")
        q=dinucleotide_preserving_shuffle(s,random.Random(987))
        self.assertEqual(len(s),len(q))
        self.assertEqual(Counter(s),Counter(q))
        self.assertEqual(Counter(zip(s,s[1:])),Counter(zip(q,q[1:])))
        self.assertEqual(s[0],q[0])
        self.assertEqual(s[-1],q[-1])

    def test_homology_candidate_requires_two_nongap_seeds_on_same_read(self):
        source=self.source["TR39766|c0_g1_i2"]
        grp="TR39766|c0_g1"
        matches=[]
        for i in range(max(1,len(source)-150)):
            seq=source[i:i+150].replace("U","T")
            if len(seq)!=150 or any(c not in "ACGT" for c in seq):continue
            r=supporting_groups(seq,b"I"*len(seq),self.prepared["real"])
            if r["unique_confirmed_group"]==grp:
                matches.append(seq)
                break
        self.assertTrue(matches,"The real source should provide some positive synthetic exact-match fixture")
        q=supporting_groups(matches[0],b"!"*150,self.prepared["real"])
        self.assertIsNone(q["unique_confirmed_group"],"Quality threshold must not accept Q0 bases")

    def test_dinucleotide_background_does_not_overlap_real_reference_index(self):
        self.assertFalse(set(self.prepared["real"])&set(self.prepared["decoy"]))
        self.assertGreater(len(self.prepared["real"]),100)
        self.assertGreater(len(self.prepared["decoy"]),100)

    def test_invalid_manifest_source_or_excessive_budget_fails(self):
        with self.assertRaises(ValueError):
            validate_manifest(max_per_run=CAP_MAX+1)
        original=json.loads(MANIFEST.read_text())
        original["runs"][0]["R1_https_url"]="http://untrusted.example.org/r1.gz"
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/"fake.json"
            path.write_text(json.dumps(original))
            with self.assertRaisesRegex(ValueError,"allowlisted"):
                validate_manifest(path,4096)

    def test_http_200_full_multi_gb_response_must_fail_before_read(self):
        run=self.runs[0]
        data=compressed_records(["ACGT"*35]*20)
        with patch("urllib.request.urlopen",return_value=DummyResponse(
            data,200,"")):
            with self.assertRaisesRegex(ValueError,"REMOTE_RANGE_UNAVAILABLE"):
                fetch_range(run,cap=256)

    def test_genuine_partial_http_206_is_bounded(self):
        run=self.runs[0]
        raw=compressed_records(["ACGTTGCA"*20]*400)
        cap=384
        self.assertGreater(len(raw),cap)
        with patch("urllib.request.urlopen",return_value=DummyResponse(
            raw,206,f"bytes 0-{cap-1}/{run['R1_full_file_bytes']}")):
            obtained=fetch_range(run,cap=cap)
        self.assertEqual(obtained,raw[:cap])
        text,eof=partial_gzip(obtained)
        self.assertGreater(len(text),0)
        self.assertFalse(eof)
        record=list(fastq_records(text))
        self.assertGreater(len(record),0)

    def test_invalid_content_range_header_is_rejected(self):
        run=self.runs[0]
        with patch("urllib.request.urlopen",return_value=DummyResponse(
            compressed_records(["ACGT"*45]*10),206,
            f"bytes 1-255/{run['R1_full_file_bytes']}")):
            with self.assertRaisesRegex(ValueError,"REMOTE_RANGE_HEADER_INVALID"):
                fetch_range(run,cap=256)

    def test_truncated_fastq_never_invents_last_record(self):
        raw=compressed_records(["ACTG"*40]*50)
        for end in (len(raw)-3,len(raw)-7):
            prefix,eof=partial_gzip(raw[:end])
            self.assertFalse(eof)
            entries=list(fastq_records(prefix))
            self.assertLessEqual(len(entries),50)
            self.assertTrue(all(len(seq)==len(q) for seq,q in entries))

    def test_prefix_source_matches_do_not_prove_absence_or_flower_reexpression(self):
        source=self.source["TR39766|c0_g1_i2"]
        seq=next(source[i:i+150] for i in range(len(source)-149)
                 if all(c in "ACGT" for c in source[i:i+150]) and
                 supporting_groups(source[i:i+150],b"I"*150,
                                   self.prepared["real"])["unique_confirmed_group"]=="TR39766|c0_g1")
        payload=compressed_records([seq]*100)
        out=scan_prefix(payload[:-8],self.prepared["real"],self.prepared["decoy"])
        self.assertEqual(out["reads_scanned"],100)
        self.assertEqual(out["source_assembly_groups_supported"]["TR39766|c0_g1"],100)
        self.assertFalse(out["full_file_md5_verified"])
        self.assertFalse(out["gzip_end_of_stream_observed"])
        self.assertEqual(out["counts"]["reads_with_single_group_2nonoverlap_31mers"],100)

    def test_uninformative_random_reads_never_imply_gene_deletion(self):
        rand=("TGCA"*35)[:140]
        full=compressed_records([rand]*12)
        d=scan_prefix(full,self.prepared["real"],self.prepared["decoy"])
        self.assertEqual(d["reads_scanned"],12)
        self.assertEqual(d["counts"].get("reads_with_single_group_2nonoverlap_31mers",0),0)
        self.assertTrue(d["full_file_md5_verified"] is False)
        self.assertFalse(bool(d["source_assembly_groups_supported"]))

    def test_mocked_six_source_pipeline_flags_all_results_as_prefix_only(self):
        seq=("ACGT"*40)
        raw=compressed_records([seq]*200)
        with patch("ryukyu_R1_prefix_kmer_pilot_v6.fetch_range",return_value=raw[:512]):
            d=run_six(limit_bytes=512)
        self.assertEqual(d["number_of_successfully_read_R1_prefixes"],6)
        self.assertFalse(d["ancestral_pigment_reexpression_confirmed"])
        self.assertEqual(d["corolla_expression_observations"],0)
        self.assertTrue(all(r["status"]=="REAL_R1_PREFIX_SCANNED" for r in d["runs"]))
        self.assertEqual(d["total_strict_candidate_reads"],0)


if __name__=="__main__":
    unittest.main()
