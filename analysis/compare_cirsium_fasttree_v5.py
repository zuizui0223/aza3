"""Compare NJ topology against a distinct approximate-ML (FastTree WAG) tree.

Inputs are independently produced MAFFT amino-acid alignments and NJ source
result JSON from Cirsium 2020 evidence, not Ryukyu 2026 genotypes.

FastTree's local SH-like support values, if present, must NOT be interpreted
as 200-bootstrapped NJ support or a protein-function posterior. The key result
is whether the anchored family-side assignment CHANGES with inference method.

No enzyme activity, orthology, trait regain or independent genetic origins
are inferred here. Small gene-family reference sampling constrains generality.
"""
import argparse
from hashlib import sha256
from io import StringIO
import json
from pathlib import Path
import shutil
import subprocess
import sys

from Bio import Phylo, SeqIO

ROOT=Path(__file__).resolve().parents[1]
MODULE=ROOT/"analysis/anto_empirical"
sys.path.insert(0,str(MODULE))
from cirsium_flavonoid_gene_tree_v5 import (
    CASES,load_reference_proteins,safe_name,unrooted_family_split
)

SOURCE=ROOT/"data/anto_source_handoff/molecular"


def run_fasttree(alignment,timeout=90):
    exe=shutil.which("FastTree") or shutil.which("fasttree")
    if exe is None:
        raise RuntimeError("FastTree is not installed; cannot infer the method robustness")
    p=Path(alignment)
    data=p.read_bytes()
    if len(data)>5000000:
        raise ValueError("Unexpectedly large source-family protein alignment")
    aligned=list(SeqIO.parse(p,"fasta"))
    if len(aligned)<5 or len({x.id for x in aligned})!=len(aligned):
        raise ValueError("Protein alignment missing or duplicates sequence labels")
    lengths={len(x.seq) for x in aligned}
    if len(lengths)!=1 or next(iter(lengths))<100:
        raise ValueError("Protein input has unequal or inadequate aligned sequence length")
    run=subprocess.run(
        [exe,"-wag","-gamma","-quiet",str(p)],
        capture_output=True,text=True,timeout=timeout,check=False
    )
    if run.returncode or not run.stdout.strip().endswith(";"):
        raise RuntimeError(
            f"FastTree failed (return {run.returncode}): {run.stderr[-500:]}"
        )
    t=Phylo.read(StringIO(run.stdout.strip()),"newick")
    if {x.name for x in t.get_terminals()}!={x.id for x in aligned}:
        raise ValueError("FastTree labels do not match frozen MAFFT alignment")
    if len(t.get_terminals())!=len(aligned):
        raise ValueError("FastTree lost or duplicated aligned protein reference")
    return t,run.stdout.strip(),sha256(data).hexdigest()


def compare_methods(source_json,alignment_dir,reference_fasta=None,reference_metadata=None):
    source_json=Path(source_json)
    alignment_dir=Path(alignment_dir)
    previous=json.loads(source_json.read_text(encoding="utf-8"))
    if previous.get("analysis_version")!="aza3_cirsium_2020_protein_family_phylogeny_v5":
        raise ValueError("Wrong source NJ version")
    if previous.get("source_ryukyu_2026_sequences_included") or previous.get("phylogenetic_orthology_established"):
        raise ValueError("Unexpected biological inference in input")
    refs=load_reference_proteins(
        reference_fasta or SOURCE/"DFR_ANR_ANS_FLS_reference_proteins_v5.fasta",
        reference_metadata or SOURCE/"DFR_ANR_ANS_FLS_reference_provenance_v5.json"
    )
    if set(previous["trees"])!=set(CASES):
        raise ValueError("Expected separate DFR/ANR, isoform and ANS/FLS tree analyses")
    cases={}
    ncomps=0
    nunchanged=0
    for case,configuration in CASES.items():
        details=previous["trees"][case]
        fasta=alignment_dir/(case+".fasta")
        fast,newick,source_sha=run_fasttree(fasta)
        positives={safe_name("REF",r["accession"]) for r in refs.values()
                   if r["group"]==configuration["positive"]}
        related={safe_name("REF",r["accession"]) for r in refs.values()
                   if r["group"]==configuration["related"]}
        queries={}
        for key,record in details["source_translation"].items():
            if record["is_reference"]:
                continue
            nl=details["queries"][key]["observed_unrooted_NJ_family_split"]
            fl=unrooted_family_split(fast,key,positives,related)
            queries[record["source_id"]]={
                "source_transcript":record["source_id"],
                "NJ_unrooted_reference_family_placement":nl,
                "FastTree_WAG_unrooted_reference_family_placement":fl,
                "classification_stable_across_inference_method":nl==fl,
                "functional_orthology_verified":False,
            }
            ncomps+=1
            nunchanged+=(nl==fl)
        cases[case]={
            "fasttree_inference":"approximate-ML WAG+CAT tree with Gamma20 log-likelihood branch-length rescaling",
            "sh_like_support_not_nj_bootstrap":True,
            "MAFFT_alignment_SHA256":source_sha,
            "source_protein_alignment_length":len(next(iter(SeqIO.parse(fasta,"fasta"))).seq),
            "gene_family_positive_ref_count":len(positives),
            "gene_family_related_ref_count":len(related),
            "FastTree_newick":newick,
            "queries":queries,
        }
    return {
        "analysis_version":"aza3_cirsium_source_family_method_comparison_v5",
        "source":"Park 2020 Cirsium japonicum annotated transcripts, not focal 2026 Ryukyu samples",
        "n_inference_cases":len(cases),
        "n_query_case_comparisons":ncomps,
        "n_identical_NJ_vs_approxML_placements":nunchanged,
        "cases":cases,
        "allele_orthology_certified":False,
        "enzyme_function_verified":False,
        "Ryukyu_white_to_pigmented_regain_demonstrated":False,
        "caution":"Different unrooted anchored family placements indicate methodological sensitivity, not functional distinction. Identical placements can remain unresolved if queries lie between sparse positive/related family anchors."
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--tree-json",type=Path,required=True)
    p.add_argument("--alignment-directory",type=Path,required=True)
    p.add_argument("--output",type=Path)
    args=p.parse_args()
    result=compare_methods(args.tree_json,args.alignment_directory)
    serialized=json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(serialized,encoding="utf-8")
    else:
        print(serialized)
    print("TREE_METHOD_SUMMARY_JSON:"+json.dumps({
        "source":"2020_Cirsium_japonicum",
        "n_inference_cases":result["n_inference_cases"],
        "comparison_count":result["n_query_case_comparisons"],
        "same_placements":result["n_identical_NJ_vs_approxML_placements"],
        "trait_reexpression":result["Ryukyu_white_to_pigmented_regain_demonstrated"],
    },separators=(",",":")))
    for case,details in result["cases"].items():
        for source,r in details["queries"].items():
            print("TREE_METHOD_QUERY_JSON:"+json.dumps({"case":case,**r},separators=(",",":")))
    return result


if __name__=="__main__":
    main()
