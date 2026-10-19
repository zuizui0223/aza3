#!/usr/bin/env python3
"""Actual heldout Cirsium orientation labels x Moreyra 2025 S1 nuclear samples.

Original data are pinned EAzami (c9870e91604ed1514344ddc8f765eae5237dbf27):
- heldout_east_asia_trait_state_audit_v1.csv, 19 botanical candidate taxa
- moreyra2025_cirsium_reconciliation_v1/manifest.json + eight source CSV parts

All 258 nuclear-sample metadata rows are SHA256-validated before crosswalk.
Name-only matching can be FALSE (e.g. C. lidjiangense, voucher Turkey and
discordant sequence name). This is a source-availability result, NOT a
reconstructed dated tree or trait-distance memory coefficient.
"""
from __future__ import annotations
import argparse
import copy
import csv
import hashlib
import json
from pathlib import Path

FOCAL={
    "Cirsium argyracanthum":{
        "biosample":"SAMN44017818", "region":"China",
        "orientation":"upward_or_erect"},
    "Cirsium vlassovianum":{
        "biosample":"SAMN34240350", "region":"Mongolia",
        "orientation":"upward_or_erect"},
}
PSEUDOMATCH={
    "taxon":"Cirsium lidjiangense","biosample":"SAMN34240312",
    "region":"Outside_target_region","name":"different_submitted_or_published_name",
}
OK_REGIONS={"China","Mongolia","Russian_Far_East","NE_China_Korea"}
OK_SCOPES={"core_east_asia","northeast_asia_bridge"}
SOURCE_COLUMNS={
    "tree_code","biosample","run","experiment","voucher_and_herbarium",
    "region_class","geographic_location","sra_link_status","scope_class",
    "tree_code_vs_sra_name","name_reconciliation_priority",
}
HELD_COLUMNS={"analysis_taxon","orientation_state","source_url",
              "discovery_overlap","orientation_pool_status"}
SOURCE_SHA="c9870e91604ed1514344ddc8f765eae5237dbf27"

def demand(test,reason):
    if not test:raise ValueError(reason)

def load_csv(path,columns):
    with Path(path).open(encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f)
        demand(columns.issubset(reader.fieldnames or []),
               "ORIGINAL_SOURCE_FIELDS_DRIFTED")
        return list(reader)

def load_manifest_parts(folder):
    folder=Path(folder)
    meta=json.loads((folder/"manifest.json").read_text(encoding="utf-8"))
    demand(meta["contract_version"]=="moreyra2025_cirsium_reconciliation_v1",
           "ORIGINAL_ARCHIVE_VERSION_DRIFT")
    demand(meta["rows"]==258 and len(meta["parts"])==8,
           "ORIGINAL_S1_SPLIT_ROW_EXPECTATION_CHANGED")
    allrows=[]
    for item in meta["parts"]:
        path=folder/item["file"]
        data=path.read_bytes()
        demand(len(data)==item["bytes"] and
               hashlib.sha256(data).hexdigest()==item["sha256"],
               "RAW_NUCLEAR_SAMPLE_PART_CHECKSUM_INVALID")
        rows=load_csv(path,SOURCE_COLUMNS)
        demand(len(rows)==item["rows"],"S1_PART_ROW_COUNT_CHANGED")
        allrows.extend(rows)
    demand(len(allrows)==258,"S1_TOTAL_ROW_COUNT_CHANGED")
    return allrows,meta

def classify_match(row):
    """An accepted source link is NOT an inferred species tree placement."""
    # Use actual herbarium/voucher country, not just a mutable region
    # metadata label. The false C. lidjiangense tree code comes from Turkey.
    voucher=row["voucher_and_herbarium"].lstrip()
    vcountry=(
        (row["region_class"]=="China" and voucher.startswith("China"))
        or (row["region_class"]=="Mongolia" and voucher.startswith("Mongolia"))
        or (row["region_class"]=="Russian_Far_East" and voucher.startswith("Russia"))
        or (row["region_class"]=="NE_China_Korea"
            and voucher.startswith(("China","Korea","South Korea","Russia")))
    )
    conditions=(vcountry
        and row["tree_code_vs_sra_name"]=="exact"
        and row["sra_link_status"]=="linked_runinfo"
        and row["scope_class"] in OK_SCOPES
        and row["region_class"] in OK_REGIONS
        and not row["geographic_location"].startswith("Turkey")
        and bool(row["biosample"])
        and bool(row["voucher_and_herbarium"])
    )
    return "EXACT_PUBLIC_NUCLEAR_SAMPLE" if conditions else "HOLD_NAME_OR_VOUCHER_DISCORDANCE"

def audit(heldout,samples):
    eligible={}
    for row in heldout:
        tax=(row.get("analysis_taxon") or "").strip()
        demand(bool(tax) and tax not in eligible,
               "HELDOUT_TAXON_DUPLICATED_OR_EMPTY" if tax in eligible or not tax else "")
        # Include other heldout state rows in a separate screened count.
        eligible[tax]=row
    demand(len(eligible)==19,"ORIGINAL_EAST_ASIAN_TRAIT_SCREEN_NOT_19")
    candidates={
        k:v for k,v in eligible.items()
        if v["orientation_pool_status"]=="candidate_primary"
    }
    demand(len(candidates)==11,"HELDOUT_ORIENTATION_CANDIDATE_N_CHANGED")
    demand(all(v["discovery_overlap"]=="no"
               and v["orientation_state"] in {
                   "upward_or_erect","downward_or_nodding"}
               for v in candidates.values()),
           "LEAKED_DISCOVERY_OR_AMBIGUOUS_ORIENTATION_CANDIDATE")
    lookup={}
    for r in samples:
        k=r["tree_code"].strip()
        lookup.setdefault(k,[]).append(r)

    by_taxon=[]
    n_clean=0
    n_misleading=0
    n_missing=0
    for tax,h in sorted(candidates.items()):
        found=lookup.get(tax,[])
        verified=[r for r in found
                  if classify_match(r)=="EXACT_PUBLIC_NUCLEAR_SAMPLE"]
        held=[r for r in found
              if classify_match(r)!="EXACT_PUBLIC_NUCLEAR_SAMPLE"]
        if verified: status="EXACT_PUBLIC_NUCLEAR_SAMPLE"
        elif held: status="NAME_MATCH_ONLY_BIOLOGICALLY_UNSAFE"
        else: status="NO_EXACT_TREE_CODE_SAMPLE_RECORD"
        n_clean+=bool(verified)
        n_misleading+=bool(held and not verified)
        n_missing+=not found
        by_taxon.append({
            "taxon":tax,
            "orientation_state":h["orientation_state"],
            "status":status,
            "n_matching_tree_code_rows":len(found),
            "verified_public_nuclear_samples":[
                {"biosample":r["biosample"],"run":r["run"],
                 "voucher":r["voucher_and_herbarium"],
                 "region":r["region_class"],"scope":r["scope_class"],
                 "tree_code_vs_SRA":r["tree_code_vs_sra_name"]}
                for r in verified],
            "rejected_name_only_samples":[
                {"biosample":r["biosample"],"region":r["region_class"],
                 "geographic_location":r["geographic_location"],
                 "tree_code_vs_SRA":r["tree_code_vs_sra_name"],
                 "scope":r["scope_class"]}
                for r in held],
            "n_patristic_distances_measured":0,
        })
    demand((n_clean,n_misleading,n_missing)==(2,1,8),
           "FROZEN_NUCLEAR_TIP_STATE_CROSSWALK_COUNT_CHANGED")
    known={
        x["taxon"]:x for x in by_taxon
        if x["status"]=="EXACT_PUBLIC_NUCLEAR_SAMPLE"
    }
    demand(set(known)==set(FOCAL),"VERIFIED_NUCLEAR_TAXON_IDENTITIES_CHANGED")
    for name,expect in FOCAL.items():
        rec=known[name]
        v=rec["verified_public_nuclear_samples"]
        demand(len(v)==1 and v[0]["biosample"]==expect["biosample"]
               and v[0]["region"]==expect["region"]
               and rec["orientation_state"]==expect["orientation"],
               "VERIFIED_NUCLEAR_BIOSAMPLE_ORIENTATION_DISAGREES")
    bad=next(x for x in by_taxon if x["taxon"]==PSEUDOMATCH["taxon"])
    demand(bad["status"]=="NAME_MATCH_ONLY_BIOLOGICALLY_UNSAFE"
           and len(bad["rejected_name_only_samples"])==1,
           "FALSE_EXTERNAL_TIP_HAS_BEEN_ADMITTED")
    x=bad["rejected_name_only_samples"][0]
    demand(x["biosample"]==PSEUDOMATCH["biosample"]
           and x["region"]==PSEUDOMATCH["region"]
           and x["tree_code_vs_SRA"]==PSEUDOMATCH["name"],
           "SOURCE_NAME_VOUCHER_CONFLICT_RECORD_CHANGED")
    return {
        "status":"REAL_TWO_NUCLEAR_SAMPLE_LINKS_ONE_FALSE_NAME_MATCH_PHYLOGENETIC_MEMORY_HOLD",
        "source_repository":"zuizui0223/EAzami",
        "source_commit":SOURCE_SHA,
        "source_linked_Cirsium_nuclear_sample_rows":len(samples),
        "east_asia_source_trait_rows":len(heldout),
        "heldout_orientation_species":len(candidates),
        "heldout_with_exact_voucher_and_submitted_name_consistent_nuclear_samples":n_clean,
        "heldout_with_false_or_discordant_tree_code_name_only_match":n_misleading,
        "heldout_without_matching_exact_nuclear_sample_record":n_missing,
        "accepted_external_species":sorted(known),
        "accepted_external_upward_or_erect_species":sum(
            x["orientation_state"]=="upward_or_erect" for x in known.values()),
        "accepted_external_downward_or_nodding_species":sum(
            x["orientation_state"]=="downward_or_nodding" for x in known.values()),
        "rejected_false_match":bad["taxon"],
        "by_taxon":by_taxon,
        "memory": {
            "source_type":"sampled public nuclear S1 table, NOT recovered branch-length phylogeny",
            "exact_nuclear_voucher_capitulum_orientation_measured":False,
            "taxon_authority_orientation_not_individual_voucher_trait":True,
            "complete_11_species_nuclear_branch_length_comparison":False,
            "number_external_trusted_tip_label_species":n_clean,
            "independent_external_patristic_tip_pairs_at_most":1,
            "external_2_species_pair_has_trait_state_contrast":False,
            "two_verified_species_both_upward_do_not_identify_memory_gradient":True,
            "external_phylogenetic_memory_rho":None,
            "genus_out_of_lineage_transfer_gain":None,
            "inferential_status":"HOLD_ONLY_TWO_VERIFIED_EXTERNAL_TAXA_NO_DATED_TREE_OR_REPLICATED_BRANCH_DISTANCES",
            "newly_resolved":"Nuclear metadata coverage is nonzero for 2 of 11 orientation candidates, and one apparent third is a name/geography conflict.",
            "not_shown":"Two trusted species do not imply two independent origins of their upward head phenotype.",
        },
        "no_biological_field_or_gene_function_measurement_created":True,
    }

def tamper_checks(heldout,samples):
    base=audit(heldout,samples)
    demand(base["heldout_with_exact_voucher_and_submitted_name_consistent_nuclear_samples"]==2,
           "ORIGINAL_SOURCE_NOT_TWO")
    def must_fail(h,s,needle):
        try:audit(h,s)
        except ValueError as err:
            demand(needle in str(err),"UNEXPECTED_SOURCE_TAMPER_ERROR:"+str(err))
            return
        raise AssertionError("SOURCE_TAMPER_ACCEPTED:"+needle)
    z=copy.deepcopy(samples)
    row=next(x for x in z if x["tree_code"]=="Cirsium argyracanthum")
    row["tree_code_vs_sra_name"]="different_submitted_or_published_name"
    must_fail(heldout,z,"FROZEN_NUCLEAR_TIP_STATE_CROSSWALK_COUNT_CHANGED")
    z=copy.deepcopy(samples)
    row=next(x for x in z if x["tree_code"]=="Cirsium lidjiangense")
    row.update(region_class="China",scope_class="core_east_asia",
               tree_code_vs_sra_name="exact",geographic_location="China")
    demand(classify_match(row)=="HOLD_NAME_OR_VOUCHER_DISCORDANCE",
           "FALSE_GEO_EDIT_OVERRIDE_WOULD_HIDE_TURKISH_VOUCHER")
    row["voucher_and_herbarium"]="China: fabricated source"
    must_fail(heldout,z,"FROZEN_NUCLEAR_TIP_STATE_CROSSWALK_COUNT_CHANGED")
    z=copy.deepcopy(heldout)
    one=next(r for r in z if r["analysis_taxon"]=="Cirsium argyracanthum")
    one["discovery_overlap"]="yes"
    must_fail(z,samples,"LEAKED_DISCOVERY_OR_AMBIGUOUS_ORIENTATION_CANDIDATE")
    z=copy.deepcopy(heldout);z.append(copy.deepcopy(z[0]))
    must_fail(z,samples,"HELDOUT_TAXON_DUPLICATED_OR_EMPTY")
    return {
        "status":"PASS_SYNTHETIC_FALSE_MATCH_AND_MISSING_PHYLOGENY_GATES",
        "n_corruption_cases":4,
        "n_real_nuclear_samples_read_for_original_source":len(samples),
        "new_analyses_of_sequencing_reads":0,
        "new_actual_trait_or_fitness_individuals":0,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-reconciliation-folder",type=Path,required=True)
    p.add_argument("--original-heldout-states",type=Path,required=True)
    p.add_argument("--out",type=Path)
    a=p.parse_args()
    states=load_csv(a.original_heldout_states,HELD_COLUMNS)
    rows,manifest=load_manifest_parts(a.original_reconciliation_folder)
    record={
        "version":"cirsium_moreyra_heldout_nuclear_sample_x_orientation_v1",
        "as_of":"2026-10-19",
        "source_manifest":{"source_artifact_sha256":
            manifest["source_full_reconciliation_sha256"],
            "n_parts":len(manifest["parts"]),
            "n_rows":manifest["rows"]},
        "source_integrity_tests":tamper_checks(states,rows),
        "actual_source_matching":audit(states,rows),
    }
    if a.out:
        a.out.parent.mkdir(parents=True,exist_ok=True)
        a.out.write_text(json.dumps(record,indent=2,ensure_ascii=False)+"\n",
                         encoding="utf-8")
    v=record["actual_source_matching"]
    print(json.dumps({
        "accepted_external_species":v["accepted_external_species"],
        "rejected_false_name_match":v["rejected_false_match"],
        "unmatched_external_species":v["heldout_without_matching_exact_nuclear_sample_record"],
        "memory":v["memory"]["inferential_status"],
        "negative_controls":record["source_integrity_tests"]["status"],
    },indent=2,ensure_ascii=False))

if __name__=="__main__":
    main()
