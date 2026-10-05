#!/usr/bin/env python3
from pathlib import Path
import csv, json

ROOT=Path(__file__).resolve().parents[1]
DOC=ROOT/"docs"/"R1B1_OWN_WGS_TRANSFERABILITY_PILOT_V1.md"
CON=ROOT/"data"/"contracts"/"aza3_r1b1_wgs_transferability_pilot_v1.json"
SHEET=ROOT/"data"/"templates"/"r1b1_wgs_pilot_sample_sheet_v1.csv"

def need(t,x):
    if x not in t:
        raise AssertionError(f"missing: {x}")

def main():
    d=DOC.read_text(encoding="utf-8")
    c=json.loads(CON.read_text(encoding="utf-8"))
    header=next(csv.reader(SHEET.open(encoding="utf-8")))

    assert c["status"]=="DESIGN_FROZEN_NOT_COLLECTION_AUTHORIZATION"
    assert c["focal_taxon"]=="Cirsium sieboldii"
    assert c["sample_target"]=={"core":8,"preferred":10}
    assert c["minimum_design"]["populations"]>=2
    assert c["minimum_design"]["within_population_contrast_pairs_target"]>=2
    assert c["sequencing"]["nominal_depth_x"]=="8-10"

    ids=[x["id"] for x in c["references"]]
    assert ids==[
        "heterophyllum_hap1","heterophyllum_hap2",
        "dissectum_hap1","dissectum_hap2","nipponicum"
    ]
    assert c["staged_unblinding"]==[
        "A_phenotype_blind_technical_metrics",
        "B_population_labels",
        "C_capitulum_phenotype_labels"
    ]

    required={
        "individual_id","population_id","voucher_or_image","direct_colour",
        "presentation_angle_deg","phyllary_state","stickiness_state",
        "genome_size_2C","cytotype","leaf_DNA_sample","authorization_record",
        "phenotype_blind_code"
    }
    assert required.issubset(set(header))

    for x in (
        "mixed populations",
        "8–10×",
        "same WGS reads are mapped to hap1 and hap2",
        "Stage A — phenotype-blind",
        "Stage B — population labels",
        "Stage C — phenotype labels",
        "WGS_TRANSFER_GREEN",
        "WGS_TRANSFER_AMBER",
        "WGS_TRANSFER_RED",
        "Do not scale this pilot into the full Phase-A panel",
    ):
        need(d,x)

    forbidden_claims={
        "no_causal_variant","no_adaptive_locus","no_haplotype_age",
        "no_trait_introgression","no_recombination_block_reuse",
        "no_SV_reuse","no_genomic_modularity","no_genotype_phenotype_association"
    }
    assert forbidden_claims.issubset(set(c["claim_ceiling"]))

    print(json.dumps({
        "status":"ok",
        "sample_core":c["sample_target"]["core"],
        "sample_preferred":c["sample_target"]["preferred"],
        "references":ids,
        "unblinding":c["staged_unblinding"],
        "authorization":"design_only",
    },indent=2))

if __name__=="__main__":
    main()
