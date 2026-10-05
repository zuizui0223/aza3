#!/usr/bin/env python3
from __future__ import annotations
import csv, json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
DOC=ROOT/"docs"/"ADAPTIVE_VARIANT_IDENTIFICATION_PROGRAM_V1.md"
LADDER=ROOT/"data"/"planning"/"adaptive_variant_evidence_ladder_v1.csv"
SYSTEMS=ROOT/"data"/"planning"/"adaptive_variant_focal_systems_v1.csv"
MASTER=ROOT/"docs"/"AZA3_NATURE_SCALE_MASTER_PLAN_V1.md"
README=ROOT/"README.md"

def need(t,x):
    if x not in t:
        raise AssertionError(f"missing: {x}")

def read_csv(p):
    with p.open(encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f))

def main():
    d=DOC.read_text(encoding="utf-8")
    m=MASTER.read_text(encoding="utf-8")
    r=README.read_text(encoding="utf-8")
    ladder=read_csv(LADDER)
    systems=read_csv(SYSTEMS)

    assert [x["level_id"] for x in ladder]==["AV0","AV1","AV2","AV3","AV4","AV5","AV6"]
    assert ladder[1]["allowed_term"]=="trait-associated_locus"
    assert "adaptive_variant_candidate" in ladder[5]["allowed_term"]
    assert ladder[6]["allowed_term"]=="reusable_adaptive_genomic_module"

    sys={x["system_id"]:x for x in systems}
    assert set(sys)=={"C1","C2","C3","C4"}
    assert "20-30 W + 20-30 BP" in sys["C1"]["new_sampling_or_data"]
    assert "75 brevicaule + 75 irumtiense" in sys["C2"]["new_sampling_or_data"]
    assert "n=60" in sys["C3"]["new_sampling_or_data"]
    assert "n=40" in sys["C4"]["new_sampling_or_data"]

    for token in (
        "For the word **adaptive**, AV1 + AV2 + AV5 are mandatory",
        "reusable adaptive genomic module",
        "C. japonicum var. takaoense W/BP polymorphism",
        "C. brevicaule / C. irumtiense Ryukyu system",
        "Do not perform a pooled colour GWAS across unrelated species",
        "Presence of DFR/ANS or another pathway gene is not sufficient",
        "same locus is reused",
        "many independent molecular routes",
    ):
        need(d,token)

    need(m,"Adaptive-variant evidence ladder")
    need(m,"selection scan")
    need(m,"adaptive-variant molecular flagship")
    need(r,"docs/ADAPTIVE_VARIANT_IDENTIFICATION_PROGRAM_V1.md")
    need(r,"adaptive_variant_evidence_ladder_v1.csv")

    forbidden=(
        "candidate gene = adaptive gene",
        "selection scan proves adaptation",
        "DFR/ANS proves adaptation",
    )
    low=(d+"\n"+m).casefold()
    for x in forbidden:
        if x.casefold() in low:
            raise AssertionError(f"overclaim: {x}")

    print(json.dumps({
        "status":"ok",
        "adaptive_variant_levels":7,
        "molecular_flagship":"takaoense_W_BP",
        "phaseA_counts_aligned":True,
        "adaptive_claim_requires_fitness":True,
        "reusable_module_requires_replication":True
    },indent=2))

if __name__=="__main__":
    main()
