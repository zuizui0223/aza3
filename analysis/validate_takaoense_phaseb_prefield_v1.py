#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
DOC=ROOT/"docs"/"TAKAOENSE_PHASEB_PREFIELD_CANDIDATE_REDUCTION_V1.md"
CON=ROOT/"data"/"contracts"/"takaoense_phaseb_prefield_contract_v1.json"
README=ROOT/"README.md"

def need(t,x):
    if x not in t:
        raise AssertionError(f"missing: {x}")

def main():
    d=DOC.read_text(encoding="utf-8")
    c=json.loads(CON.read_text(encoding="utf-8"))
    r=README.read_text(encoding="utf-8")

    expected={
        "FC":"SRR35152718",
        "WY":"SRR35152717",
        "FB":"SRR35152738",
        "TJ":"SRR35152736",
        "NH":"SRR35152735",
        "LT":"SRR35152734",
    }
    assert {x["code"]:x["run"] for x in c["samples"]}==expected
    assert sum(x["morph"]=="W" for x in c["samples"])==3
    assert sum(x["morph"]=="BP" for x in c["samples"])==3
    assert c["av_ladder_role"]=="pre_AV1_candidate_reduction_only"
    assert c["next_population_panel"]["minimum_populations_per_state"]==2
    assert c["next_population_panel"]["mixed_populations_preferred"] is True

    for token in (
        "mean altitude BP = 1,160.67 m",
        "mean altitude W = 357.00 m",
        "BP − W mean difference = 803.67 m",
        "no differential-expression claim from these young-leaf libraries",
        "no causal-gene claim from a 3-vs-3 fixed difference",
        "pre-AV1 candidate set",
        "STRUCTURAL_CODING_CANDIDATE",
        "REGULATORY_ROUTE_FAVOURED",
        "ANCESTRY_CONFOUNDED",
        "NOT_IDENTIFIABLE",
    ):
        need(d,token)

    need(r,"TAKAOENSE_PHASEB_PREFIELD_CANDIDATE_REDUCTION_V1.md")

    print(json.dumps({
        "status":"ok",
        "takaoense_public_samples":6,
        "morph_balance":"3W_3BP",
        "public_panel_role":"pre_AV1_candidate_reduction",
        "altitude_geography_confounding":"explicit",
        "new_population_replication_required":True
    },indent=2))

if __name__=="__main__":
    main()
