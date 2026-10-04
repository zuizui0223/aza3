#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
README=ROOT/"README.md"
DOC=ROOT/"docs"/"AZAMI_EAZAMI_TO_AZA3_NATURE_ORIGIN_HANDOFF_V1.md"
CON=ROOT/"data"/"contracts"/"aza3_nature_origin_handoff_v1.json"

def need(text, token):
    if token not in text:
        raise AssertionError(f"missing: {token}")

def main():
    r=README.read_text(encoding="utf-8")
    d=DOC.read_text(encoding="utf-8")
    c=json.loads(CON.read_text(encoding="utf-8"))

    need(r,"Why could a young thistle radiation repeatedly generate so much reproductive-phenotype diversity in so little evolutionary time?")
    need(r,"docs/AZAMI_EAZAMI_TO_AZA3_NATURE_ORIGIN_HANDOFF_V1.md")
    need(r,"repeated de novo mutation")
    need(r,"ancestral standing variation")
    need(r,"introgressive reuse")
    need(r,"polyploid/homeolog/regulatory reuse")

    for token in (
        "46,276 observations",
        "259 taxa",
        "orientation: **4–6**",
        "phyllary posture: **exactly 3**",
        "stickiness: **exactly 5**",
        "**0/3 trait pairs**",
        "**99/376 = 26.3%**",
        "**0/324**",
        "**0/21**",
        "How could so much phenotypic change happen so quickly?",
        "Where does the capacity for rapid phenotypic innovation come from?",
        "repeated de novo mutation",
        "ancestral standing variation",
        "introgressive reuse",
        "polyploid/homeolog/regulatory reuse",
        "phenotype space → repeated evolutionary history → genomic source of evolvability → fitness consequence",
    ):
        need(d,token)

    assert c["contract_version"]=="aza3_nature_origin_handoff_v1"
    assert c["source_repositories"]["azami"]["source_main_sha"]=="b463dbc938e454baa54b2e480f7f1d8badd7795d"
    assert c["source_repositories"]["eazami"]["source_main_sha"]=="f406c32ac6d4a8da49777d7973a92f48146bd6df"
    assert c["source_repositories"]["aza3_previous_main"]["source_main_sha"]=="d8becd9d2ca221e88b37832bc573ab0072733b18"
    assert c["source_repositories"]["azami"]["frozen_context"]["observations"]==46276
    assert c["source_repositories"]["azami"]["frozen_context"]["taxa"]==259
    assert c["source_repositories"]["eazami"]["frozen_context"]["dominant_radiation_concepts"]==36
    assert c["source_repositories"]["eazami"]["frozen_context"]["orientation_minimum_changes"]=="4-6"
    assert c["source_repositories"]["eazami"]["frozen_context"]["robust_shared_transition_localization"]=="0/3"
    assert c["source_repositories"]["eazami"]["frozen_context"]["orientation_historical_regime_match"]=="99/376"
    assert set(c["competing_models"])=={
        "repeated_de_novo_mutation",
        "ancestral_standing_variation",
        "introgressive_reuse",
        "polyploid_homeolog_or_regulatory_reuse",
    }
    assert len(c["claim_boundaries"])>=7

    banned=[
        "standing variation caused repeated capitulum states",
        "reticulation caused the radiation",
        "ploidy increases evolvability",
        "the Japanese Cirsium radiation is an adaptive radiation",
    ]
    low=(r+"\n"+d).lower()
    for b in banned:
        if b.lower() in low:
            raise AssertionError(f"overclaim present: {b}")

    print(json.dumps({
        "status":"ok",
        "project_entrypoint":"nature_scale_evolvability",
        "source_commits_frozen":True,
        "competing_genomic_models":4,
        "claim_boundaries_preserved":True
    },indent=2))

if __name__=="__main__":
    main()
