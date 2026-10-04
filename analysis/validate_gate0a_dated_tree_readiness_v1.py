#!/usr/bin/env python3
from __future__ import annotations
import csv, json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
DOC=ROOT/"docs"/"GATE0A_DATED_TREE_READINESS_AND_RECOVERY_V1.md"
RECHECK=ROOT/"docs"/"GATE0A_ROUTE_A_RECHECK_20261004.md"
CON=ROOT/"data"/"contracts"/"aza3_gate0a_dated_tree_readiness_v1.json"
ACC=ROOT/"data"/"contracts"/"aza3_gate0a_route_b_acceptance_v1.json"
ROUTES=ROOT/"data"/"planning"/"aza3_gate0a_recovery_routes_v1.csv"
MASTER=ROOT/"docs"/"AZA3_NATURE_SCALE_MASTER_PLAN_V1.md"

def need(t,x):
    if x not in t:
        raise AssertionError(f"missing: {x}")

def main():
    d=DOC.read_text(encoding="utf-8")
    rr=RECHECK.read_text(encoding="utf-8")
    c=json.loads(CON.read_text(encoding="utf-8"))
    a=json.loads(ACC.read_text(encoding="utf-8"))
    m=MASTER.read_text(encoding="utf-8")
    with ROUTES.open(encoding="utf-8-sig",newline="") as f:
        rows=list(csv.DictReader(f))

    assert c["current_state"]["japan38_scaffold"]=="undated_phylogram_substitutions_per_site"
    assert c["current_state"]["public_exact_moreyra_dated_tree_recovered"] is False
    assert c["current_state"]["gate0b_rate_test_authorized"] is False
    assert c["pass_states"]==["DATED_ENSEMBLE_RECOVERED","DATED_ENSEMBLE_REBUILT_VALIDATED"]

    assert [r["route_id"] for r in rows]==["A","B","C"]
    assert rows[0]["status"]=="ACTIVE_FIRST"
    assert rows[1]["status"]=="READY_IF_A_FAILS"
    assert rows[2]["status"]=="FALLBACK"

    need(rr,"AUTHOR_SOURCE_DATED_TREE_NOT_RECOVERED")
    need(rr,"Route B preparation")
    need(d,"Never rescale the substitutions/site Japan38 tree")
    need(d,"Never assign published ages to guessed nodes")
    need(m,"Gate 0A")
    need(m,"Only a recovered or independently validated dated ensemble opens Gate 0B.")

    landmarks={x["landmark"]:x for x in a["temporal_landmark_validation"]}
    assert landmarks["Cirsium crown"]["published_interval_ma"]==[7.2,12.2]
    assert landmarks["dominant Japanese radiation / jump-dispersal node"]["published_interval_ma"]==[1.7,3.6]
    assert landmarks["Cirsium dipsacolepis separate arrival"]["published_interval_ma"]==[0.4,2.2]
    assert landmarks["Cirsium lineare East Asia-to-Japan expansion"]["published_interval_ma"]==[0.7,2.7]
    assert a["gate0a_pass_rule"]["gate0b_rate_test_authorized_only_after_pass"] is True

    print(json.dumps({
        "status":"ok",
        "routeA":"not_recovered",
        "routeB":"prepared_with_acceptance_contract",
        "routeC":"fallback_only",
        "gate0b_authorized":False,
        "exceptional_evolvability_claim":"not_authorized"
    },indent=2))

if __name__=="__main__":
    main()
