#!/usr/bin/env python3
from pathlib import Path
import csv, json

ROOT=Path(__file__).resolve().parents[1]
DOC=ROOT/"docs"/"GATE0B_PHENOTYPE_RATE_READINESS_V1.md"
CON=ROOT/"data"/"contracts"/"aza3_gate0b_rate_estimand_v1.json"
REC=ROOT/"data"/"planning"/"aza3_gate0b_trait_recovery_priority_v1.csv"
AUDIT=ROOT/"data"/"evidence"/"aza3_gate0b_japan_trait_coverage_audit_v1.json"
QUEUE=ROOT/"data"/"planning"/"aza3_gate0b_lowcost_recovery_queue_v1.csv"
MASTER=ROOT/"docs"/"AZA3_NATURE_SCALE_MASTER_PLAN_V1.md"

def need(t,x):
    if x not in t:
        raise AssertionError(f"missing: {x}")

def main():
    d=DOC.read_text(encoding="utf-8")
    m=MASTER.read_text(encoding="utf-8")
    c=json.loads(CON.read_text(encoding="utf-8"))
    a=json.loads(AUDIT.read_text(encoding="utf-8"))
    with QUEUE.open(encoding="utf-8-sig",newline="") as f:
        qrows=list(csv.DictReader(f))
    with REC.open(encoding="utf-8-sig",newline="") as f:
        rows=list(csv.DictReader(f))

    cov=c["current_coverage"]
    assert cov["japan38_paper_concepts"]==38
    assert cov["exact_japan38_trait_concepts_strict_spatial"]==14
    assert cov["paper_concepts_represented_at_binomial_level"]==17
    assert cov["distinct_japan38_trait_taxa_n_ge_10"]==6
    assert cov["distinct_japan38_binomials_absent_from_current_azami_source_pool"]==18

    assert c["primary_axes"]==[
        "presentation_angle",
        "floral_lightness",
        "floral_chroma",
        "head_elongation",
        "head_compactness"
    ]
    assert c["primary_estimand"]["name"]=="total_multivariate_rate_ratio"
    assert c["primary_estimand"]["formula"]=="R_total = sum_m sigma2_J,m / sum_m sigma2_BG,m"
    assert c["primary_estimand"]["null"]=="R_total = 1"

    need(d,"only 6 distinct Japan38 trait taxa with at least 10 strict-spatial observations")
    need(d,"R_total = sum_m sigma2_J,m / sum_m sigma2_BG,m")
    need(d,"COVERAGE_NOT_IDENTIFIABLE")
    need(d,"Do not claim acceleration from the six >=10-observation Japanese taxa alone.")
    need(m,"Gate 0B is currently **not authorized**")
    need(m,"R_total = sum_m sigma2_J,m / sum_m sigma2_BG,m")

    taxa={r["taxon"] for r in rows}
    assert "Cirsium dipsacolepis" in taxa
    assert "Cirsium alpicola" in taxa
    assert "Cirsium yuki-uenoanum" in taxa
    assert "Cirsium effusum" in taxa
    assert len(rows)>=20
    assert a["detector_positive_japan38_binomials"]==18
    assert a["replication_threshold_counts"]=={
        "n_ge_1":18,"n_ge_2":13,"n_ge_3":11,"n_ge_5":11,
        "n_ge_10":9,"n_ge_20":6,"n_ge_50":5
    }
    assert len(a["low_replication_taxa"])==7
    assert qrows[0]["queue_id"]=="Q0"
    assert any(x["taxon"]=="Cirsium dipsacolepis" and x["queue_id"]=="Q1" for x in qrows)

    print(json.dumps({
        "status":"ok",
        "gate0b_authorized":False,
        "exact_japan38_current":14,
        "n_ge_10_current":6,
        "absent_source_pool":18,
        "primary_axes":5,
        "primary_estimand":"R_total",
        "recovery_rows":len(rows),
        "exhaustive_n_ge_10":9,
        "low_replication_taxa":7
    },indent=2))

if __name__=="__main__":
    main()
