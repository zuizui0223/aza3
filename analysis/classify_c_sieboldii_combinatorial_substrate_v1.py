#!/usr/bin/env python3
"""Classify the primary C × O C. sieboldii combinatorial-substrate gate.

Primary gate only:
- A1/A2 anthesis observations
- colour white/non-white descriptive state
- continuous floral chroma
- gravity-referenced orientation (0 up, 90 horizontal, 180 down)

This script does not open secondary phyllary pairs. Those require a separately
frozen technical-error/variation rule before use.
"""
from __future__ import annotations
import argparse, csv, json, math
from collections import defaultdict
from pathlib import Path

VALID_STAGES={"A1","A2","A1_early_anthesis","A2_full_anthesis"}

def fnum(x):
    try:
        v=float(x)
        return v if math.isfinite(v) else None
    except Exception:
        return None

def bval(x):
    s=(x or "").strip().casefold()
    if s in {"1","true","t","yes","y","white","w"}: return True
    if s in {"0","false","f","no","n","nonwhite","non-white","coloured","colored","c"}: return False
    return None

def ranks(xs):
    idx=sorted(range(len(xs)), key=lambda i: xs[i])
    out=[0.0]*len(xs)
    i=0
    while i<len(idx):
        j=i+1
        while j<len(idx) and xs[idx[j]]==xs[idx[i]]:
            j+=1
        r=(i+1+j)/2.0
        for k in range(i,j):
            out[idx[k]]=r
        i=j
    return out

def pearson(x,y):
    if len(x)<3: return None
    mx=sum(x)/len(x); my=sum(y)/len(y)
    dx=[a-mx for a in x]; dy=[b-my for b in y]
    den=math.sqrt(sum(a*a for a in dx)*sum(b*b for b in dy))
    if den==0: return None
    return sum(a*b for a,b in zip(dx,dy))/den

def spearman(x,y):
    return pearson(ranks(x),ranks(y))

def classify_group(rows, group_id, group_type):
    complete=[]
    for r in rows:
        if r.get("flower_stage","").strip() not in VALID_STAGES:
            continue
        o=fnum(r.get("orientation_deg_gravity"))
        ch=fnum(r.get("colour_chroma"))
        w=bval(r.get("colour_white_descriptive"))
        if o is None or ch is None or w is None:
            continue
        complete.append((r,o,ch,w))
    n=len(complete)
    orient=[x[1] for x in complete]
    chroma=[x[2] for x in complete]
    white=[x[3] for x in complete]
    orng=(max(orient)-min(orient)) if orient else None
    nwhite=sum(white)
    nnon=n-nwhite
    classes=[]
    cells=defaultdict(int)
    for _,o,_,w in complete:
        oc="upright" if o<=60 else ("nodding" if o>=120 else "intermediate")
        classes.append(oc)
        if oc!="intermediate":
            cells[("white" if w else "nonwhite",oc)]+=1
    nu=classes.count("upright"); nn=classes.count("nodding")
    occupied={f"{c}|{o}":v for (c,o),v in sorted(cells.items()) if v>=3}
    rho=spearman(chroma,orient) if n>=3 else None

    tests={
        "complete_A1_A2_n_pass":n>=20,
        "orientation_range_pass":orng is not None and orng>=30,
        "white_count_pass":nwhite>=3,
        "nonwhite_count_pass":nnon>=3,
        "upright_extreme_count_pass":nu>=3,
        "nodding_extreme_count_pass":nn>=3,
        "three_of_four_cells_pass":len(occupied)>=3,
        "continuous_coupling_pass":rho is not None and abs(rho)<0.8,
    }
    primary_pass=all(tests.values())
    return {
        "group_id":group_id,
        "group_type":group_type,
        "n_complete_A1_A2":n,
        "orientation_range_deg":orng,
        "white_n":nwhite,
        "nonwhite_n":nnon,
        "upright_extreme_n":nu,
        "nodding_extreme_n":nn,
        "intermediate_orientation_n":classes.count("intermediate"),
        "occupied_cells_min3":occupied,
        "occupied_cells_min3_n":len(occupied),
        "spearman_chroma_orientation":rho,
        "tests":tests,
        "primary_CO_pass":primary_pass,
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--input",required=True,type=Path)
    ap.add_argument("--output",required=True,type=Path)
    args=ap.parse_args()
    rows=list(csv.DictReader(args.input.open(encoding="utf-8-sig")))
    if not rows:
        raise SystemExit("empty census")

    bypop=defaultdict(list)
    byset=defaultdict(list)
    for r in rows:
        p=(r.get("population_id") or "").strip()
        s=(r.get("local_set_id") or "").strip()
        if p: bypop[p].append(r)
        if s: byset[s].append(r)

    popres=[classify_group(v,k,"population") for k,v in sorted(bypop.items())]
    setres=[classify_group(v,k,"local_set") for k,v in sorted(byset.items())]

    green=[x for x in popres if x["primary_CO_pass"]]
    # Local-set GREEN requires primary pass AND representation by >=2 population IDs,
    # with each colour and each orientation extreme represented in >=2 populations.
    local_green=[]
    for result in setres:
        if not result["primary_CO_pass"]: continue
        group=[r for r in rows if (r.get("local_set_id") or "").strip()==result["group_id"] and r.get("flower_stage","").strip() in VALID_STAGES]
        pop_white=defaultdict(set); pop_orient=defaultdict(set)
        for r in group:
            p=(r.get("population_id") or "").strip()
            o=fnum(r.get("orientation_deg_gravity")); w=bval(r.get("colour_white_descriptive"))
            if not p or o is None or w is None: continue
            pop_white["white" if w else "nonwhite"].add(p)
            if o<=60: pop_orient["upright"].add(p)
            elif o>=120: pop_orient["nodding"].add(p)
        result["local_set_population_support"]={
            "white_populations":len(pop_white["white"]),
            "nonwhite_populations":len(pop_white["nonwhite"]),
            "upright_populations":len(pop_orient["upright"]),
            "nodding_populations":len(pop_orient["nodding"]),
        }
        result["local_set_confounding_pass"]=all(v>=2 for v in result["local_set_population_support"].values())
        if result["local_set_confounding_pass"]:
            local_green.append(result)

    if green:
        decision="CS_GREEN"
    elif local_green:
        decision="CS_GREEN_LOCAL_SET"
    else:
        # This script cannot distinguish geography-only variation from absence/phenology
        # without complete field disposition; fail closed to NOT_IDENTIFIABLE.
        decision="NOT_IDENTIFIABLE_PRIMARY_CO"

    out={
        "classifier_version":"classify_c_sieboldii_combinatorial_substrate_v1",
        "primary_pair":"floral_colour_x_anthesis_orientation",
        "orientation_convention":"0=up;90=horizontal;180=down; upright<=60; nodding>=120",
        "population_results":popres,
        "local_set_results":setres,
        "decision":decision,
        "decision_boundary":"Automated classifier can positively authorize only CS_GREEN or CS_GREEN_LOCAL_SET for the primary C×O pair. AMBER/RED/secondary-pair decisions require the complete field-disposition audit and cannot be inferred from missing data.",
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(out,indent=2))

if __name__=="__main__":
    main()
