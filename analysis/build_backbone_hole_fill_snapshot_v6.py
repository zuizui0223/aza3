#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json,re
from pathlib import Path
import pandas as pd

INFRA={"var.","subsp.","ssp.","f."}

def clean(x):
    s="" if x is None else str(x)
    s=s.replace("（","(").replace("）",")")
    return re.sub(r"\s+"," ",s).strip()

def binomial(x):
    t=clean(x).replace("C. ","Cirsium ",1).split()
    return " ".join(t[:2]) if len(t)>=2 else clean(x)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--universe",required=True)
    p.add_argument("--moreyra-ledger",required=True)
    p.add_argument("--repair-ledger",required=True)
    p.add_argument("--out-dir",required=True)
    a=p.parse_args()

    u=pd.read_csv(a.universe,dtype=str,keep_default_na=False)
    m=pd.read_csv(a.moreyra_ledger,dtype=str,keep_default_na=False)
    r=pd.read_csv(a.repair_ledger,dtype=str,keep_default_na=False)

    m["species_binomial"]=m["paper_taxon_concept"].map(binomial)
    r["species_binomial"]=r["taxon"].map(binomial)
    moreyra_bins=set(m["species_binomial"])
    repair_bins=set(r["species_binomial"])

    # Collapse NMNS infraspecific rows to species-binomial screening units. This is a
    # candidate census, not a final accepted-species list.
    rows=[]
    for b,g in u.groupby("species_binomial",sort=True):
        concepts=sorted(set(g["authority_taxon_concept"]))
        jnames=sorted(x for x in set(g["japanese_name"]) if x)
        distributions=sorted(x for x in set(g["distribution_summary"]) if x)
        represented=b in moreyra_bins
        if not represented:
            slot="H1_ABSENT_FROM_MOREYRA_SPECIES_SCREEN"
            new_n=2
        elif b in repair_bins:
            slot="H2_MOREYRA_REPAIR_SPECIES_SCREEN"
            new_n=2
        else:
            slot="H3_ADEQUATE_PUBLIC_SPECIES_SCREEN"
            new_n=1
        rows.append({
            "species_binomial":b,
            "slot_class":slot,
            "new_samples_nominal":new_n,
            "n_nmns_records":len(g),
            "authority_concepts":" | ".join(concepts),
            "japanese_names":" | ".join(jnames),
            "distribution_summaries":" | ".join(distributions),
            "represented_in_moreyra_binomial_screen":represented,
            "taxonomic_status":"UNRECONCILED_SPECIES_BINOMIAL_SCREEN",
            "sampling_status":"DO_NOT_COLLECT_UNTIL_OPERATIONAL_SPECIES_CENSUS_FREEZES",
        })
    out=Path(a.out_dir); out.mkdir(parents=True,exist_ok=True)
    frame=pd.DataFrame(rows)
    frame.to_csv(out/"nmns_moreyra_species_hole_fill_screen_v6.csv",index=False,encoding="utf-8")
    h1=frame[frame.slot_class.str.startswith("H1")]
    h2=frame[frame.slot_class.str.startswith("H2")]
    h3=frame[frame.slot_class.str.startswith("H3")]
    summary={
        "contract_version":"nmns_moreyra_species_hole_fill_screen_v6",
        "n_nmns_authority_records":int(len(u)),
        "n_unique_nmns_species_binomials":int(frame.species_binomial.nunique()),
        "n_unique_moreyra_species_binomials":int(len(moreyra_bins)),
        "n_repair_species_binomials":int(len(repair_bins)),
        "h1_absent_species_binomial_screen":int(len(h1)),
        "h2_repair_species_binomial_screen":int(len(h2)),
        "h3_adequate_species_binomial_screen":int(len(h3)),
        "nominal_new_samples_if_every_screen_unit_admitted":int(frame.new_samples_nominal.astype(int).sum()),
        "claim_boundary":"Species-binomial hole-fill SCREEN only. NMNS records include infraspecific and unpublished/new-name entries; this output must not be called the final Japanese species list until taxonomic reconciliation is frozen."
    }
    (out/"nmns_moreyra_species_hole_fill_screen_v6.json").write_text(json.dumps(summary,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(summary,indent=2,ensure_ascii=False))

if __name__=="__main__": main()
