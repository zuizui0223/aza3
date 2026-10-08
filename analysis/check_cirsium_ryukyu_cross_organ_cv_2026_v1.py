#!/usr/bin/env python3
"""Source-table validation of 2026 Ryukyu/Taiwan Cirsium organ CV contrasts.

Published mean and SD are encoded directly from Table 1 and Table 2;
these are not same-plant raw measurements, phylogenetic replicates, or
independently established trait-development effects.
Source DOI 10.1186/s12870-026-08097-6, PMC13020037.
"""
import argparse
import json
from pathlib import Path

SOURCE="https://pmc.ncbi.nlm.nih.gov/articles/PMC13020037/"
DOI="10.1186/s12870-026-08097-6"
# name, table, rosette leaf length mean,SD, leaf width mean,SD,
# capitulum length mean,SD, head width mean,SD, floret length mean,SD.
RAW=[
("C. japonicum var. albescens","Table 1",12.50,5.34,4.34,1.97,2.66,.60,1.36,.32,2.36,.19),
("C. brevicaule","Table 1",32.12,10.36,10.25,2.84,3.68,.18,1.79,.14,2.94,.08),
("C. irumtiense","Table 1",31.55,17.91,11.10,5.86,3.40,.65,1.93,.36,2.59,.05),
("C. japonicum var. takaoense white","Table 2",26.99,10.71,7.86,1.87,3.73,.54,1.43,.24,3.15,.35),
("C. japonicum var. takaoense purple","Table 2",35.26,5.85,10.73,1.52,3.19,.72,1.26,.38,3.50,.09),
("C. japonicum var. australe","Table 2",19.15,6.26,5.92,1.67,3.51,.82,1.56,.55,2.38,.31),
]

def build():
    records=[]
    for (name,table,ll,llsd,lw,lwsd,hl,hlsd,hw,hwsd,fl,flsd) in RAW:
        assert min(ll,lw,hl,hw,fl)>0
        v={"taxon_or_morph":name,"source_table":table,
           "rosette_leaf_length_cv":llsd/ll,
           "rosette_leaf_width_cv":lwsd/lw,
           "capitulum_length_cv":hlsd/hl,
           "capitulum_width_cv":hwsd/hw,
           "floret_length_cv":flsd/fl}
        v.update({
          "leaf_head_length_ratio":v["rosette_leaf_length_cv"]/v["capitulum_length_cv"],
          "leaf_head_width_ratio":v["rosette_leaf_width_cv"]/v["capitulum_width_cv"],
          "head_floret_length_ratio":v["capitulum_length_cv"]/v["floret_length_cv"],
          "leaf_floret_length_ratio":v["rosette_leaf_length_cv"]/v["floret_length_cv"]})
        records.append(v)
    count=lambda key:sum(v[key]>1 for v in records)
    summary={"source_doi":DOI,"source_url":SOURCE,
             "n_named_taxa_or_color_morph_groups":len(records),
             "leaf_gt_head_length_groups":count("leaf_head_length_ratio"),
             "leaf_gt_head_width_groups":count("leaf_head_width_ratio"),
             "head_gt_floret_length_groups":count("head_floret_length_ratio"),
             "leaf_gt_floret_length_groups":count("leaf_floret_length_ratio"),
             "source_measurements":"Mean ± SD taken from publication Tables 1 and 2, no individual-level matching or uncertainty/variance re-estimation.",
             "inferential_ceiling":"Descriptive out-of-sample contrast, NOT head adaptation, canalization proof, or lineage independence.",
             "rows":records}
    assert (summary["leaf_gt_head_length_groups"],
            summary["leaf_gt_head_width_groups"],
            summary["head_gt_floret_length_groups"],
            summary["leaf_gt_floret_length_groups"]) == (5,4,6,6)
    return summary

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--out",type=Path)
    args=parser.parse_args()
    summary=build()
    if args.out:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        args.out.write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps({k:v for k,v in summary.items() if k!="rows"},indent=2))

if __name__=="__main__": main()
