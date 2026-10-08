#!/usr/bin/env python3
"""Enforce Tofts source-footnote and non-pollination interpretation boundaries."""
import csv
import json
from collections import Counter
from pathlib import Path

PATH = Path(__file__).resolve().parents[1] / "data/evidence/eriophorum_microhabitat_role_audit_v1.csv"

def main():
    rows=list(csv.DictReader(PATH.open(encoding="utf-8-sig", newline="")))
    assert rows and all(r["focal_species"] == "Cirsium eriophorum" for r in rows)
    needed=("partner_taxon","life_stage","organ_or_microhabitat",
            "documented_activity","source_precision","admitted_role","notes")
    assert all(all(r.get(k) for k in needed) for r in rows)
    seen=set()
    for r in rows:
        key=(r["partner_taxon"], r["documented_activity"])
        assert key not in seen, "duplicate partner-role assertion"
        seen.add(key)
        if "footnote_three_stars" in r["source_precision"]:
            assert "DO_NOT_TREAT_AS_PRESENT_IN_ERI_OPHORUM" in r["admitted_role"]
            assert "C_palustre_only" in r["documented_activity"]
        if "visitor" in r["admitted_role"]:
            assert "not_verified_pollinator" in r["admitted_role"]
        if r["partner_taxon"]=="Forficula auricularia":
            assert r["documented_activity"]=="shelter"
            assert "no_predation" in r["admitted_role"]
        if r["partner_taxon"]=="Palloptera modesta":
            assert "not_confirmed_predator" in r["admitted_role"]
        if r["partner_taxon"]=="Pteromalus vibulenus":
            assert "Rhinocyllus" in r["organ_or_microhabitat"]
    prohibited=sum("footnote_three_stars" in r["source_precision"] for r in rows)
    assert prohibited==2
    print(json.dumps({
      "status":"PASS_TOFTS_FOOTNOTE_GUILD_GATES",
      "n_source_role_rows":len(rows),
      "prohibited_as_eriophorum_present":prohibited,
      "compilation_or_outside_UK_rows":sum(
          r["source_precision"].startswith("compilation") for r in rows),
      "scope":"bibliographic microhabitat associations, not co-located community census",
      "distinct_roles":dict(Counter(r["admitted_role"] for r in rows))
    },indent=2))

if __name__=="__main__": main()
