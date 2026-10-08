#!/usr/bin/env python3
"""Validate curated source anchors without mistaking article evidence for GloBI."""
import csv
import json
import collections
import hashlib
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PATH=ROOT/"data/evidence/cirsium_capitulum_interaction_source_anchors_v1.csv"
GLOBI_TYPE="versioned_globi_source_review"
LIT_TYPE="primary_literature"

def main():
    content=PATH.read_bytes()
    rows=list(csv.DictReader(content.decode("utf-8-sig").splitlines()))
    assert len(rows)>0
    keys=set()
    for r in rows:
        for k in ("focal_cirsium","partner","role","organ","relation",
                  "evidence_source","source_id","evidence_scope"):
            assert r.get(k), (k,r)
        assert r["focal_cirsium"].startswith("Cirsium ")
        assert "://" in r["source_id"] or r["source_id"].startswith("10.")
        key=tuple(r[k] for k in ("focal_cirsium","partner","role",
                                  "evidence_scope","source_id"))
        assert key not in keys, ("duplicate non-independent label",key)
        keys.add(key)
        if r["globi_exact_record_confirmed"]=="yes":
            assert r["evidence_scope"]==GLOBI_TYPE
            assert r["source_id"].startswith("https://zenodo.org/records/")
        else:
            assert r["globi_exact_record_confirmed"]=="no"
            assert r["evidence_scope"]==LIT_TYPE
        if r["role"]=="head_seed_feeder_parasitoid":
            assert r["organ"]=="capitulum"
            assert r["parasitoid_host"]
            assert r["relation"] in {"parasitoid_of",
                         "parasitizes_head_tephritid_assemblage"}
        if r["partner"]=="Urophora cardui":
            assert r["organ"]=="stem" and r["role"]=="stem_gall_feeder"
        if r["partner"]=="Harpalus rufipes":
            assert r["role"]=="seed_consumption_nonhead"
            assert r["organ"]!="capitulum"
        if r["relation"]=="interactsWith":
            assert r["role"]=="interaction_only"
    assert sum(r["globi_exact_record_confirmed"]=="yes" for r in rows)==7
    # Ensure the two-hop side never treats a parasitoid's host as the plant.
    assert any(r["parasitoid_host"].startswith("tephritid assemblage") for r in rows)
    assert any(r["parasitoid_host"]=="Terellia ruficauda" for r in rows)
    eligible={r["focal_cirsium"] for r in rows if r["organ"]=="capitulum"}
    assert "Cirsium arvense" in eligible and "Cirsium palustre" in eligible
    assert any(r["focal_cirsium"]=="Cirsium sp." and r["role"]=="interaction_only"
               for r in rows), "unresolved genus-level plant records must stay unresolved"
    assert any(r["focal_cirsium"]=="Cirsium vulgare"
               and r["partner"]=="Megachile inermis" and r["role"]=="interaction_only"
               for r in rows), "do not classify generic interactsWith as pollination"
    out={
      "validation_status":"PASS",
      "scientific_label":"SOURCE_BACKED_NONEXHAUSTIVE__GLOBI_RECORDS_SEPARATE_FROM_PRIMARY_LITERATURE",
      "rows":len(rows),
      "source_sha256":hashlib.sha256(content).hexdigest(),
      "globi_confirmed_source_review_rows":sum(r["globi_exact_record_confirmed"]=="yes" for r in rows),
      "primary_literature_rows":sum(r["evidence_scope"]==LIT_TYPE for r in rows),
      "capitulum_specific_rows":sum(r["organ"]=="capitulum" for r in rows),
      "capitulum_species_count":len(eligible),
      "all_plant_names":sorted({r["focal_cirsium"] for r in rows}),
      "head_species_names":sorted(eligible),
      "functional_class_counts":dict(collections.Counter(r["role"] for r in rows)),
      "limitation":"Source-selected, not global GloBI species list. Effects on angle, spine, stickiness, or adaptation are NOT estimated."
    }
    print(json.dumps(out,indent=2))

if __name__=="__main__":
    main()
