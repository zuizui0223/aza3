#!/usr/bin/env python3
"""Build the Chapter 3 one-slot-per-sample acquisition/site-freeze queue.

The public repository stores range sectors and deidentified locality IDs only.
Exact coordinates and sensitive locality descriptions are deliberately prohibited.
The full 128-species national conservation screen determines acquisition mode;
local/prefectural protection and land-manager review remain additional gates.
"""
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PLANNING = ROOT / "data" / "planning"
OUT = PLANNING / "chapter3_exact_site_freeze_queue_v8.csv"
CONSERVATION = PLANNING / "chapter3_core_conservation_screen_v10.csv"

SECTOR_FILES = {
    1: PLANNING / "chapter3_holefill_priority1_sample_sectors_v7.csv",
    2: PLANNING / "chapter3_holefill_priority2_sample_sectors_v7.csv",
    3: PLANNING / "chapter3_holefill_priority3_sample_sectors_v7.csv",
    4: PLANNING / "chapter3_holefill_priority4_sample_sectors_v7.csv",
    5: PLANNING / "chapter3_holefill_priority5_sample_sectors_v7.csv",
    6: PLANNING / "chapter3_holefill_priority6_sample_sectors_v7.csv",
    7: PLANNING / "chapter3_holefill_priority7_sample_sectors_v7.csv",
}

FIELDS = [
    "slot_id","priority","species_binomial","japanese_names","sample_slot","required_for_tree",
    "acquisition_mode","national_redlist_categories","range_sector","nmns_distribution_source","sector_design",
    "public_locality_audit_required","public_locality_audit_status","current_occurrence_status",
    "current_occurrence_evidence_id","current_occurrence_checked_date","conservation_review_status",
    "land_manager_status","collection_permission_status","tissue_permission_status","target_locality_id",
    "private_exact_site_record_status","field_window_status","collaborator_or_collector","freeze_status",
    "exclusion_or_block_reason",
]

def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))

def safe_id(species: str) -> str:
    return species.replace(" ", "_").replace(".", "")

def load_conservation() -> dict[str, dict[str, str]]:
    if not CONSERVATION.exists():
        raise AssertionError("run build_core_conservation_screen_v10.py before building site queue")
    rows=read_rows(CONSERVATION)
    if len(rows)!=128:
        raise AssertionError(f"conservation screen must contain 128 species, got {len(rows)}")
    return {r['species_binomial']:r for r in rows}

def base_row(priority: int, row: dict[str, str], slot: str, required: bool, sector: str, design: str,
             public_audit_required: bool = False) -> dict[str, str]:
    return {
        "slot_id": f"P{priority}_{safe_id(row['species_binomial'])}_{slot}",
        "priority": f"P{priority}","species_binomial": row["species_binomial"],
        "japanese_names": row.get("japanese_names", ""),"sample_slot": slot,
        "required_for_tree": "true" if required else "false",
        "acquisition_mode": "STANDARD_WITH_LOCAL_REVIEW","national_redlist_categories":"NOT_SCREENED",
        "range_sector": sector,"nmns_distribution_source": row.get("nmns_distribution_source", ""),
        "sector_design": design,"public_locality_audit_required": "true" if public_audit_required else "false",
        "public_locality_audit_status": "NOT_REQUIRED" if not public_audit_required else "NOT_CHECKED",
        "current_occurrence_status": "NOT_CHECKED","current_occurrence_evidence_id": "",
        "current_occurrence_checked_date": "","conservation_review_status": "NOT_CHECKED",
        "land_manager_status": "NOT_CHECKED","collection_permission_status": "NOT_CHECKED",
        "tissue_permission_status": "NOT_CHECKED","target_locality_id": "",
        "private_exact_site_record_status": "NOT_CREATED","field_window_status": "NOT_CHECKED",
        "collaborator_or_collector": "","freeze_status": "NOT_FROZEN","exclusion_or_block_reason": "",
    }

def apply_conservation(r: dict[str, str], c: dict[str, str]) -> dict[str, str]:
    mode=c['acquisition_mode']; r['acquisition_mode']=mode
    r['national_redlist_categories']=c['national_redlist_categories']
    if mode=='HISTORICAL_MATERIAL_ONLY':
        r['current_occurrence_status']='NOT_APPLICABLE_EXTINCT_OR_HISTORICAL_MODE'
        r['land_manager_status']='NOT_APPLICABLE_FOR_HISTORICAL_MODE'
        r['collection_permission_status']='WILD_COLLECTION_PROHIBITED_BY_PLAN'
        r['tissue_permission_status']='MUSEUM_OR_HERBARIUM_PERMISSION_REQUIRED'
        r['field_window_status']='NOT_APPLICABLE';r['freeze_status']='WILD_COLLECTION_BLOCKED'
        r['conservation_review_status']='NATIONAL_EX_STATUS_CONFIRMED'
        r['exclusion_or_block_reason']='National Red List EX: use authenticated historical/herbarium material only'
    elif mode!='STANDARD_WITH_LOCAL_REVIEW':
        r['freeze_status']='CONSERVATION_GATE_REQUIRED'
        r['conservation_review_status']='NATIONAL_SCREEN_NONSTANDARD_LOCAL_REVIEW_REQUIRED'
        r['exclusion_or_block_reason']=f'National conservation acquisition mode: {mode}'
    else:
        r['conservation_review_status']='NATIONAL_SCREEN_STANDARD_LOCAL_REVIEW_STILL_REQUIRED'
    return r

def build() -> list[dict[str, str]]:
    out=[]; conservation=load_conservation()
    for priority in range(1,6):
        for row in read_rows(SECTOR_FILES[priority]):
            for slot,key in (("A","sample_A_sector"),("B","sample_B_sector")):
                r=base_row(priority,row,slot,True,row[key],row['sector_design'])
                out.append(apply_conservation(r,conservation[row['species_binomial']]))
    for row in read_rows(SECTOR_FILES[6]):
        r=base_row(6,row,'OWN',True,row['own_sample_sector_rule'],row['sector_design'],public_audit_required=True)
        out.append(apply_conservation(r,conservation[row['species_binomial']]))
    for row in read_rows(SECTOR_FILES[7]):
        r=base_row(7,row,'OPTIONAL_OWN',False,row['optional_own_trait_link_sector_rule'],
                   'OPTIONAL_COMPLEMENT_AFTER_PUBLIC_LOCALITY_AUDIT',public_audit_required=True)
        out.append(apply_conservation(r,conservation[row['species_binomial']]))
    return out

def main() -> int:
    rows=build()
    if len(rows)!=229: raise AssertionError(f"expected 229 slots, got {len(rows)}")
    required=sum(r['required_for_tree']=='true' for r in rows)
    if required!=228: raise AssertionError(f"expected 228 required, got {required}")
    if not any(r['freeze_status']=='WILD_COLLECTION_BLOCKED' for r in rows):
        raise AssertionError('expected at least one EX wild-collection block')
    with OUT.open('w',encoding='utf-8',newline='') as handle:
        w=csv.DictWriter(handle,fieldnames=FIELDS);w.writeheader();w.writerows(rows)
    print('chapter3_exact_site_freeze_queue_v8_built=true')
    print(f'queue_rows={len(rows)}');print(f'required_tree_slots={required}')
    print('optional_trait_link_slots=1')
    print(f"wild_collection_blocked_slots={sum(r['freeze_status']=='WILD_COLLECTION_BLOCKED' for r in rows)}")
    print(f"conservation_gate_slots={sum(r['freeze_status']=='CONSERVATION_GATE_REQUIRED' for r in rows)}")
    print('exact_coordinates_in_public_repo=false')
    return 0
if __name__=='__main__': raise SystemExit(main())
