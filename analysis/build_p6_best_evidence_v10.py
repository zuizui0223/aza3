#!/usr/bin/env python3
from __future__ import annotations
import csv
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
P=ROOT/'data'/'planning'
NCBI=P/'chapter3_p6_biosample_locality_audit_v10.csv'
B1=P/'chapter3_p6_public_locality_audit_batch1_v8.csv'
CONS=P/'chapter3_core_conservation_screen_v10.csv'
P6=P/'chapter3_holefill_priority6_add_own_trait_link_v7.csv'
OUT=P/'chapter3_p6_best_evidence_v10.csv'

def rows(path):
 with path.open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f))

def main():
 nc={r['species_binomial']:r for r in rows(NCBI)}
 b1={r['species_binomial']:r for r in rows(B1)}
 cons={r['species_binomial']:r for r in rows(CONS)}
 p6=rows(P6)
 out=[]
 for r in p6:
  sp=r['species_binomial']; n=nc[sp]; b=b1.get(sp)
  if b and b['public_locality_resolution']=='PRECISE_SECTOR':
   resolution='PRECISE_SECTOR'; summary=b['public_locality_summary']; source=b['evidence_source']; rule=b['own_complement_rule']
  elif n['locality_resolution']=='PRECISE_SECTOR':
   resolution='PRECISE_SECTOR'; summary=n['geo_loc_name_public']; source='NCBI_BioSample_public_geo_loc_name'; rule='Choose verified wild sector complementary to public sector after NMNS-range and permission audit'
  elif b:
   resolution=b['public_locality_resolution']; summary=b['public_locality_summary']; source=b['evidence_source']; rule='Do not guess complement; resolve voucher/public locality or use two-own-sector sensitivity design'
  else:
   resolution=n['locality_resolution']; summary=n['geo_loc_name_public']; source='NCBI_BioSample_public_geo_loc_name'; rule='Do not guess complement; resolve voucher/public locality or use two-own-sector sensitivity design'
  mode=cons[sp]['acquisition_mode']
  if mode in {'HISTORICAL_MATERIAL_ONLY','HERBARIUM_OR_MINIMAL_AUTHORIZED_WILD','CONSERVATION_GATE_REQUIRED'}:
   collection='CONSERVATION_MODE_OVERRIDES_SIMPLE_COMPLEMENT'
  elif mode in {'CONSERVATION_REVIEW','LOCALITY_AND_SUBTAXON_GATE'}:
   collection='COMPLEMENT_REQUIRES_CONSERVATION_OR_SUBTAXON_REVIEW'
  elif resolution=='PRECISE_SECTOR': collection='COMPLEMENT_SECTOR_CAN_BE_RESOLVED_AFTER_OCCURRENCE_PERMISSION_AUDIT'
  else: collection='PUBLIC_LOCALITY_UNRESOLVED_NO_COMPLEMENT_GUESS'
  out.append({
   'species_binomial':sp,'moreyra_jpn':r['moreyra_jpns'],'biosample':n['biosample'],
   'best_public_locality_resolution':resolution,'best_public_locality_summary':summary,
   'best_public_locality_source':source,'national_acquisition_mode':mode,
   'own_sampling_decision':collection,'own_complement_rule':rule,
   'exact_coordinates_public_repo':'false'
  })
 with OUT.open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=out[0].keys());w.writeheader();w.writerows(out)
 from collections import Counter
 print('p6_species=26')
 print(Counter(x['best_public_locality_resolution'] for x in out))
 print(Counter(x['own_sampling_decision'] for x in out))
 return 0
if __name__=='__main__':raise SystemExit(main())
