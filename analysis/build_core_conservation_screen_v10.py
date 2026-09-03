#!/usr/bin/env python3
from __future__ import annotations
import csv
from collections import defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
P=ROOT/'data'/'planning'
SOURCE=P/'chapter3_national_redlist_cirsium_source_v10.csv'
OUT=P/'chapter3_core_conservation_screen_v10.csv'

CORE_FILES=[
 P/'chapter3_holefill_priority1_empty_blocks_v7.csv',
 P/'chapter3_holefill_priority2_broken_tips_v7.csv',
 P/'chapter3_holefill_priority3_kaganoazami_block_v7.csv',
 P/'chapter3_holefill_priority3_sawaazami_block_v7.csv',
 P/'chapter3_holefill_priority3_norikura_series_v7.csv',
 P/'chapter3_holefill_priority3_hamaazami_block_v7.csv',
 P/'chapter3_holefill_priority3_yamaazami_block_v7.csv',
 P/'chapter3_holefill_priority4_partial_blocks_v7.csv',
 P/'chapter3_holefill_priority5_local_gaps_v7.csv',
 P/'chapter3_holefill_priority6_add_own_trait_link_v7.csv',
 P/'chapter3_holefill_priority7_tree_filled_v7.csv',
]

def rows(path):
 with path.open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f))

def main():
 core={}
 for f in CORE_FILES:
  for r in rows(f): core[r['species_binomial']]=r
 if len(core)!=128: raise AssertionError(f'core species drift: {len(core)}')
 red=defaultdict(list)
 for r in rows(SOURCE): red[r['core_species_crosswalk']].append(r)
 out=[]
 for sp in sorted(core):
  hits=red.get(sp,[])
  full=[h for h in hits if h['scope_match'] in {'FULL_SPECIES','AUTHORITY_CROSSWALK_REQUIRED'}]
  infra=[h for h in hits if h['scope_match'].startswith('INFRASPECIFIC')]
  cats='|'.join(sorted({h['category'] for h in hits})) if hits else ''
  if any(h['category']=='EX' and h['scope_match']=='FULL_SPECIES' for h in hits): mode='HISTORICAL_MATERIAL_ONLY'
  elif any(h['category']=='CR' for h in full): mode='HERBARIUM_OR_MINIMAL_AUTHORIZED_WILD'
  elif any(h['category'] in {'EN','VU'} for h in full): mode='CONSERVATION_GATE_REQUIRED'
  elif any(h['category']=='NT' for h in full): mode='CONSERVATION_REVIEW'
  elif infra: mode='LOCALITY_AND_SUBTAXON_GATE'
  else: mode='STANDARD_WITH_LOCAL_REVIEW'
  out.append({
   'species_binomial':sp,
   'national_redlist_categories':cats or 'NOT_LISTED_AT_CORE_SPECIES_LEVEL',
   'full_species_or_crosswalk_hits':'|'.join(h['redlist_scientific_name'] for h in full),
   'infraspecific_hits':'|'.join(h['redlist_scientific_name'] for h in infra),
   'acquisition_mode':mode,
   'national_screen_complete':'true',
   'local_prefectural_and_land_status_still_required':'true',
  })
 with OUT.open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=out[0].keys());w.writeheader();w.writerows(out)
 from collections import Counter
 c=Counter(r['acquisition_mode'] for r in out)
 print('core_species=128')
 for k,v in sorted(c.items()): print(f'{k}={v}')
 return 0
if __name__=='__main__': raise SystemExit(main())
