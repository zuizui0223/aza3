#!/usr/bin/env python3
from __future__ import annotations
import csv, re, time, urllib.parse, urllib.request, xml.etree.ElementTree as ET
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
P=ROOT/'data'/'planning'
MEMBERSHIP_URL='https://raw.githubusercontent.com/zuizui0223/EAzami/af36bce7a42a7fcfdafde22e4a78b30d93075f23/data/evidence/moreyra2025_japan_38_membership_audit_2026-08-10.csv'
P6=P/'chapter3_holefill_priority6_add_own_trait_link_v7.csv'
OUT=P/'chapter3_p6_biosample_locality_audit_v10.csv'

def read_csv_url(url):
 with urllib.request.urlopen(url,timeout=30) as r: txt=r.read().decode('utf-8')
 return list(csv.DictReader(txt.splitlines()))

def local_rows(path):
 with path.open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f))

def fetch_biosample(acc):
 term=urllib.parse.quote(f'{acc}[Accession]')
 with urllib.request.urlopen(f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=biosample&term={term}',timeout=30) as r:
  root=ET.fromstring(r.read())
 ids=[x.text for x in root.findall('.//Id') if x.text]
 if not ids:return {}
 with urllib.request.urlopen(f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=biosample&id={ids[0]}',timeout=30) as r:
  root=ET.fromstring(r.read())
 attrs={}
 for a in root.findall('.//Attribute'):
  name=a.attrib.get('attribute_name') or a.attrib.get('harmonized_name') or ''
  if name:attrs[name]=a.text or ''
 return attrs

def resolution(geo):
 s=(geo or '').strip()
 if not s:return 'NO_PUBLIC_LOCALITY'
 if re.search(r'pref|province|city|shi|gun|island|mount|mt\.?|honto|shima',s,re.I):return 'PRECISE_SECTOR'
 if ':' in s or ',' in s:return 'REGION_ONLY'
 if s.lower() in {'japan','japan.'}:return 'COUNTRY_ONLY'
 return 'REGION_ONLY'

def main():
 membership=read_csv_url(MEMBERSHIP_URL)
 by_jpn={r['paper_japan_member_id'].replace('JPN_','JPN_'):r for r in membership}
 p6=local_rows(P6)
 out=[]
 for i,r in enumerate(p6):
  jpn=r['moreyra_jpns'].split('|')[0]
  m=by_jpn.get(jpn,{})
  acc=(m.get('biosamples') or '').split('|')[0]
  attrs={}
  error=''
  if acc:
   try: attrs=fetch_biosample(acc)
   except Exception as e:error=type(e).__name__
   time.sleep(0.35)
  geo=attrs.get('geo_loc_name') or attrs.get('geographic location') or attrs.get('geographic location (country and/or sea)') or attrs.get('country') or ''
  latlon=attrs.get('lat_lon') or attrs.get('latitude and longitude') or ''
  out.append({
   'species_binomial':r['species_binomial'],'moreyra_jpn':jpn,'biosample':acc,
   'geo_loc_name_public':geo,'lat_lon_present_public':'true' if latlon else 'false',
   'locality_resolution':resolution(geo),'audit_error':error,
   'complement_rule': 'CHOOSE_COMPLEMENT_AFTER_SECTOR_RECONCILIATION' if resolution(geo)=='PRECISE_SECTOR' else 'DO_NOT_GUESS_COMPLEMENT',
   'public_sensitive_coordinates_copied_to_repo':'false'
  })
 with OUT.open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=out[0].keys());w.writeheader();w.writerows(out)
 print(f'p6_rows={len(out)}')
 from collections import Counter
 for k,v in Counter(x['locality_resolution'] for x in out).items():print(k,v)
 return 0
if __name__=='__main__':raise SystemExit(main())
