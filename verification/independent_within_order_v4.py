#!/usr/bin/env python3
"""Independent brute observation-pair expectation check for aza3 v4.

Does NOT import category_identity_components, block_reference_expectation,
or identity_contrast. Expects SOURCE-prepared dataframe, original pinned
five source strata and the JSON output; never stores event/observer details.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import itertools
import json
import math
from pathlib import Path
import sys
from zipfile import ZipFile
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'analysis'))
from cirsium_flower_partner_regional_turnover_v1 import ARCHIVE_SHA256, OBSERVATION_MEMBER, prepare_events
from cirsium_visit_dates_v2 import add_cells,join_dates

A='Cirsium arvense';B='Cirsium vulgare'
ORDERS={'Lepidoptera','Hymenoptera','Diptera','Coleoptera','Hemiptera','Thysanoptera','Orthoptera'}

def independent_pair_reference(g):
    """Direct ordered-pair expectations for observer × order conditional null."""
    z=g.reset_index(drop=True)
    na=sum(z.focal_taxon.eq(A));nb=sum(z.focal_taxon.eq(B))
    if min(na,nb)<2:raise ValueError('INSUFFICIENT_REAL_HOST_RECORDS')
    names=z.partner_taxon.tolist();hosts=z.focal_taxon.tolist();block=(z.observer_key+'|'+z.visitor_order).tolist()
    bc=z.groupby(['observer_key','visitor_order']).focal_taxon.value_counts()
    p_a={};p_b={};blocktotal={}
    for key,h in zip(block,(z.observer_key+'|'+z.visitor_order)):
        if key in p_a:continue
        obs,order=h.rsplit('|',1)
        n=int(sum(v==key for v in block));v=int(bc.get((obs,order,A),0));
        blocktotal[key]=n;p_a[key]=v/n;p_b[key]=(n-v)/n
    aa=bb=ab=ea=eb=eab=0.
    for i in range(len(z)):
        for j in range(i+1,len(z)):
            if names[i]!=names[j]:continue
            if hosts[i]==A and hosts[j]==A:aa+=1
            elif hosts[i]==B and hosts[j]==B:bb+=1
            else:ab+=1
            first,last=block[i],block[j]
            if first==last:
                n=blocktotal[first]
                ka=int(round(n*p_a[first]));kb=n-ka
                paa=ka*(ka-1)/(n*(n-1)) if n>1 else 0
                pbb=kb*(kb-1)/(n*(n-1)) if n>1 else 0
                pab=2*ka*kb/(n*(n-1)) if n>1 else 0
            else:
                paa=p_a[first]*p_a[last]
                pbb=p_b[first]*p_b[last]
                pab=p_a[first]*p_b[last]+p_b[first]*p_a[last]
            ea+=paa;eb+=pbb;eab+=pab
    obs=.5*(2*aa/(na*(na-1))+2*bb/(nb*(nb-1)))-ab/(na*nb)
    expected=.5*(2*ea/(na*(na-1))+2*eb/(nb*(nb-1)))-eab/(na*nb)
    return obs,expected


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--archive',type=Path,required=True)
    p.add_argument('--metadata',type=Path,required=True)
    p.add_argument('--result',type=Path,required=True)
    args=p.parse_args()
    if hashlib.sha256(args.archive.read_bytes()).hexdigest()!=ARCHIVE_SHA256:
        raise ValueError('SOURCE_ARCHIVE_CHANGED')
    with ZipFile(args.archive) as f:
        raw=f.read(OBSERVATION_MEMBER);manifest=json.loads(f.read('summary.json'))
    if hashlib.sha256(raw).hexdigest()!=manifest['csv_sha256']:
        raise ValueError('ORIGINAL_CSV_CHANGED')
    source,_=prepare_events(pd.read_csv(io.BytesIO(raw),dtype=str,keep_default_na=False))
    saved=json.load(args.metadata.open())
    dated,_=join_dates(source,saved['observations'])
    dated=add_cells(dated)
    dated['visitor_order']=dated.source_taxon_path.map(lambda p: next(iter({s.strip() for s in p.split('|')}&ORDERS)))
    current=json.load(args.result.open())
    prior=current['by_stratum']
    found=[]
    for d in prior:
        z=dated[dated.cell.eq(d['cell'])&dated.month.eq(d['month'])&dated.focal_taxon.isin([A,B])]
        assert len(z)==d['n_records']
        obs,expected=independent_pair_reference(z)
        orig=d['observer_x_order']
        if not math.isclose(obs,orig['observed_u'],rel_tol=0,abs_tol=1e-12):
            raise AssertionError(('OBSERVED_U_DOES_NOT_REPLICATE',d['cell'],obs,orig['observed_u']))
        if not math.isclose(expected,orig['label_reference_mean_u'],rel_tol=0,abs_tol=1e-12):
            raise AssertionError(('REFERENCE_EXPECTATION_NOT_INDEPENDENTLY_REPLICATED',d['cell'],expected,orig['label_reference_mean_u']))
        found.append([d['cell'],len(z),obs,expected,obs-expected])
    obs=sum(r[2] for r in found)/len(found)
    exp=sum(r[3] for r in found)/len(found)
    out=current['primary_observations_and_exact_reference']
    assert abs(obs-out['species_identity_u'])<1e-12
    assert abs(exp-out['observer_plus_order_reference_mean_u'])<1e-12
    assert abs(obs-exp-out['observer_plus_order_residual_u'])<1e-12
    print('PASS INDEPENDENT REAL RECORD-PAIR ENUMERATION',{'n_strata':len(found),'n_records':sum(r[1] for r in found),
      'U':obs,'reference':exp,'excess':obs-exp})

if __name__=='__main__':main()