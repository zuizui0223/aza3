#!/usr/bin/env python3
"""Recover public source dates; compare the frozen visits within calendar strata.

No new occurrences are recruited. Month pooling is phenology alignment, not
contemporaneity. Even year-month matching is not matched survey effort.
Only minimal non-location API metadata is retained; observers are pseudonymous.
"""
from __future__ import annotations
import argparse
from collections import Counter
from datetime import date, datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import re
import time
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from zipfile import ZipFile
import numpy as np
import pandas as pd
from cirsium_flower_partner_regional_turnover_v1 import (
    ARCHIVE_SHA256, OBSERVATION_MEMBER, prepare_events, analyse_grid,
    rarefied_bray_curtis, rectangle_contrasts, summarize_rectangles)

SNAPSHOT='2026-10-08'
API='https://api.inaturalist.org/v1/observations'
UA='aza3-Cirsium-source-date-recovery/2.0'


def compact_observation(o: dict) -> dict:
    """Retain provenance without coordinates, names of observers, or photos."""
    oid=o.get('id')
    if not isinstance(oid,int) or isinstance(oid,bool) or oid<=0:
        raise ValueError('INVALID_API_OBSERVATION_ID')
    uid=(o.get('user') or {}).get('id')
    tax=o.get('taxon') or {}
    return {'event_id':str(oid),'observed_on':o.get('observed_on') or '',
      'observer_key':hashlib.sha256(('inat-observer:'+str(uid)).encode()).hexdigest() if uid is not None else '',
      'current_taxon_id':tax.get('id'),'current_taxon_name':tax.get('name') or '',
      'current_ancestor_names':[t['name'] for t in tax.get('ancestors',[]) if t.get('name')],
      'quality_grade':o.get('quality_grade') or '',
      'captive':o.get('captive'),'obscured':o.get('obscured'),
      'geoprivacy':o.get('geoprivacy'),'taxon_geoprivacy':o.get('taxon_geoprivacy')}


def fetch_metadata(ids: list[str], path: Path) -> dict:
    """At most 200 IDs/request and >=1.1s between calls; no authentication."""
    ids=sorted(set(ids),key=int)
    id_hash=hashlib.sha256(('\n'.join(ids)+'\n').encode()).hexdigest()
    if path.exists():
        saved=json.loads(path.read_text())
        if saved['requested_ids_sha256']!=id_hash:
            raise ValueError('METADATA_CACHE_COHORT_MISMATCH')
        return saved
    out={'version':'cirsium_inat_minimal_metadata_v2',
      'retrieved_at_utc':datetime.now(timezone.utc).isoformat(),
      'api':API,'requested_ids_sha256':id_hash,'requested_count':len(ids),
      'batches':[],'observations':[],'retrieval_errors':[]}
    seen=set()
    for start in range(0,len(ids),200):
        batch=ids[start:start+200]
        url=API+'?'+urlencode({'id':','.join(batch),'per_page':200,'verifiable':'any'})
        raw=None
        for attempt in range(3):
            time.sleep(1.1 if attempt==0 else 5*(attempt+1))
            try:
                with urlopen(Request(url,headers={'User-Agent':UA}),timeout=45) as r:
                    raw=r.read()
                break
            except HTTPError as err:
                if err.code not in {429,500,502,503,504} or attempt==2:
                    out['retrieval_errors'].append({'batch_start':start,'http_status':err.code})
                    break
            except (URLError,TimeoutError) as err:
                if attempt==2:out['retrieval_errors'].append({'batch_start':start,'error_type':type(err).__name__})
        if raw is None:
            # Access failures are not biological absences; retain failed IDs.
            continue
        payload=json.loads(raw)
        results=payload.get('results')
        if not isinstance(results,list):raise ValueError('API_RESULTS_SCHEMA_MISMATCH')
        found=[]
        for ob in results:
            c=compact_observation(ob);oid=c['event_id']
            if oid not in batch or oid in seen:raise ValueError('UNREQUESTED_OR_DUPLICATE_API_ID')
            seen.add(oid);found.append(oid);out['observations'].append(c)
        out['batches'].append({'batch_start':start,'n_requested':len(batch),
          'n_returned':len(found),'reported_total_results':payload.get('total_results'),
          'response_sha256':hashlib.sha256(raw).hexdigest()})
        print(f'metadata batch {start//200+1}: {len(found)}/{len(batch)}',flush=True)
    out['unreturned_ids']=sorted(set(ids)-seen,key=int)
    out['observations'].sort(key=lambda x:int(x['event_id']))
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n')
    return out


def join_dates(events: pd.DataFrame, records: list[dict]) -> tuple[pd.DataFrame,dict]:
    """Never substitute creation time or impute missing/partial dates."""
    x=events.copy();x['event_id']=x.event_id.astype(str)
    if x.event_id.duplicated().any():raise ValueError('DUPLICATE_SOURCE_EVENT')
    ids=[str(r['event_id']) for r in records]
    if len(set(ids))!=len(ids) or not set(ids).issubset(set(x.event_id)):
        raise ValueError('DUPLICATE_OR_UNREQUESTED_METADATA')
    audit=Counter(input_events=len(x),metadata_returned=len(records))
    valid=[]
    for row in records:
        s=row.get('observed_on') or ''
        if not s:
            audit['missing_observed_date']+=1;continue
        if not isinstance(s,str) or not re.fullmatch(r'\d{4}-\d{2}-\d{2}',s):
            audit['invalid_or_partial_date']+=1;continue
        try:d=date.fromisoformat(s)
        except ValueError:
            audit['invalid_or_partial_date']+=1;continue
        if s>SNAPSHOT:
            audit['after_source_snapshot']+=1;continue
        valid.append({**row,'event_id':str(row['event_id']), 'date':s,
          'year':str(d.year),'month':f'{d.month:02d}','year_month':s[:7],
          'quarter':str((d.month-1)//3+1),'observer_key':row.get('observer_key') or ''})
    audit['dated_valid']=len(valid)
    audit['metadata_unreturned']=len(x)-len(records)
    for key in ['missing_observed_date','invalid_or_partial_date','after_source_snapshot']:
        audit.setdefault(key,0)
    cols=['event_id','date','year','month','year_month','quarter','observer_key']
    v=pd.DataFrame(valid) if valid else pd.DataFrame(columns=cols)
    return x.merge(v,on='event_id',how='inner',validate='one_to_one'),dict(audit)


def add_cells(events: pd.DataFrame) -> pd.DataFrame:
    x=events.copy()
    la=np.minimum(np.floor((x.lat+90)/2),89).astype(int)
    lo=np.floor(((x.lon+180)%360)/2).astype(int)
    x['cell']=la.astype(str)+'_'+lo.astype(str)
    return x


def reduce_observer_days(events: pd.DataFrame) -> pd.DataFrame:
    """One record per observer x host x visitor x cell x day; no zero filling."""
    x=events.copy()
    key=x.observer_key.fillna('').astype(str)
    x['_observer_unit']=key.where(key.ne(''),'missing-event-'+x.event_id.astype(str))
    x=x.sort_values('event_id').drop_duplicates(
        ['_observer_unit','focal_taxon','partner_taxon','cell','date'])
    return x.drop(columns='_observer_unit')


def season_summary(events: pd.DataFrame, stratum: str|None, draw: int=10,
                   restrict_cells: list[str]|None=None) -> dict:
    """Keep the same two hosts, two cells AND time stratum within rectangles."""
    x=events.copy()
    if restrict_cells is not None:x=x[x.cell.isin(restrict_cells)]
    x['_stratum']=x[stratum].astype(str) if stratum else 'all'
    frames=[];coverage=[]
    for label,group in x.groupby('_stratum',sort=True):
        c=pd.crosstab([group.focal_taxon,group.cell],group.partner_taxon)
        c=c[c.sum(axis=1)>=draw]
        if c.empty:continue
        d=rarefied_bray_curtis(c.to_numpy(),draw)
        rect=rectangle_contrasts(d,list(c.index),min_shared_cells=3)
        coverage.append({'stratum':label,'eligible_host_cells':len(c),
          'eligible_records':int(c.to_numpy().sum()),'n_rectangles':len(rect)})
        if not rect.empty:
            rect['stratum']=label;frames.append(rect)
    out={'stratum':stratum or 'none','draw':draw,'input_records':len(x),
      'stratum_coverage':coverage,'n_host_pairs':0,'n_reused_rectangles':0,
      'equal_host_pair_mean_delta':None,'host_pairs':[],
      'interpretation':'DESCRIPTIVE_SOURCE_RECORD_PATTERN_NOT_CAUSAL_OR_STANDARDIZED_SURVEY'}
    if not frames:return out
    r=pd.concat(frames,ignore_index=True)
    # Within pair: equal weight across admitted strata, not rectangle count.
    means=r.groupby(['host_a','host_b','stratum'],sort=True)[
       ['same_host_across_regions','different_hosts_within_regions','regional_minus_host']].mean().reset_index()
    ps=[]
    for (a,b),g in means.groupby(['host_a','host_b'],sort=True):
        rr=r[(r.host_a==a)&(r.host_b==b)]
        cells=sorted(set(rr.cell_a)|set(rr.cell_b));loo={}
        for cell in cells:
            z=rr[(rr.cell_a!=cell)&(rr.cell_b!=cell)]
            # This deletion sensitivity reuses eligible original comparisons.
            loo[cell]=float(z.groupby('stratum').regional_minus_host.mean().mean()) if len(z) else None
        ps.append({'host_a':a,'host_b':b,'n_strata':len(g),
           'strata':g.stratum.tolist(),'n_shared_regions':len(cells),
           'n_reused_rectangles':len(rr),
           'same_host_across_regions':float(g.same_host_across_regions.mean()),
           'different_hosts_within_regions':float(g.different_hosts_within_regions.mean()),
           'delta':float(g.regional_minus_host.mean()),
           'by_stratum':g.to_dict('records'),'leave_one_region_out_delta':loo})
    out.update(n_host_pairs=len(ps),n_reused_rectangles=len(r),host_pairs=ps,
        equal_host_pair_mean_delta=float(np.mean([p['delta'] for p in ps])))
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--archive',type=Path,required=True)
    p.add_argument('--metadata',type=Path,required=True)
    p.add_argument('--out-dir',type=Path,required=True)
    args=p.parse_args()
    if hashlib.sha256(args.archive.read_bytes()).hexdigest()!=ARCHIVE_SHA256:
        raise ValueError('ARCHIVE_SHA256_MISMATCH')
    with ZipFile(args.archive) as z:
        raw=z.read(OBSERVATION_MEMBER);source=json.loads(z.read('summary.json'))
    if hashlib.sha256(raw).hexdigest()!=source['csv_sha256']:
        raise ValueError('SOURCE_CSV_HASH_MISMATCH')
    events,admission=prepare_events(pd.read_csv(io.BytesIO(raw),dtype=str,keep_default_na=False))
    if len(events)!=3243:raise ValueError('FROZEN_COHORT_DRIFT')
    meta=fetch_metadata(events.event_id.tolist(),args.metadata)
    dated,recovery=join_dates(events,meta['observations'])
    events=add_cells(events);dated=add_cells(dated)
    baseline,_,_=analyse_grid(events,2,10)
    core=baseline['complete_primary_host_core']
    core_cells=core.get('cells',[]);core_hosts=core.get('hosts',[])
    scenarios={}
    for name,df,stratum,draw in [
      ('original_all_dates_missing',events,None,10),
      ('recovered_dates_unstratified',dated,None,10),
      ('primary_same_calendar_month',dated,'month',10),
      ('strict_same_year_month',dated,'year_month',10),
      ('same_calendar_quarter',dated,'quarter',10),
      ('same_month_daily_observer_dedup',reduce_observer_days(dated),'month',10),
      ('same_month_draw5_sensitivity',dated,'month',5),
      ('same_year_month_draw5_sensitivity',dated,'year_month',5)]:
        scenarios[name]=season_summary(df,stratum,draw)
    q=dated[dated.cell.isin(core_cells)&dated.focal_taxon.isin(core_hosts)].copy()
    common_core={name:season_summary(q,stratum,10) for name,stratum in [
      ('dated_unstratified',None),('same_calendar_month','month'),('same_year_month','year_month')]}
    seasonal=(q.groupby(['focal_taxon','cell','month'],sort=True).size().rename('n').reset_index())
    if len(seasonal):
        seasonal['record_fraction_within_host_cell']=seasonal.n/seasonal.groupby(['focal_taxon','cell']).n.transform('sum')
    observer_counts=dated.observer_key[dated.observer_key.ne('')].value_counts()
    concordant=dated.apply(lambda x:x.partner_taxon==x.get('current_taxon_name') or x.partner_taxon in x.get('current_ancestor_names',[]),axis=1) if len(dated) else pd.Series([],dtype=bool)
    scenarios['same_month_current_visitor_concordant']=season_summary(dated.loc[concordant],'month',10)
    result={'version':'cirsium_visit_dates_season_matched_v2',
      'source_snapshot':SNAPSHOT,'archive_sha256':ARCHIVE_SHA256,
      'minimal_metadata_sha256':hashlib.sha256(args.metadata.read_bytes()).hexdigest(),
      'metadata_retrieved_at_utc':meta['retrieved_at_utc'],
      'original_admission':admission,'date_recovery':recovery,
      'retrieval_errors':meta['retrieval_errors'],
      'dated_year_counts':dated.year.value_counts().sort_index().to_dict(),
      'dated_month_counts':dated.month.value_counts().sort_index().to_dict(),
      'current_visitor_name_or_ancestor_concordant':int(concordant.sum()),
      'observer_summary':{'n_nonmissing_observers':len(observer_counts),
        'records_with_missing_observer':int(dated.observer_key.eq('').sum()),
        'top_one_share':float(observer_counts.iloc[0]/len(dated)) if len(observer_counts) else None,
        'top_ten_share':float(observer_counts.head(10).sum()/len(dated)) if len(dated) else None,
        'records_after_daily_observer_dedup':len(reduce_observer_days(dated))},
      'scenarios':scenarios,'original_common_core':{'cells':core_cells,'hosts':core_hosts,
        'n_dated_records':len(q),'comparisons':common_core},
      'claim_limits':[
        'The comparison was declared before recovering dates, but remains exploratory.',
        'Month matching pools different years; year-month still pools different days and observers.',
        'Unreturned or invalid-date observations are not absences; current API labels can differ from the archived snapshot.',
        'Daily observer deduplication is not a correction for sampling effort, visibility or nondetection.',
        'Calendar strata can change eligible hosts and regions; changed deltas are not a causal adjustment effect.',
        'Reused rectangles and shared hosts are not independent replicates; no significance test is claimed.',
        'Visits are not effective pollination; no head morphology, fitness, re-expression, or collection approval is measured.'
      ]}
    args.out_dir.mkdir(parents=True,exist_ok=True)
    (args.out_dir/'result.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    seasonal.to_csv(args.out_dir/'original_core_month_composition.csv',index=False)
    print(json.dumps({'recovery':recovery,'scenarios':{k:{f:v[f] for f in ['input_records','n_host_pairs','n_reused_rectangles','equal_host_pair_mean_delta']} for k,v in scenarios.items()}},indent=2))

if __name__=='__main__':main()
