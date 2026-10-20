#!/usr/bin/env python3
"""Source-pinned, descriptive Cirsium floral-visitor turnover across regions.

Uses only explicit adult insect visitsFlowersOf records in the 2026-10-08
GloBI/iNaturalist archive. Observation IDs are the sampling unit; ambiguous
host/visitor/coordinate assignments are excluded. No feeding, pollination
claim, morphology, fitness, contemporaneous-community or causal inference.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

ARCHIVE_SHA256='ce0513fde35d9c1be42cef87da7a6f4b262e9b4fe237ba0b78e9bf31771b58a4'
OBSERVATION_MEMBER='cirsium_globi_observations.csv'
ADULT={'adult','adult(s)','post-juvenile adult stage','imago'}
INSECT=r'(?:^| \| )(?:Insecta|Lepidoptera|Diptera|Hymenoptera|Coleoptera|Hemiptera|Thysanoptera|Orthoptera)(?: \||$)'
RECT_COLUMNS=['host_a','host_b','cell_a','cell_b','same_host_across_regions',
              'different_hosts_within_regions','regional_minus_host']


def prepare_events(rows: pd.DataFrame) -> tuple[pd.DataFrame,dict]:
    """Retain observation-level explicit adult flower visits, not efficacy."""
    required={'focal_taxon','partner_taxon','focal_position','interaction_type',
      'relation_direction_status','study_title','source_taxon_path',
      'source_specimen_occurrence_id','source_specimen_life_stage',
      'event_date','latitude','longitude'}
    if not required.issubset(rows.columns):
        raise ValueError('MISSING_SOURCE_COLUMNS:'+','.join(sorted(required-set(rows))))
    r=rows.fillna('').copy()
    mask=(r.interaction_type.eq('visitsFlowersOf') & r.focal_position.eq('target')
      & r.relation_direction_status.eq('DIRECTION_COMPATIBLE_NOT_FUNCTION_VERIFIED')
      & r.study_title.str.contains('inaturalist.org/observations/globi-observations-resource-relationships-dwca.zip',regex=False,case=False))
    x=r.loc[mask].copy()
    audit={'raw_archive_rows':len(r),'incoming_iNaturalist_flower_visit_rows':len(x)}
    x['event_id']=x.source_specimen_occurrence_id.str.extract(
        r'^https?://(?:www\.)?inaturalist\.org/observations/([1-9][0-9]*)$',expand=False)
    audit['unresolved_observation_id_rows']=int(x.event_id.isna().sum())
    x=x[x.event_id.notna()].copy()
    x['lat']=pd.to_numeric(x.latitude,errors='coerce')
    x['lon']=pd.to_numeric(x.longitude,errors='coerce')
    # Check ambiguity before named-species/stage filters could conceal conflicts.
    sig=x.groupby('event_id')[['focal_taxon','partner_taxon','lat','lon']].nunique(dropna=False)
    bad=set(sig.index[(sig>1).any(axis=1)])
    audit['ambiguous_observation_ids_excluded']=len(bad)
    x=x[~x.event_id.isin(bad)].copy()
    named=x.focal_taxon.str.fullmatch(r'Cirsium [a-z][a-z-]+') & x.partner_taxon.str.fullmatch(r'[A-Z][a-z]+ [a-z][a-z-]+')
    insect=x.source_taxon_path.str.contains(INSECT,regex=True)
    adult=x.source_specimen_life_stage.str.strip().str.lower().isin(ADULT)
    coord=(np.isfinite(x.lat)&np.isfinite(x.lon)&x.lat.between(-90,90)
      &x.lon.between(-180,180)&~(x.lat.eq(0)&x.lon.eq(0)))
    audit['excluded_flags_nonexclusive']={
       'not_two_binomial_source_names':int((~named).sum()),
       'insect_path_not_supported':int((~insect).sum()),
       'explicit_adult_stage_not_supported':int((~adult).sum()),
       'coordinates_invalid':int((~coord).sum())}
    x=x[named&insect&adult&coord].copy()
    before=len(x)
    x=x.sort_values('event_id').drop_duplicates('event_id')
    audit['duplicate_rows_collapsed_after_admission']=before-len(x)
    audit.update(retained_unique_observations=len(x),
       retained_host_name_strings=int(x.focal_taxon.nunique()),
       retained_visitor_name_strings=int(x.partner_taxon.nunique()),
       retained_rows_without_event_date=int(x.event_date.eq('').sum()),
       taxonomic_names_are_source_strings_not_revalidated_accepted_species=True)
    return x,audit


def rarefied_bray_curtis(counts: np.ndarray,sample_size: int) -> np.ndarray:
    """Exact E[BC] for independent, without-replacement draws of equal size.

    E[min(X,Y)]=sum_{k=1}^r Pr(X>=k)Pr(Y>=k), with each marginal
    hypergeometric. BC=1-sum_s min(X_s,Y_s)/r. No Monte Carlo or p-value.
    Diagonal describes two independent draws, not identity of one draw.
    """
    a=np.asarray(counts)
    if (a.ndim!=2 or a.shape[0]==0 or not np.isfinite(a).all()
        or (a<0).any() or not np.equal(a,np.floor(a)).all()
        or not isinstance(sample_size,(int,np.integer)) or sample_size<1):
        raise ValueError('REQUIRE_NONNEGATIVE_INTEGER_COUNTS_AND_POSITIVE_DRAW_SIZE')
    totals=a.sum(axis=1).astype(int)
    if np.any(totals<sample_size):raise ValueError('INSUFFICIENT_RECORDS_FOR_RAREFACTION')
    tails=np.stack([hypergeom.sf(k-1,totals[:,None],a,sample_size)
                    for k in range(1,sample_size+1)],axis=-1)
    flat=tails.reshape(len(a),-1)
    return np.clip(1-flat@flat.T/sample_size,0,1)


def rectangle_contrasts(distances: np.ndarray,units: list[tuple[str,str]],
                        min_shared_cells: int=3) -> pd.DataFrame:
    """Same two host species in the same two coarse regions (2x2 matching)."""
    d=np.asarray(distances)
    if len(set(units))!=len(units) or d.shape!=(len(units),len(units)):
        raise ValueError('DUPLICATE_UNIT_OR_DISTANCE_SHAPE')
    lookup={u:i for i,u in enumerate(units)}
    hosts=sorted({h for h,c in units})
    cells={h:{c for hh,c in units if hh==h} for h in hosts}
    out=[]
    for a,b in itertools.combinations(hosts,2):
        shared=sorted(cells[a]&cells[b])
        if len(shared)<min_shared_cells:continue
        for c,e in itertools.combinations(shared,2):
            ac,ae,bc,be=[lookup[v] for v in ((a,c),(a,e),(b,c),(b,e))]
            region=(d[ac,ae]+d[bc,be])/2
            host=(d[ac,bc]+d[ae,be])/2
            out.append([a,b,c,e,float(region),float(host),float(region-host)])
    return pd.DataFrame(out,columns=RECT_COLUMNS)


def summarize_rectangles(rect: pd.DataFrame) -> dict:
    if rect.empty:return {'n_host_pairs':0,'n_reused_rectangles':0,'equal_host_pair_mean_delta':None,'host_pairs':[]}
    pairs=rect.groupby(['host_a','host_b'],sort=True).agg(
        n_reused_rectangles=('regional_minus_host','size'),
        same_host_across_regions=('same_host_across_regions','mean'),
        different_hosts_within_regions=('different_hosts_within_regions','mean'),
        delta=('regional_minus_host','mean')).reset_index()
    return {'n_host_pairs':len(pairs),'n_reused_rectangles':len(rect),
      'equal_host_pair_mean_delta':float(pairs.delta.mean()),
      'rectangle_weighted_mean_delta':float(rect.regional_minus_host.mean()),
      'n_host_pairs_positive':int((pairs.delta>1e-12).sum()),
      'n_host_pairs_negative':int((pairs.delta< -1e-12).sum()),
      'host_pairs':pairs.to_dict('records')}


def analyse_grid(events: pd.DataFrame,grid_degrees: int,sample_size: int,
                 offset_degrees: int=0) -> tuple[dict,pd.DataFrame,pd.DataFrame]:
    x=events.copy()
    lat=np.minimum(np.floor((x.lat+90)/grid_degrees),180//grid_degrees-1).astype(int)
    lon=np.floor((((x.lon+180-offset_degrees)%360))/grid_degrees).astype(int)
    x['cell']=lat.astype(str)+'_'+lon.astype(str)
    counts=pd.crosstab([x.focal_taxon,x.cell],x.partner_taxon)
    counts=counts[counts.sum(axis=1)>=sample_size]
    units=[tuple(v) for v in counts.index]
    if len(units):
        distances=rarefied_bray_curtis(counts.to_numpy(),sample_size)
        rect=rectangle_contrasts(distances,units)
    else:rect=pd.DataFrame(columns=RECT_COLUMNS)
    summary=summarize_rectangles(rect)
    common={'Cirsium arvense','Cirsium vulgare'}
    without=rect[~rect.host_a.isin(common)&~rect.host_b.isin(common)]
    summary.update(grid_degrees=grid_degrees,longitude_offset_degrees=offset_degrees,
      records_rarefied_per_host_cell=sample_size,
      n_eligible_host_cells=len(counts),
      n_eligible_hosts=len(set(h for h,c in units)),
      n_eligible_regions=len(set(c for h,c in units)),
      n_records_in_eligible_units=int(counts.to_numpy().sum()),
      excluding_arvense_and_vulgare=summarize_rectangles(without))
    summary['complete_primary_host_core']=complete_core(counts,rect,sample_size)
    core=summary['complete_primary_host_core']
    if core['status']=='COMPLETE_ALL_HOST_COMMON_REGIONS':
        q=x[x.focal_taxon.isin(core['hosts']) & x.cell.isin(core['cells'])].copy()
        orders={'Lepidoptera','Hymenoptera','Diptera','Coleoptera','Hemiptera','Thysanoptera','Orthoptera'}
        def source_order(path):
            found={t.strip() for t in path.split('|')} & orders
            return next(iter(found)) if len(found)==1 else 'UNRESOLVED'
        q['order']=q.source_taxon_path.map(source_order)
        by_order=pd.crosstab([q.focal_taxon,q.cell],q.order)
        order_means=by_order.div(by_order.sum(axis=1),axis=0).groupby(level=0).mean()
        core['source_order_mean_record_fractions']=order_means.to_dict('index')
        by_species=pd.crosstab([q.focal_taxon,q.cell],q.partner_taxon)
        spp_means=by_species.div(by_species.sum(axis=1),axis=0).groupby(level=0).mean()
        core['top_five_recorded_visitors_equal_region_weight']={h:spp_means.loc[h].sort_values(ascending=False).head(5).to_dict() for h in core['hosts']}
    compact=counts.stack().reset_index(name='record_count')
    compact=compact[compact.record_count>0]
    return summary,rect,compact


def complete_core(counts: pd.DataFrame,rectangles: pd.DataFrame,sample_size: int) -> dict:
    """Use ALL hosts in the primary admitted pairs, not a selected subset.

    This is an exploratory common-support sensitivity, not a new primary
    test selected because of a favorable sign. Every common region enters.
    """
    hosts=sorted(set(rectangles.host_a)|set(rectangles.host_b))
    cells=set.intersection(*({c for h,c in counts.index if h==host}
                            for host in hosts)) if hosts else set()
    if len(hosts)<2 or len(cells)<3:
        return {'status':'INSUFFICIENT_ALL_HOST_COMMON_REGIONS',
                'hosts':hosts,'cells':sorted(cells)}
    sub=counts.loc[[u for u in counts.index if u[0] in hosts and u[1] in cells]]
    d=rarefied_bray_curtis(sub.to_numpy(),sample_size)
    rect=rectangle_contrasts(d,list(sub.index))
    return {'status':'COMPLETE_ALL_HOST_COMMON_REGIONS','hosts':hosts,
      'cells':sorted(cells),'n_recorded_visits':int(sub.to_numpy().sum()),
      'n_host_region_units':len(sub),'contrasts':summarize_rectangles(rect),
      'leave_one_common_region_out':{c:summarize_rectangles(
          rect[(rect.cell_a!=c)&(rect.cell_b!=c)]) for c in sorted(cells)}}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--archive',type=Path,required=True)
    p.add_argument('--out-dir',type=Path,required=True)
    args=p.parse_args()
    if hashlib.sha256(args.archive.read_bytes()).hexdigest()!=ARCHIVE_SHA256:
        raise ValueError('PINNED_GLOBI_ARCHIVE_SHA256_MISMATCH')
    with ZipFile(args.archive) as z:
        source=json.loads(z.read('summary.json'))
        raw=z.read(OBSERVATION_MEMBER)
    expected=source['csv_sha256']
    if hashlib.sha256(raw).hexdigest()!=expected:
        raise ValueError('INTERNAL_OBSERVATION_CSV_SHA256_MISMATCH')
    import io
    rows=pd.read_csv(io.BytesIO(raw),dtype=str,keep_default_na=False)
    events,audit=prepare_events(rows)
    result={'version':'cirsium_flower_partner_regional_turnover_v1',
      'computed_date':'2026-10-20',
      'source':{'repository':'zuizui0223/aza3','workflow_run_id':37743990447,
        'artifact_id':11534639952,'snapshot_date':'2026-10-08',
        'archive_sha256':ARCHIVE_SHA256,'observation_csv_sha256':expected,
        'GloBI_citation_doi':'10.5281/zenodo.3950589'},
      'source_admission':audit,
      'design':{'primary_grid_degrees':2,'primary_records_per_unit':10,
        'min_shared_regions_per_host_pair':3,
        'metric':'Exact expected Bray-Curtis after equal-size record rarefaction',
        'contrast':'Same hosts across regions minus different hosts within those same regions',
        'primary_aggregation':'Mean within each host pair, then equal weight across host pairs',
        'analysis_status':'EXPLORATORY_SOURCE_RECORD_PATTERN_NOT_PREREGISTERED_CONFIRMATORY_TEST'},
      'scenarios':{},
      'claim_limits':[
        'Flower visit records are not measured pollination effectiveness.',
        'Plant and visitor names are original binomial strings, not new taxonomic validation.',
        'All retained event dates must be reported; missing dates prevent contemporaneous-community claims.',
        'Citizen observations lack measured effort, observer coverage and nondetection controls.',
        'Regions use degree grids with unequal physical area; grain sensitivity is descriptive.',
        'Rectangles and host pairs reuse species and regions; they are NOT independent experiments.',
        'Rarefaction balances number of recorded events, not biological observation effort.',
        'No causal geographic driver, head-trait benefit, selection, field permission or re-expression is inferred.'
      ]}
    args.out_dir.mkdir(parents=True,exist_ok=True)
    for grid,r,offset,label in [(2,10,0,'primary'),(1,10,0,'grid_1deg'),
            (5,10,0,'grid_5deg'),(2,5,0,'records_5'),(2,20,0,'records_20'),
            (2,10,1,'longitude_shift_1deg')]:
        s,rect,matrix=analyse_grid(events,grid,r,offset)
        result['scenarios'][label]=s
        if label=='primary':
            rect.to_csv(args.out_dir/'matched_rectangles.csv',index=False)
            matrix.to_csv(args.out_dir/'host_region_visitor_counts.csv',index=False)
    (args.out_dir/'result.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'admission':audit,'results':{k:{
       'host_pairs':v['n_host_pairs'],'rectangles':v['n_reused_rectangles'],
       'delta':v['equal_host_pair_mean_delta'],
       'positive_pairs':v.get('n_host_pairs_positive'),
       'without_common_hosts_delta':v['excluding_arvense_and_vulgare']['equal_host_pair_mean_delta']}
       for k,v in result['scenarios'].items()}},indent=2))

if __name__=='__main__':main()
