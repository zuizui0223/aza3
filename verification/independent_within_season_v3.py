#!/usr/bin/env python3
"""Independent record-pair enumeration and combinatorial reference verification.

Imports the previous frozen source admission, but NONE of the v3 statistics.
"""
import argparse
from collections import Counter,defaultdict
import hashlib,io,itertools,json,math,sys
from pathlib import Path
from zipfile import ZipFile
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'analysis'))
from cirsium_flower_partner_regional_turnover_v1 import prepare_events


def explicit(a,b):
    aa=sum(x==y for i,x in enumerate(a) for j,y in enumerate(a) if i!=j)/(len(a)*(len(a)-1))
    bb=sum(x==y for i,x in enumerate(b) for j,y in enumerate(b) if i!=j)/(len(b)*(len(b)-1))
    ab=sum(x==y for x in a for y in b)/(len(a)*len(b))
    return (aa+bb)/2-ab


def combinatorial_reference(g,a,b,keys):
    cats=sorted(g.partner_taxon.unique());means=Counter();variances=Counter()
    for _,z in g.groupby(keys,sort=True):
        counts=Counter(z.partner_taxon);n=len(z);take=sum(z.focal_taxon==a)
        denom=math.comb(n,take)
        for species,k in counts.items():
            mu=0.;second=0.
            for x in range(max(0,take-(n-k)),min(take,k)+1):
                p=math.comb(k,x)*math.comb(n-k,take-x)/denom
                mu+=p*x;second+=p*x*x
            means[species]+=mu;variances[species]+=second-mu*mu
    n_a=sum(g.focal_taxon==a);n_b=sum(g.focal_taxon==b);pool=Counter(g.partner_taxon)
    same_a=same_b=cross=0.
    for species in cats:
        ea=means[species];eaa=ea*ea+variances[species];k=pool[species]
        same_a+=(eaa-ea)/(n_a*(n_a-1))
        same_b+=(k*k-2*k*ea+eaa-k+ea)/(n_b*(n_b-1))
        cross+=(k*ea-eaa)/(n_a*n_b)
    return (same_a+same_b)/2-cross


def region_mean(items):
    d=defaultdict(list)
    for cell,v in items:d[cell].append(v)
    return sum(sum(v)/len(v) for v in d.values())/len(d)


def close(a,b):
    assert math.isclose(a,b,rel_tol=0,abs_tol=1e-12),(a,b)


def main():
    p=argparse.ArgumentParser();p.add_argument('--archive',type=Path,required=True)
    p.add_argument('--metadata',type=Path,required=True);p.add_argument('--result',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);args=p.parse_args()
    r=json.loads(args.result.read_text());source=r['source']
    assert hashlib.sha256(args.archive.read_bytes()).hexdigest()==source['archive_sha256']
    assert hashlib.sha256(args.metadata.read_bytes()).hexdigest()==source['metadata_sha256']
    with ZipFile(args.archive) as z:raw=z.read('cirsium_globi_observations.csv')
    assert hashlib.sha256(raw).hexdigest()==source['csv_sha256']
    d,admission=prepare_events(pd.read_csv(io.BytesIO(raw),dtype=str,keep_default_na=False))
    m=json.loads(args.metadata.read_text())
    dates=pd.DataFrame([{'event_id':o['event_id'],'date':o['observed_on'],
        'observer_key':o['observer_key']} for o in m['observations']])
    d=d.merge(dates,on='event_id',validate='one_to_one');assert len(d)==3243
    d['month']=d.date.str[5:7];d['year']=d.date.str[:4]
    d['cell']=d.apply(lambda x:f'{min(math.floor((x.lat+90)/2),89)}_{math.floor(((x.lon+180)%360)/2)}',axis=1)
    new_units=defaultdict(list)
    for (cell,month),g in d.groupby(['cell','month']):
        n=g.groupby('focal_taxon').size()
        for a,b in itertools.combinations(sorted(n[n>=10].index),2):
            new_units[(a,b)].append((cell,month))
    primary=r['scenario_results']['primary_month10']
    assert len(new_units)==primary['n_host_pairs']==4
    checked=0
    for q in primary['host_pairs']:
        a,b=q['host_a'],q['host_b'];expected=new_units[(a,b)]
        assert sorted(expected)==sorted((s['cell'],s['period']) for s in q['by_stratum'])
        observed=[];reference=[]
        for s in q['by_stratum']:
            g=d[(d.cell==s['cell'])&(d.month==s['period'])&d.focal_taxon.isin([a,b])]
            x=g[g.focal_taxon==a].partner_taxon.tolist();y=g[g.focal_taxon==b].partner_taxon.tolist()
            u=explicit(x,y);ref=combinatorial_reference(g,a,b,['observer_key'])
            close(u,s['observed_u']);close(ref,s['observer_reference_mean_u'])
            observed.append((s['cell'],u));reference.append((s['cell'],ref));checked+=1
        close(region_mean(observed),q['observed_u']);close(region_mean(reference),q['observer_reference_mean_u'])
    extra_count=0
    for extra in r['post_primary_exploratory_time_block_references']:
        a,b=extra['host_a'],extra['host_b']
        for name,ref in extra['references'].items():
            means=[]
            for s in ref['by_stratum']:
                g=d[(d.cell==s['cell'])&(d.month==s['period'])&d.focal_taxon.isin([a,b])]
                expected=combinatorial_reference(g,a,b,ref['block_keys'])
                close(expected,s['reference_mean_u']);means.append((s['cell'],expected))
                extra_count+=1
            close(region_mean(means),ref['exact_reference_mean_u'])
    # Aggregate files must not publish personal identifiers or point locations.
    def inspect(x):
        if isinstance(x,dict):
            assert not ({'lat','lon','latitude','longitude','event_id','observer_key'}&set(x))
            for y in x.values():inspect(y)
        elif isinstance(x,list):
            for y in x:inspect(y)
    inspect(r)
    out={'status':'PASS_INDEPENDENT_EXPLICIT_RECORD_PAIRS_AND_COMBINATORIAL_HYPERGEOMETRIC_MEANS',
      'source_records':len(d),'primary_host_pairs_verified':4,'primary_strata_verified':checked,
      'post_primary_reference_strata_verified':extra_count,'new_biological_observations':0,
      'output_contains_precise_coordinates_or_observer_identifiers':False,
      'result_sha256':hashlib.sha256(args.result.read_bytes()).hexdigest()}
    args.out.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))


if __name__=='__main__':main()
