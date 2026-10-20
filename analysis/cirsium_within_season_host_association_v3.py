#!/usr/bin/env python3
"""Fixed-cohort Cirsium visitor identity association within region/calendar strata.

U is the mean within-host identity-match probability minus the cross-host
match probability, with no self-pairs. Half equality-kernel unbiased MMD^2.
Nulls reassign recorded labels; they are NOT randomized ecological experiments.
Source observation effort, local floral availability, and morphology are absent.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
import hashlib
import io
import itertools
import json
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd
from cirsium_flower_partner_regional_turnover_v1 import ARCHIVE_SHA256, OBSERVATION_MEMBER, prepare_events
from cirsium_visit_dates_v2 import join_dates, add_cells, reduce_observer_days

METADATA_SHA256='313142b1e4fcb1a179533c9b985ea8364b221122ed158b9e2d312ee213a731c9'
SOURCE_MAIN='cd94a991dcc86b5e109f4def22efc7755eaa22fb'
DESIGN_COMMIT='7eea08a9d067bf3217a1dd0a9aa2ec4c598391c3'
ORDERS={'Lepidoptera','Hymenoptera','Diptera','Coleoptera','Hemiptera','Thysanoptera','Orthoptera'}


def _counts(a,b,minimum=2):
    a,b=np.asarray(a),np.asarray(b)
    if (a.ndim!=1 or b.shape!=a.shape or len(a)==0 or
        not np.isfinite(a).all() or not np.isfinite(b).all() or
        (a<0).any() or (b<0).any() or (a!=np.floor(a)).any() or
        (b!=np.floor(b)).any() or a.sum()<minimum or b.sum()<minimum):
        raise ValueError('INVALID_ALIGNED_NONNEGATIVE_INTEGER_COUNTS')
    return a.astype(np.int64),b.astype(np.int64)


def identity_contrast(a,b):
    """May be negative; under iid sampling target is .5*sum((pA-pB)^2)."""
    a,b=_counts(a,b)
    na,nb=a.sum(),b.sum()
    return float(.5*((a*(a-1)).sum()/(na*(na-1))+
                      (b*(b-1)).sum()/(nb*(nb-1)))-(a*b).sum()/(na*nb))


def _blocks(blocks):
    if not blocks:raise ValueError('EMPTY_REFERENCE_BLOCKS')
    normalized=[_counts(a,b,0) for a,b in blocks]
    k=len(normalized[0][0])
    if any(len(a)!=k or a.sum()+b.sum()==0 for a,b in normalized):
        raise ValueError('INCOMPATIBLE_OR_EMPTY_REFERENCE_BLOCK')
    a=sum(x for x,y in normalized);b=sum(y for x,y in normalized)
    _counts(a,b)
    return normalized,a,b


def block_reference_expectation(blocks):
    """Exact E[U] under independent fixed-margin hypergeometric label blocks.

    E[A_s^2]=E[A_s]^2+Var(A_s). The U polynomial needs sum of these
    marginal second moments, not cross-category covariances.
    """
    blocks,a,b=_blocks(blocks)
    na,nb=int(a.sum()),int(b.sum());total=a+b
    mu=np.zeros(len(a));var=np.zeros(len(a))
    for x,y in blocks:
        k=x+y;n=int(k.sum());take=int(x.sum());p=k/n
        mu+=take*p
        if n>1:var+=take*p*(1-p)*(n-take)/(n-1)
    second=float((mu*mu+var).sum());linear=float(total@mu)
    ea=(second-na)/(na*(na-1))
    eb=(float(total@total)-2*linear+second-nb)/(nb*(nb-1))
    cross=(linear-second)/(na*nb)
    return float(.5*(ea+eb)-cross)


def reference_draws(blocks,reps,seed):
    """Fixed-margin host-label reassignment; no new biological data."""
    if not isinstance(reps,int) or reps<1:raise ValueError('INVALID_REPLICATES')
    blocks,a,b=_blocks(blocks);total=a+b;na,nb=int(a.sum()),int(b.sum())
    rng=np.random.default_rng(seed);draw=np.zeros((reps,len(a)),dtype=np.int64)
    for x,y in blocks:
        k=x+y;take=int(x.sum());n=int(k.sum())
        if take==n:draw+=k
        elif take>0:draw+=rng.multivariate_hypergeometric(k,take,size=reps)
    other=total-draw
    return .5*((draw*(draw-1)).sum(1)/(na*(na-1))+
                (other*(other-1)).sum(1)/(nb*(nb-1)))-(draw*other).sum(1)/(na*nb)


def region_balanced_mean(values):
    if not values:raise ValueError('EMPTY_REGION_MEAN')
    groups=defaultdict(list)
    for region,v in values:groups[region].append(v)
    return float(np.mean([np.mean(v) for v in groups.values()]))


def _observer_units(d):
    key=d.observer_key.fillna('').astype(str)
    return key.where(key.ne(''),'missing-event-'+d.event_id.astype(str))


def _unit(group,host_a,host_b,visitor='partner_taxon'):
    categories=sorted(group[visitor].unique())
    def counts(z,host):
        return z.loc[z.focal_taxon.eq(host),visitor].value_counts().reindex(categories,fill_value=0).to_numpy(dtype=np.int64)
    a,b=counts(group,host_a),counts(group,host_b)
    _counts(a,b)
    blocks=[];movable=0;exchangeable=0
    for _,z in group.groupby(_observer_units(group),sort=True):
        x,y=counts(z,host_a),counts(z,host_b);blocks.append((x,y))
        if x.sum()>0 and y.sum()>0:movable+=len(z);exchangeable+=1
    na,nb=int(a.sum()),int(b.sum())
    wa=float((a*(a-1)).sum()/(na*(na-1)))
    wb=float((b*(b-1)).sum()/(nb*(nb-1)))
    ca=float((a*b).sum()/(na*nb))
    year_a=set(group.loc[group.focal_taxon.eq(host_a),'year'])
    year_b=set(group.loc[group.focal_taxon.eq(host_b),'year'])
    details={'n_host_a':na,'n_host_b':nb,'n_visitor_names':len(categories),
       'within_host_a_match':wa,'within_host_b_match':wb,'mean_within_host_match':(wa+wb)/2,
       'cross_host_match':ca,'observed_u':identity_contrast(a,b),
       'observer_reference_mean_u':block_reference_expectation(blocks),
       'unrestricted_reference_mean_u':block_reference_expectation([(a,b)]),
       'n_observers':int(group.loc[group.observer_key.ne(''),'observer_key'].nunique()),
       'n_exchangeable_observer_blocks':exchangeable,'n_exchangeable_observer_block_records':movable,
       'n_calendar_years':len(year_a|year_b),'n_shared_calendar_years':len(year_a&year_b),
       'bray_curtis_raw_proportions':float(np.abs(a/na-b/nb).sum()/2)}
    details['observer_excess_u']=details['observed_u']-details['observer_reference_mean_u']
    return details,blocks,a,b


def _reference_summary(draws,observed,analytic):
    return {'replicates':len(draws),'exact_mean_u':analytic,
        'monte_carlo_mean_u':float(np.mean(draws)),
        'reference_quantiles_025_975':np.quantile(draws,[.025,.975]).tolist(),
        'upper_tail_label_reference_probability':float((1+np.sum(draws>=observed-1e-12))/(len(draws)+1)),
        'is_biological_confidence_interval':False,
        'is_field_randomization_p_value':False}


def analyse(events,stratum,min_n,permutations=9999,visitor='partner_taxon',diagnostics=False):
    required={'event_id','focal_taxon',visitor,'cell',stratum,'year','observer_key'}
    if not required.issubset(events.columns):raise ValueError('MISSING_ANALYSIS_COLUMNS')
    if events.event_id.duplicated().any():raise ValueError('DUPLICATE_SOURCE_EVENT')
    if min_n<2:raise ValueError('MINIMUM_TWO_FOR_IDENTITY_COMPARISON')
    pairs=defaultdict(list)
    for (cell,period),g in events.groupby(['cell',stratum],sort=True):
        counts=g.focal_taxon.value_counts()
        hosts=sorted(counts.index[counts>=min_n])
        for a,b in itertools.combinations(hosts,2):
            z=g[g.focal_taxon.isin([a,b])].copy()
            pairs[(a,b)].append((str(cell),str(period),z))
    results=[]
    for (a,b),units in sorted(pairs.items()):
        stats=[];memo=[]
        for cell,period,z in units:
            s,blocks,aa,bb=_unit(z,a,b,visitor)
            s.update(cell=cell,period=period);stats.append(s);memo.append((blocks,aa,bb))
        cells=sorted({s['cell'] for s in stats})
        mean=lambda key:region_balanced_mean([(s['cell'],s[key]) for s in stats])
        all_records=pd.concat([z for cell,period,z in units],ignore_index=True)
        # Every event belongs to one stratum within this pair, not across pairs.
        if all_records.event_id.duplicated().any():raise ValueError('PAIR_EVENT_REUSE_ACROSS_STRATA')
        q={'host_a':a,'host_b':b,'n_regions':len(cells),'n_strata':len(stats),
           'n_unique_records':len(all_records),
           'n_unique_observers':int(all_records.loc[all_records.observer_key.ne(''),'observer_key'].nunique()),
           'status':'REPLICATED_THREE_OR_MORE_REGIONS' if len(cells)>=3 else 'EXPLORATORY_FEWER_THAN_THREE_REGIONS',
           'observed_u':mean('observed_u'),'observer_reference_mean_u':mean('observer_reference_mean_u'),
           'observer_excess_u':mean('observer_excess_u'),
           'mean_within_host_match':mean('mean_within_host_match'),'cross_host_match':mean('cross_host_match'),
           'n_exchangeable_observer_block_records':sum(s['n_exchangeable_observer_block_records'] for s in stats),
           'n_exchangeable_observer_blocks':sum(s['n_exchangeable_observer_blocks'] for s in stats),
           'by_stratum':stats,'leave_one_region_out':{},'label_reference':None}
        for cell in cells:
            keep=[s for s in stats if s['cell']!=cell]
            q['leave_one_region_out'][cell]={key:region_balanced_mean([(s['cell'],s[key]) for s in keep]) if keep else None
                for key in ['observed_u','observer_excess_u']}
        if len(cells)>=3 and permutations:
            weighted={'unrestricted':np.zeros(permutations),'observer_conditioned':np.zeros(permutations)}
            for i,(s,(blocks,aa,bb)) in enumerate(zip(stats,memo)):
                weight=1/len(cells)/sum(t['cell']==s['cell'] for t in stats)
                stable=int(hashlib.sha256(f'{a}|{b}|{s["cell"]}|{s["period"]}|{visitor}'.encode()).hexdigest()[:8],16)
                weighted['unrestricted']+=weight*reference_draws([(aa,bb)],permutations,stable)
                weighted['observer_conditioned']+=weight*reference_draws(blocks,permutations,stable+1)
            q['label_reference']={kind:_reference_summary(v,q['observed_u'],0. if kind=='unrestricted' else q['observer_reference_mean_u'])
                for kind,v in weighted.items()}
        if diagnostics and len(cells)>=3:
            vals=[];invalid=0;under_threshold=0
            keys=sorted(set(_observer_units(all_records)))
            for key in keys:
                v=[];vi=[];bad=False
                for cell,period,z in units:
                    reduced=z[_observer_units(z)!=key]
                    ns=reduced.focal_taxon.value_counts().reindex([a,b],fill_value=0)
                    if (ns<2).any():bad=True;break
                    if (ns<min_n).any():under_threshold+=1
                    s,_,_,_=_unit(reduced,a,b,visitor)
                    v.append((cell,s['observed_u']));vi.append((cell,s['observer_excess_u']))
                if bad:invalid+=1
                else:vals.append((region_balanced_mean(v),region_balanced_mean(vi)))
            q['leave_one_observer_out_frozen_strata']={
              'n_deletions':len(keys),'n_uncomputable_deletions':invalid,
              'n_stratum_deletions_below_original_minimum':under_threshold,
              'observed_u_range':[min(v[0] for v in vals),max(v[0] for v in vals)] if vals else None,
              'observer_excess_u_range':[min(v[1] for v in vals),max(v[1] for v in vals)] if vals else None,
              'changes_support_threshold':'Original strata retained; n>=2 needed for U after deletion. Not a new primary cohort.'}
        results.append(q)
    replicated=[q for q in results if q['n_regions']>=3]
    return {'stratum':stratum,'minimum_records_per_host_stratum':min_n,'visitor_resolution':visitor,
        'input_records':len(events),'n_host_pairs':len(results),'n_replicated_host_pairs':len(replicated),
        'equal_replicated_pair_mean_u':float(np.mean([q['observed_u'] for q in replicated])) if replicated else None,
        'host_pairs':results,'interpretation':'SOURCE_RECORD_ASSOCIATION_NOT_MORPHOLOGICAL_FILTERING_OR_PREFERENCE'}


def additional_reference_for_pair(events,pair,block_keys,reps=9999):
    """Post-primary diagnostic on FIXED primary units, not a new matched cohort.

    This was added after observing primary v3 results, so is explicitly
    exploratory. Date blocking cannot create concurrent botanical surveys.
    """
    a,b=pair['host_a'],pair['host_b'];details=[];draws=np.zeros(reps)
    cells={s['cell'] for s in pair['by_stratum']};movable=0;nblocks=0
    for s in pair['by_stratum']:
        z=events[events.cell.eq(s['cell']) & events.month.eq(s['period']) &
                 events.focal_taxon.isin([a,b])].copy()
        if len(z)!=s['n_host_a']+s['n_host_b']:raise ValueError('FROZEN_PRIMARY_UNIT_DRIFT')
        z['observer_key']=_observer_units(z)
        cats=sorted(z.partner_taxon.unique());blocks=[];unit_move=0;unit_blocks=0
        for _,g in z.groupby(block_keys,sort=True):
            x=g.loc[g.focal_taxon.eq(a),'partner_taxon'].value_counts().reindex(cats,fill_value=0).to_numpy()
            y=g.loc[g.focal_taxon.eq(b),'partner_taxon'].value_counts().reindex(cats,fill_value=0).to_numpy()
            blocks.append((x,y))
            if x.sum()>0 and y.sum()>0:unit_move+=len(g);unit_blocks+=1
        ref=block_reference_expectation(blocks)
        details.append({'cell':s['cell'],'period':s['period'],
            'observed_u':s['observed_u'],'reference_mean_u':ref,
            'excess_u':s['observed_u']-ref,'n_exchangeable_records':unit_move,
            'n_exchangeable_label_blocks':unit_blocks})
        movable+=unit_move;nblocks+=unit_blocks
        weight=1/len(cells)/sum(t['cell']==s['cell'] for t in pair['by_stratum'])
        seed=int(hashlib.sha256(f'{a}|{b}|{s["cell"]}|{s["period"]}|{block_keys}'.encode()).hexdigest()[:8],16)
        draws+=weight*reference_draws(blocks,reps,seed)
    avg=region_balanced_mean([(s['cell'],s['reference_mean_u']) for s in details])
    return {'block_keys':block_keys,'post_primary_exploratory':True,
        'n_total_primary_records':pair['n_unique_records'],'n_exchangeable_records':movable,
        'n_exchangeable_label_blocks':nblocks,'exact_reference_mean_u':avg,
        'observed_minus_reference_u':pair['observed_u']-avg,'by_stratum':details,
        'leave_one_region_out_excess':{cell:region_balanced_mean([(s['cell'],s['excess_u']) for s in details if s['cell']!=cell])
            if len(cells)>1 else None for cell in sorted(cells)},
        'reference':_reference_summary(draws,pair['observed_u'],avg)}


def load_sources(archive,metadata):
    if hashlib.sha256(archive.read_bytes()).hexdigest()!=ARCHIVE_SHA256:raise ValueError('ARCHIVE_HASH_MISMATCH')
    if hashlib.sha256(metadata.read_bytes()).hexdigest()!=METADATA_SHA256:raise ValueError('METADATA_HASH_MISMATCH')
    with ZipFile(archive) as z:
        raw=z.read(OBSERVATION_MEMBER);manifest=json.loads(z.read('summary.json'))
    if hashlib.sha256(raw).hexdigest()!=manifest['csv_sha256']:raise ValueError('CSV_HASH_MISMATCH')
    events,admission=prepare_events(pd.read_csv(io.BytesIO(raw),dtype=str,keep_default_na=False))
    meta=json.loads(metadata.read_text());dated,recovery=join_dates(events,meta['observations'])
    if len(dated)!=3243 or len(events)!=3243:raise ValueError('FIXED_3243_COHORT_DRIFT')
    return add_cells(dated),{'archive_sha256':ARCHIVE_SHA256,'metadata_sha256':METADATA_SHA256,
        'csv_sha256':manifest['csv_sha256'],'admission':admission,'date_recovery':recovery}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--archive',type=Path,required=True);p.add_argument('--metadata',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--permutations',type=int,default=9999)
    args=p.parse_args()
    if not 99<=args.permutations<=99999:raise ValueError('PERMUTATION_BUDGET_OUT_OF_BOUNDS')
    d,source=load_sources(args.archive,args.metadata)
    concordant=d.apply(lambda r:r.partner_taxon==r.current_taxon_name or r.partner_taxon in r.current_ancestor_names,axis=1)
    scenarios={}
    for name,frame,time,minimum in [
        ('primary_month10',d,'month',10),
        ('month10_observer_day_dedup',reduce_observer_days(d),'month',10),
        ('month10_current_visitor_concordant',d[concordant],'month',10),
        ('month5_exploratory',d,'month',5),
        ('strict_year_month10',d,'year_month',10),
        ('strict_year_month5',d,'year_month',5)]:
        scenarios[name]=analyse(frame,time,minimum,permutations=args.permutations if name=='primary_month10' else 0,
            diagnostics=name=='primary_month10')
    def get_order(path):
        found={s.strip() for s in path.split('|')}&ORDERS
        return next(iter(found)) if len(found)==1 else 'UNRESOLVED'
    d['visitor_order']=d.source_taxon_path.map(get_order)
    scenarios['month10_insect_order_same_support']=analyse(d,'month',10,permutations=0,visitor='visitor_order')
    extra=[]
    for pair in scenarios['primary_month10']['host_pairs']:
        if pair['n_regions']>=3:
            extra.append({'host_a':pair['host_a'],'host_b':pair['host_b'],
                'references':{name:additional_reference_for_pair(d,pair,keys,args.permutations)
                 for name,keys in [('year',['year']),('observer_year',['observer_key','year']),
                                   ('observer_date',['observer_key','date'])]}})
    result={'version':'cirsium_within_season_host_association_v3','computed_date':'2026-10-20',
        'source_main':SOURCE_MAIN,'design_commit':DESIGN_COMMIT,'source':source,
        'estimand':'Within same cell/time stratum, mean within-host visitor identity match minus across-host match; no self-pairs. Equal cells then equal periods within cells.',
        'observer_reference':'Exact expectation under observer-specific host totals and pooled visitor counts inside each region/time stratum.',
        'scenario_results':scenarios,'post_primary_exploratory_time_block_references':extra,
        'claim_limits':[
            'This changes the estimand from geographic dominance to within-stratum host association; original empty time-matched rectangles remain empty.',
            'Calendar month pools years; year-month controls still do not measure common sampling effort or local flower abundance.',
            'Host labels were not randomized. Label-reference probabilities and ranges are conditional diagnostics, not field causal tests or population confidence intervals.',
            'Observer conditioning holds allocations fixed; its residual is not a causal observer correction. Single-host observers cannot exchange labels.',
            'U is a finite-sample signed estimator; negative estimates are not negative biological distances and are not truncated.',
            'Source names and interactions are not new taxonomic or pollination-efficacy validation.',
            'Repeated observations, unmeasured within-cell locations, dates, detection, and effort can confound host association.',
            'No head structure, floral preferences conditional on availability, viable achene output, selection, genetic reuse, or field permissions are measured.'
        ]}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    for name,s in scenarios.items():
        print(name,'pairs',s['n_host_pairs'],'replicated',s['n_replicated_host_pairs'],flush=True)
        for q in s['host_pairs']:
            print(q['host_a'],q['host_b'],q['n_regions'],q['n_unique_records'],q['observed_u'],q['observer_reference_mean_u'],q['observer_excess_u'],flush=True)


if __name__=='__main__':main()
