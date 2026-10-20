#!/usr/bin/env python3
"""Exploratory v4 audit of GloBI Cirsium host/visitor-record association.

New question (not the v3 predeclared estimand): can the v3 same-region,
same-calendar-month 191-observation apparent visitor-name association be
explained only by recorded insect ORDER and observer-specific host records?

The observation unit is an iNaturalist photograph relationship, NOT a bee
visit rate, physical head access, head phenotype, site replicate or fitness.
Source IDs, dates and observer hashes are used transiently and not persisted
in result JSON. Diagnostic conditional label shuffles are NOT ecological
randomization, a selectivity preference test or genotype evidence.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import pandas as pd
from cirsium_within_season_host_association_v3 import (
    load_sources, identity_contrast, block_reference_expectation,
    reference_draws, _unit,
)

HOST_A='Cirsium arvense'
HOST_B='Cirsium vulgare'
INSECT_ORDERS={'Lepidoptera','Hymenoptera','Diptera','Coleoptera',
               'Hemiptera','Thysanoptera','Orthoptera'}
EXPECTED_COARSE_ORDERS={'Lepidoptera','Hymenoptera','Diptera','Coleoptera'}
REDACT={'lat','lon','latitude','longitude','event_id','observer_key',
        'source_specimen_occurrence_id','source_specimen_life_stage',
        'observed_on','date','locality','geoprivacy','user','username',
        'site_precise_coords','observer_hashes','photos'}


def source_insect_order(path:str)->str:
    found=set(t.strip() for t in path.split('|')) & INSECT_ORDERS
    if len(found)!=1:raise ValueError('ORDER_AMBIGUOUS_OR_MISSING_IN_SOURCE_PATH')
    return next(iter(found))


def _checked_block_arrays(blocks):
    if not blocks:raise ValueError('NO_SOURCE_OBSERVER_BLOCKS')
    k=len(blocks[0][0]);clean=[]
    for a,b in blocks:
        a=np.asarray(a);b=np.asarray(b)
        if (a.ndim!=1 or b.ndim!=1 or a.shape!=(k,) or b.shape!=(k,)
            or not np.isfinite(a).all() or not np.isfinite(b).all()
            or (a<0).any() or (b<0).any()
            or not np.equal(a,np.floor(a)).all()
            or not np.equal(b,np.floor(b)).all()
            or a.sum()+b.sum()==0):
            raise ValueError('BAD_FIXED_SOURCE_COUNTS')
        clean.append((a.astype(np.int64),b.astype(np.int64)))
    return clean


def exact_category_moments(blocks):
    """E[A_s], Var[A_s] under observer(-order) fixed-margin shuffles."""
    blocks=_checked_block_arrays(blocks)
    k=len(blocks[0][0]);mean=np.zeros(k);var=np.zeros(k)
    for a,b in blocks:
        totals=a+b;N=int(totals.sum());take=int(a.sum());p=totals/N
        mean+=take*p
        if N>1:var+=take*p*(1-p)*(N-take)/(N-1)
    return mean,var


def category_identity_components(blocks):
    """Observed, exact-reference, excess contribution of EACH visitor species.

    Sum over species is unbiased equality-kernel U; signed contributions may
    be negative and MUST NOT be clamped to zero. Analytic expectation uses
    hypergeometric second moments within all provided observer blocks.
    """
    blocks=_checked_block_arrays(blocks)
    a=sum((x for x,y in blocks));b=sum((y for x,y in blocks))
    na,nb=int(a.sum()),int(b.sum())
    if min(na,nb)<2:raise ValueError('MINIMUM_TWO_PER_HOST_FOR_U')
    t=a+b;mu,var=exact_category_moments(blocks)
    second=mu*mu+var
    obs=.5*(a*(a-1)/(na*(na-1))+b*(b-1)/(nb*(nb-1)))-a*b/(na*nb)
    expect=.5*((second-mu)/(na*(na-1))+
             ((t-mu)**2+var-(t-mu))/(nb*(nb-1))) -(t*mu-second)/(na*nb)
    if not np.isclose(obs.sum(),identity_contrast(a,b),atol=1e-12,rtol=0):
        raise ValueError('OBSERVED_SPECIES_PARTITION_INCONSISTENT')
    if not np.isclose(expect.sum(),block_reference_expectation(blocks),atol=1e-12,rtol=0):
        raise ValueError('EXPECTED_SPECIES_PARTITION_INCONSISTENT')
    return obs,expect,obs-expect


def block_vectors(g,host_a,host_b,visitor_col,keys):
    taxa=sorted(g[visitor_col].unique())
    if len(taxa)<1:raise ValueError('NO_VISITOR_NAMES_IN_GLOBI_STRATUM')
    if not isinstance(keys,(tuple,list)) or not keys:raise ValueError('MISSING_SOURCE_BLOCK_KEYS')
    if not set(keys).issubset(set(g.columns)):raise ValueError('MISSING_LABEL_REFERENCE_KEYS')
    blocks=[];exchangeable_n=0;exchangeable_blocks=0
    for _,z in g.groupby(list(keys),sort=True,dropna=False):
        aa=z.loc[z.focal_taxon.eq(host_a),visitor_col].value_counts().reindex(taxa,fill_value=0).to_numpy()
        bb=z.loc[z.focal_taxon.eq(host_b),visitor_col].value_counts().reindex(taxa,fill_value=0).to_numpy()
        blocks.append((aa,bb))
        if aa.sum()>0 and bb.sum()>0:
            exchangeable_n+=len(z);exchangeable_blocks+=1
    return taxa,blocks,exchangeable_n,exchangeable_blocks


def order_conditioned_units(g,host_a,host_b,keys=('observer_key','visitor_order'),
                            visitor_col='partner_taxon'):
    taxa,blocks,exchangeable_n,exchangeable_blocks=block_vectors(
        g,host_a,host_b,visitor_col,keys)
    a=sum((x for x,y in blocks));b=sum((y for x,y in blocks))
    observed=identity_contrast(a,b)
    expected=block_reference_expectation(blocks)
    return {'n_records':len(g),'n_visitor_names':len(taxa),
        'n_host_a':int(a.sum()),'n_host_b':int(b.sum()),
        'observed_u':observed,'label_reference_mean_u':expected,
        'reference_excess_u':observed-expected,
        # Alias deliberately kept only to support machine checks for this v4.
        'observer_order_excess_u':observed-expected,
        'n_exchangeable_records':exchangeable_n,
        'n_exchangeable_blocks':exchangeable_blocks}


def observer_fixed_order_difference(g,host_a,host_b,target_order='Hymenoptera'):
    """Difference in proportion of records at target order, B minus A.

    Analytic expectation holds each observer's host record totals and order
    distribution fixed within each actual region-month. No survey effort
    or floral availability interpretation is allowed.
    """
    aa=g[g.focal_taxon.eq(host_a)];bb=g[g.focal_taxon.eq(host_b)]
    na,nb=len(aa),len(bb)
    if min(na,nb)<2:raise ValueError('INVALID_GROUP_FOR_ORDER_DIFFERENCE')
    order=g.visitor_order.eq(target_order)
    total=int(order.sum());obs=float(bb.visitor_order.eq(target_order).mean()-aa.visitor_order.eq(target_order).mean())
    expected_a=0.; var_a=0.;exchange_n=0
    for _,z in g.groupby('observer_key',sort=True,dropna=False):
        N=len(z);a=int(z.focal_taxon.eq(host_a).sum());b=N-a
        H=int(z.visitor_order.eq(target_order).sum());p=H/N
        expected_a+=a*p
        if N>1:var_a+=a*p*(1-p)*b/(N-1)
        if a>0 and b>0:exchange_n+=N
    expected=(total-expected_a)/nb-expected_a/na
    return {'n_host_a':na,'n_host_b':nb,'n_target_order_a':int(aa.visitor_order.eq(target_order).sum()),
        'n_target_order_b':int(bb.visitor_order.eq(target_order).sum()),
        'observed_host_b_minus_a':obs,'observer_reference_difference':expected,
        'observer_excess_difference':obs-expected,
        'observer_reference_variance':var_a*(1/na+1/nb)**2,
        'n_exchangeable_records':exchange_n}


def frozen_cohort(events,expected):
    """Reconstruct exactly the BEFORE-v4 v3 outcome-independent source support."""
    required={'event_id','cell','month','focal_taxon'}
    if not required.issubset(events.columns):raise ValueError('SOURCE_PRIMARY_COLUMNS_MISSING')
    if events.event_id.duplicated().any():raise ValueError('SOURCE_PRIMARY_DUPLICATE_ID')
    hosts={HOST_A,HOST_B}
    if any((not isinstance(k,tuple) or len(k)!=2 or len(v)!=2 or min(v)<2)
           for k,v in expected.items()):raise ValueError('INVALID_FROZEN_PRIMARY_SUPPORT_CONTRACT')
    subset=events[events.focal_taxon.isin(hosts)].copy()
    subset=subset[subset.apply(lambda r:(r.cell,r.month) in expected,axis=1)]
    actual={(c,m):(int((g.focal_taxon==HOST_A).sum()),int((g.focal_taxon==HOST_B).sum()))
       for (c,m),g in subset.groupby(['cell','month'],sort=True)}
    if actual!=expected:raise ValueError('SOURCE_PRIMARY_COUNT_DRIFT')
    return subset


def summarize_audited_cohort(j):
    def walk(node):
        if isinstance(node,dict):
            for k,v in node.items():
                if k in REDACT:raise ValueError('REDACTED_FIELD_PRESENT:'+k)
                walk(v)
        elif isinstance(node,list):
            for v in node:walk(v)
    walk(j)
    return j


def primary_strata_from_previous_receipt(j):
    if (j.get('source_cohort')!=3243 or
        j.get('primary_host_pair')!=[HOST_A,HOST_B] or
        j['primary']['n_regions']!=5 or
        j['primary']['n_unique_records']!=191 or
        len(j['primary_units'])!=5):
        raise ValueError('FROZEN_V3_SUPPORT_OR_COHORT_DRIFT')
    expected={}
    for cell,month,na,nb,*extras in j['primary_units']:
        k=(cell,month)
        if k in expected:raise ValueError('DUPLICATED_PRIOR_SUPPORT_CELL')
        expected[k]=(int(na),int(nb))
    return expected


def _reference_prob(g,blocks,n,seed):
    return reference_draws(blocks,n,seed)


def analyse_source(archive:Path,metadata:Path,previous_receipt:Path,reps:int):
    if not 999<=reps<=99999:raise ValueError('INVALID_FIXED_REFERENCE_REPLICATION_BUDGET')
    old=json.loads(previous_receipt.read_text(encoding='utf-8'))
    events,source=load_sources(archive,metadata)
    if (old['source']['archive_sha256']!=source['archive_sha256'] or
        old['source']['metadata_sha256']!=source['metadata_sha256'] or
        old['source']['csv_sha256']!=source['csv_sha256']):
        raise ValueError('PREVIOUS_REAL_SOURCE_PROVENANCE_CHANGED')
    expected=primary_strata_from_previous_receipt(old)
    group=frozen_cohort(events,expected)
    if group.source_taxon_path.isna().any():raise ValueError('VISITOR_ORDER_MISSING')
    group['visitor_order']=group.source_taxon_path.map(source_insect_order)
    if set(group.visitor_order)!=EXPECTED_COARSE_ORDERS:
        raise ValueError('ORIGINAL_SOURCE_VISITOR_ORDER_COVERAGE_CHANGED')
    obs_by_host=group.groupby(['focal_taxon','visitor_order']).size()
    source_order_counts={h:{c:int(obs_by_host.get((h,c),0)) for c in sorted(EXPECTED_COARSE_ORDERS)}
                         for h in [HOST_A,HOST_B]}
    if sum(source_order_counts[HOST_A].values())!=107 or sum(source_order_counts[HOST_B].values())!=84:
        raise ValueError('FROZEN_PRIMARY_HOST_VISITOR_COUNTS_DRIFT')
    strata=[];type_contribution=defaultdict(float);name_contribution=defaultdict(float)
    sim=np.zeros(reps);within_lep=[];order_fractions={h:defaultdict(list) for h in [HOST_A,HOST_B]}
    for cell,mon in sorted(expected):
        z=group[group.cell.eq(cell)&group.month.eq(mon)].copy()
        prior=next(q for q in old['primary_units'] if q[0]==cell and q[1]==mon)
        base,_bb,_cc,_dd=_unit(z,HOST_A,HOST_B)
        if abs(base['observed_u']-prior[4])>1e-12 or abs(base['observer_reference_mean_u']-prior[5])>1e-12:
            raise ValueError('V3_ORIGINAL_STUDY_RESULT_NOT_REPRODUCED')
        obs=order_conditioned_units(z,HOST_A,HOST_B,('observer_key',))
        strict=order_conditioned_units(z,HOST_A,HOST_B,('observer_key','visitor_order'))
        org=order_conditioned_units(z,HOST_A,HOST_B,('observer_key','visitor_order','year'))
        day=order_conditioned_units(z,HOST_A,HOST_B,('observer_key','visitor_order','date'))
        bee=observer_fixed_order_difference(z,HOST_A,HOST_B,'Hymenoptera')
        if abs(obs['observer_order_excess_u']-base['observer_excess_u'])>1e-12:
            raise ValueError('ORIGINAL_OBSERVER_REFERENCE_CONFLICT')
        taxa,blocks,n_move,n_blocks=block_vectors(z,HOST_A,HOST_B,'partner_taxon',('observer_key','visitor_order'))
        o_c,r_c,d_c=category_identity_components(blocks)
        if not np.isclose(d_c.sum(),strict['observer_order_excess_u'],atol=1e-12):
            raise ValueError('SOURCE_VISITOR_DECOMPOSITION_FAILED')
        order_lookup=z.drop_duplicates('partner_taxon').set_index('partner_taxon').visitor_order.to_dict()
        if z.groupby('partner_taxon').visitor_order.nunique().gt(1).any():
            raise ValueError('ONE_VISITOR_HAS_CONFLICTING_INSECT_ORDERS')
        for name,value in zip(taxa,d_c):
            name_contribution[name]+=float(value)/len(expected)
            type_contribution[order_lookup[name]]+=float(value)/len(expected)
        lept=z[z.visitor_order.eq('Lepidoptera')].copy()
        n_lep=lept.focal_taxon.value_counts().reindex([HOST_A,HOST_B],fill_value=0)
        if min(n_lep)<5:raise ValueError('NOT_ENOUGH_LEPIDOPTERA_TO_CHECK_WITHIN_ORDER')
        lep=order_conditioned_units(lept,HOST_A,HOST_B,('observer_key',))
        within_lep.append(lep)
        for host in [HOST_A,HOST_B]:
            h=z[z.focal_taxon.eq(host)]
            for insect_order in EXPECTED_COARSE_ORDERS:
                order_fractions[host][insect_order].append(float(h.visitor_order.eq(insect_order).mean()))
        seed=int(hashlib.sha256(f'{cell}|{mon}|aza3_v4_observer_order'.encode()).hexdigest()[:8],16)
        sim+=_reference_prob(z,blocks,reps,seed)/len(expected)
        strata.append({'cell':cell,'month':mon,'n_records':len(z),'prior_observer_only':{
            'observed_u':obs['observed_u'],'reference_mean_u':obs['label_reference_mean_u'],
            'excess_u':obs['reference_excess_u']},
            'observer_x_order':strict,
            'observer_x_order_x_year':org,
            'observer_x_order_x_date':day,
            'hymenoptera_versus_other_observer_fixed':bee,
            'lepidoptera_only':lep})
    mean=lambda path:float(np.mean([path(q) for q in strata]))
    original_observer_excess=mean(lambda q:q['prior_observer_only']['excess_u'])
    strict_excess=mean(lambda q:q['observer_x_order']['reference_excess_u'])
    raw_u=mean(lambda q:q['observer_x_order']['observed_u'])
    strict_exp=mean(lambda q:q['observer_x_order']['label_reference_mean_u'])
    order_additional_reference=strict_exp-mean(lambda q:q['prior_observer_only']['reference_mean_u'])
    if abs(strict_excess-(original_observer_excess-order_additional_reference))>1e-12:
        raise ValueError('REFERENCE_ACCOUNTING_NOT_CONSERVED')
    if abs(sum(type_contribution.values())-strict_excess)>1e-12:
        raise ValueError('VISITOR_ORDER_CONTRIBUTIONS_NOT_ADDING_UP')
    cells=sorted(expected)
    leave_one_cell={cell:float(np.mean([q['observer_x_order']['reference_excess_u']
                              for q in strata if q['cell']!=cell])) for cell,month in cells}
    visitor_drop=[]
    for name in sorted(group.partner_taxon.unique()):
        z=group[group.partner_taxon.ne(name)]
        parts=[];eligible=True
        for cell,mon in cells:
            q=z[z.cell.eq(cell)&z.month.eq(mon)]
            if min(q.focal_taxon.value_counts().reindex([HOST_A,HOST_B],fill_value=0))<10:
                eligible=False
            parts.append(order_conditioned_units(q,HOST_A,HOST_B,
                                                ('observer_key','visitor_order'))['reference_excess_u'])
        visitor_drop.append({'taxon':name,'removed_records':int((group.partner_taxon==name).sum()),
            'new_fixed_strata_mean_excess':float(np.mean(parts)),
            'all_original_min_n_still_met':bool(eligible)})
    overall={
      'n_original_source_observation_ids':source['date_recovery']['dated_valid'],
      'n_same_cohort_records':len(group),'n_host_a':107,'n_host_b':84,
      'n_matched_region_months':len(cells),'n_original_observers':int(group.observer_key.nunique()),
      'species_identity_u':raw_u,
      'previous_observer_only_reference_mean_u':mean(lambda q:q['prior_observer_only']['reference_mean_u']),
      'previous_observer_only_excess_u':original_observer_excess,
      'observer_plus_order_reference_mean_u':strict_exp,
      'added_reference_from_insect_order_conditioning':order_additional_reference,
      'observer_plus_order_residual_u':strict_excess,
      'residual_fraction_of_previous_observer_excess':strict_excess/original_observer_excess,
      'n_observer_order_exchangeable_records':sum(q['observer_x_order']['n_exchangeable_records'] for q in strata),
      'n_observer_order_exchangeable_blocks':sum(q['observer_x_order']['n_exchangeable_blocks'] for q in strata),
      'post_hoc_observer_year_order_residual':mean(lambda q:q['observer_x_order_x_year']['reference_excess_u']),
      'post_hoc_observer_day_order_residual':mean(lambda q:q['observer_x_order_x_date']['reference_excess_u']),
      'n_post_hoc_observer_day_order_exchangeable_records':sum(q['observer_x_order_x_date']['n_exchangeable_records'] for q in strata),
      'observer_x_order_label_reference':{
        'n_replicates':reps,
        'mean_mc_u':float(np.mean(sim)),
        'exact_mean_u':strict_exp,
        'prob_upper_label_reference':float((1+(sim>=raw_u-1e-12).sum())/(reps+1)),
        'quantiles_025_50_975':np.quantile(sim,[.025,.5,.975]).tolist(),
        'randomized_field_experiment':False,
        'prospectively_preregistered':False},
      'hymenoptera_record_share':{host:source_order_counts[host]['Hymenoptera']/sum(source_order_counts[host].values()) for host in [HOST_A,HOST_B]},
      'hymenoptera_share_vulgare_minus_arvense_by_equal_cell':mean(
        lambda q:q['hymenoptera_versus_other_observer_fixed']['observed_host_b_minus_a']),
      'hymenoptera_share_observer_fixed_expected_by_equal_cell':mean(
        lambda q:q['hymenoptera_versus_other_observer_fixed']['observer_reference_difference']),
      'hymenoptera_share_observer_fixed_excess':mean(
        lambda q:q['hymenoptera_versus_other_observer_fixed']['observer_excess_difference']),
      'within_lepidoptera_name_identity':{'n_recorded_lepidoptera':sum(q['n_records'] for q in within_lep),
        'n_strata':len(within_lep),
        'mean_u':float(np.mean([q['observed_u'] for q in within_lep])),
        'mean_observer_reference':float(np.mean([q['label_reference_mean_u'] for q in within_lep])),
        'mean_excess':float(np.mean([q['reference_excess_u'] for q in within_lep]))},
      'leave_one_region_residual':leave_one_cell,
      'visitor_name_contribution_top_positive':sorted(name_contribution.items(),key=lambda x:-x[1])[:12],
      'visitor_name_contribution_top_negative':sorted(name_contribution.items(),key=lambda x:x[1])[:12],
      'species_residual_components_grouped_by_order':dict(sorted(type_contribution.items())),
      'leave_one_visitor':{
        'n_tested':len(visitor_drop),'n_still_meeting_original_minimum':sum(x['all_original_min_n_still_met'] for x in visitor_drop),
        'n_positive_fixed_strata_excess':sum(x['new_fixed_strata_mean_excess']>0 for x in visitor_drop),
        'all_fixed_strata_excess_range':[min(x['new_fixed_strata_mean_excess'] for x in visitor_drop),
            max(x['new_fixed_strata_mean_excess'] for x in visitor_drop)],
        'top_loss':sorted(visitor_drop,key=lambda x:x['new_fixed_strata_mean_excess'])[:7],
        'changes_source_support':True},
      'source_order_record_counts_by_host':source_order_counts,
      'source_order_record_share_equal_cell':{h:{o:float(np.mean(order_fractions[h][o])) for o in sorted(EXPECTED_COARSE_ORDERS)} for h in [HOST_A,HOST_B]},
    }
    result={'version':'cirsium_within_order_visitor_identity_v4',
       'source':{'archive_sha256':source['archive_sha256'],
          'metadata_sha256':source['metadata_sha256'],
          'prior_v3_repo_git_blob_sha':'b5ce6ee6c394f1769e2d612ca74d396d48aca040',
          'source_first_collected':'2026-10-08','analysis_date':'2026-10-20',
          'prior_support_is_same_5_regions_191_records':True},
       'status':'SOURCE_OBSERVER_ORDER_CONDITIONAL_LABEL_AUDIT_NOT_MORPHOLOGICAL_FITNESS',
       'source_support':{'n_hosts':2,'host_names':[HOST_A,HOST_B],
          'n_original_visited_insect_orders':len(EXPECTED_COARSE_ORDERS),
          'n_original_source_visitor_names':int(group.partner_taxon.nunique()),
          'n_original_source_observation_ids':len(group),
          'region_months':[[x,y] for x,y in cells]},
       'primary_observations_and_exact_reference':overall,
       'by_stratum':strata,
       'biological_claim_limit':[
         'These are recorded photographs, not measured visitation rates, pollination success, or available head resources.',
         'Hymenoptera and Lepidoptera are taxonomic orders, not verified functional pollination guilds.',
         'Same coarse 2-degree cell and month pool years, dates, local habitat and flower abundance.',
         'Observers recording only one host cannot exchange labels within observer/order and are not causal controls.',
         'The observer × order label null has been selected after seeing previous host association; P-tail is exploratory, conditional on fixed photo records, not a confirmatory trial.',
         'Leave-one-name retained frozen strata even if original >=10 threshold fails; support changes are separately flagged.',
         'Conditioning on order is a change in hypothetical label assignment, not identification of a causal mediated contribution.',
         'No capitulum head phenotype, physical contact, seed predator entry, viable achene fitness, ancestral molecular re-expression or permitted field manipulation measured.',
       ]}
    return summarize_audited_cohort(result)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--archive',required=True,type=Path)
    p.add_argument('--metadata',required=True,type=Path)
    p.add_argument('--prior-receipt',required=True,type=Path)
    p.add_argument('--out',required=True,type=Path)
    p.add_argument('--permutations',type=int,default=19999)
    args=p.parse_args()
    j=analyse_source(args.archive,args.metadata,args.prior_receipt,args.permutations)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(j,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
    s=j['primary_observations_and_exact_reference']
    print(json.dumps({'n':s['n_same_cohort_records'],
      'observer_excess':s['previous_observer_only_excess_u'],
      'observer_order_excess':s['observer_plus_order_residual_u'],
      'order_conditional_label_tail':s['observer_x_order_label_reference']['prob_upper_label_reference'],
      'hym_record_share_delta':s['hymenoptera_share_vulgare_minus_arvense_by_equal_cell']},indent=2))

if __name__=='__main__':main()