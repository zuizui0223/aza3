#!/usr/bin/env python3
"""Original-source Cirsium capitulum CORE vs outer phyllary-tip spread.

Real 2023 publisher XLSX, DOI 10.1007/s00606-023-01854-2 (Michalkova
et al., Online Resource 6). Frozen 114 original *nonhybrid same-plant*
records were exported with artifact_tool from the author source XLSX
verified against SHA256 c06423cc6f9800ba4775a4428afa072ba1b3dc6fb48ccd05fdb7e42a2f03ee79.

No external botanical DNA, fieldhead activity, pollination, true caliper
spine dimension, reproductive fitness or sampling authorization is inferred.
'Involucre spinosity or fuzzy' in the publisher means the ratio of involucre
width *including bract tips* to involucre width. It is *not* an anatomical
measurement of puncture force, spine length or actual guild accessibility.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

SOURCE_XLSX_SHA256="c06423cc6f9800ba4775a4428afa072ba1b3dc6fb48ccd05fdb7e42a2f03ee79"
SOURCE_CSV_SHA256="e69b20badd340e08b366ac51596661a2f8951b020994bc2f3154d838f1c2c9dc"
SOURCE_DOI="10.1007/s00606-023-01854-2"
SOURCE_URL="https://media.springernature.com/original/springer-static/esm/art%3A10.1007%2Fs00606-023-01854-2/MediaObjects/606_2023_1854_MOESM6_ESM.xlsx"
SOURCE_COLUMNS=("excel_row","taxon","sex","Cor_L","Inv_L","Inv_W","Ibr_W","Cap_W","Nec_P","Nec_W","Inv_Sp")
FIELDS=("Cor_L","Inv_L","Inv_W","Ibr_W","Cap_W","Nec_P","Nec_W")
PARENTAL={"Cirsium acaulon","Cirsium bertolonii - Monte Sagro","Cirsium erisithales"}
METRICS=("r_core","r_spread","delta","adjusted_r_core","adjusted_r_spread","adjusted_delta")
SEED=20261020


def require(ok, message):
    if not ok:
        raise ValueError(message)


def _as_row(r):
    """Minimal typed biological row; don't silently fill unknowns or zeros."""
    row=dict(r)
    require(row.get("taxon") not in (None,""),"MISSING_TAXON")
    require(row.get("sex") not in (None,""),"MISSING_SOURCE_SEX")
    for field in FIELDS+("Inv_Sp",):
        try: val=float(row[field])
        except (ValueError,TypeError,KeyError): raise ValueError("INVALID_BOTANICAL_NUMBER_"+field) from None
        require(math.isfinite(val) and val>0,"NONPOSITIVE_HEAD_MEASURE")
        row[field]=val
    require(math.isclose(row['Ibr_W']/row['Inv_W'],row['Inv_Sp'],rel_tol=0,abs_tol=1e-10),
            'BRAC_TIP_RATIO_MISMATCH')
    return row


def read_pinned_csv(path):
    path=Path(path)
    require(hashlib.sha256(path.read_bytes()).hexdigest()==SOURCE_CSV_SHA256,
            'FROZEN_SOURCE_CSV_SHA256_MISMATCH')
    with path.open(encoding='utf-8',newline='') as f:
        reader=csv.DictReader(f)
        require(tuple(reader.fieldnames or ())==SOURCE_COLUMNS,'SOURCE_COLUMNS_MISMATCH')
        rows=[_as_row(r) for r in reader]
    require(len(rows)==114,"REAL_SOURCE_114_PLANT_ROWS_DRIFT")
    require({r['taxon'] for r in rows}==PARENTAL and {r['sex'] for r in rows}=={'f','h'},
            'REAL_3_TAXON_2_SEX_PANEL_DRIFT')
    require(len({str(r['excel_row']) for r in rows})==114,'DUPLICATED_PUBLISHER_ROW')
    n=Counter((r['taxon'],r['sex']) for r in rows)
    require(sorted(n.values())==[10,10,10,24,30,30],"SOURCE_STRATUM_COUNTS_DRIFT")
    return rows


def residualized_log_columns(rows):
    parsed=[_as_row(r) for r in rows]
    require(len(parsed)>=6,"TOO_FEW_INDIVIDUALS")
    groups=defaultdict(list)
    for i,r in enumerate(parsed):groups[(r['taxon'],r['sex'])].append(i)
    require(all(len(v)>=3 for v in groups.values()),"TOO_FEW_PER_STRATUM")
    x={k:np.log([r[k] for r in parsed]).astype(float) for k in FIELDS}
    x['spread']=np.log([r['Ibr_W']/r['Inv_W'] for r in parsed]).astype(float)
    for a in x.values():
        for ids in groups.values():
            a[ids]-=a[ids].mean()
    return x,groups


def _cor(a,b,controls=()):
    aa,bb=a.copy(),b.copy()
    if controls:
        X=np.column_stack(controls)
        aa=aa-X@np.linalg.lstsq(X,aa,rcond=None)[0]
        bb=bb-X@np.linalg.lstsq(X,bb,rcond=None)[0]
    require(min(np.std(aa),np.std(bb))>1e-13,'DEGENERATE_VARIANCE_IN_CONTRAST')
    return float(np.corrcoef(aa,bb)[0,1])


def compute_metrics(rows):
    """Same-individual phenotypic covariance, not genetic modularity.

    All log variables demeaned separately within 3 species x 2 sex strata.
    Partial correlations are descriptive allometry sensitivity only; same
    width variable enters the tip-ratio denominator, so note denominator
    coupling and confirm adjusted direct-tip identity analytically.
    """
    x,_=residualized_log_columns(rows)
    core=_cor(x['Cor_L'],x['Inv_L'])
    outer=_cor(x['Cor_L'],x['spread'])
    controls=(x['Inv_W'],x['Cap_W'])
    adjcore=_cor(x['Cor_L'],x['Inv_L'],controls)
    adjouter=_cor(x['Cor_L'],x['spread'],controls)
    result={
        'n':len(rows),'n_strata':len({(r['taxon'],r['sex']) for r in rows}),
        'r_core':core,'r_spread':outer,'delta':core-outer,
        'adjusted_r_core':adjcore,'adjusted_r_spread':adjouter,
        'adjusted_delta':adjcore-adjouter,
        'r_tip_width_unadjusted':_cor(x['Cor_L'],x['Ibr_W']),
        'adjusted_r_tip_width':_cor(x['Cor_L'],x['Ibr_W'],controls),
    }
    require(abs(result['adjusted_r_tip_width']-adjouter)<1e-10,
            'RATIO_DENOMINATOR_PARTIAL_CORRELATION_IDENTITY_FAILED')
    return result


def bootstrap_metrics(rows, n_boot=5000, seed=SEED):
    require(isinstance(n_boot,int) and 100<=n_boot<=50000,'BOOTSTRAP_REPLICATE_BOUNDS')
    groups=defaultdict(list)
    for r in rows:groups[(r['taxon'],r['sex'])].append(r)
    require(len(groups)>=2,'MULTIPLE_STRATA_NEEDED_FOR_BOOTSTRAP')
    rng=np.random.default_rng(seed)
    draws={k:[] for k in METRICS}
    for _ in range(n_boot):
        smp=[]
        for group in sorted(groups):
            items=groups[group]
            smp.extend(items[int(j)] for j in rng.integers(len(items),size=len(items)))
        metric=compute_metrics(smp)
        for key in draws:draws[key].append(metric[key])
    result={
       'n_draws':n_boot,'seed':seed,
       'stratum_n':{f'{k[0]}|{k[1]}':len(v) for k,v in sorted(groups.items())},
       'ci95':{key:np.quantile(vals,[.025,.975]).tolist() for key,vals in draws.items()},
       'fraction_positive_difference':{
           'unadjusted':float(np.mean(np.array(draws['delta'])>0)),
           'size_adjusted':float(np.mean(np.array(draws['adjusted_delta'])>0))},
       'not_phylogenetic_or_independent_population_uncertainty':True,
    }
    return result


def analyze(rows,n_boot=5000,seed=SEED):
    require(len(rows)==114 and {r['taxon'] for r in rows}==PARENTAL,
            'REAL_ORIGINAL_PANEL_REQUIRED_FOR_PUBLIC_RECEIPT')
    obs=compute_metrics(rows)
    sensitivity={
        'leave_one_taxon_out':{taxon:compute_metrics([r for r in rows if r['taxon']!=taxon])
                               for taxon in sorted(PARENTAL)},
        'sex_specific':{sex:compute_metrics([r for r in rows if r['sex']==sex]) for sex in ('f','h')},
        'source_exact_spread_one_fraction':float(np.mean([r['Inv_Sp']==1.0 for r in rows])),
        'raw_ratio_per_taxon':{
            taxon:{'n':sum(r['taxon']==taxon for r in rows),
                   'median':float(np.median([r['Inv_Sp'] for r in rows if r['taxon']==taxon])),
                   'min':float(min(r['Inv_Sp'] for r in rows if r['taxon']==taxon)),
                   'max':float(max(r['Inv_Sp'] for r in rows if r['taxon']==taxon))}
            for taxon in sorted(PARENTAL)}
    }
    return {
        'version':'cirsium_head_core_external_source_real_morphometric_v1',
        'status':'DESCRIPTIVE_WITHIN_HEAD_PHENOTYPIC_COVARIANCE_ONLY',
        'as_of':'2026-10-20',
        'source':{'doi':SOURCE_DOI,'supplement':'Online Resource 6',
                  'title':'Hybridization may endanger the rare North Apennine endemic Cirsium bertolonii',
                  'publisher_xlsx_sha256':SOURCE_XLSX_SHA256,
                  'frozen_114head_derived_CSV_sha256':SOURCE_CSV_SHA256,
                  'n_real_nonhybrid_individuals':114,'sex_taxon_strata':6,
                  'source_csv_from_artifact_tool':True,
                  'row_selection':'all nonhybrid source rows; original Excel row indices retained',
                  'license_attribution':'Michálková et al. 2023, CC BY 4.0; derivative exports only published head dimensions'},
        'primary_observed':obs,'paired_plant_stratified_bootstrap':bootstrap_metrics(rows,n_boot,seed),
        'robustness':sensitivity,
        'estimand':'within-six-taxon-sex centered Pearson correlations of log(Cor_L) with log(Inv_L) vs log(Ibr_W/Inv_W); paired difference, with descriptive partial log Inv_W and log Cap_W controls',
        'claim_ceiling':{
            'one_physical_reproductive_head_is_not_one_synchronous_phenotype':True,
            'within_head_corr_not_developmental_pleiotropy':True,
            'true_anatomical_spine_length_measured':False,
            'functional_spine_defence_identified':False,
            'guild_specific_physical_contact_measured':False,
            'mature_viable_achenes_measured':False,
            'causal_reproductive_fitness_effect':False,
            'evolutionary_reexpression_or_reused_genetic_parts':False,
        },
        'caveats':[
            'Three species and fixed source localities; plant bootstrap is conditional on original design.',
            'Pooled within-stratum correlation can be dominated by hermaphrodite strata (84/114); female strata only 30.',
            'Ibr_W/Inv_W in author source is a dimensionless outer bract-tip width proxy, not true sharp spine length, contact force or physical accessibility.',
            'Cirsium acaulon ratio has near-floor values; leave-one-taxon-out and sex-specific estimates are descriptive robustness only.',
            'Shared size, growth stage, head rank, sex, botanical species ancestry and measurement methods can shape phenotypic correlation.',
            'Controlling ratio denominator Inv_W or head size can induce mathematical and causal model artifacts; adjustment is a sensitivity not proof.',
            'No observed pollinator, seed-feeder, parasitoid interaction, plant fitness, genomic architecture, or new field authorization.'
        ]
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--source-csv',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--bootstrap',type=int,default=5000)
    args=parser.parse_args()
    data=analyze(read_pinned_csv(args.source_csv),n_boot=args.bootstrap)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({'source_n':data['source']['n_real_nonhybrid_individuals'],
                      'primary':data['primary_observed'],
                      'ci95':data['paired_plant_stratified_bootstrap']['ci95'],
                      'status':data['status']},indent=2))

if __name__=='__main__':main()
