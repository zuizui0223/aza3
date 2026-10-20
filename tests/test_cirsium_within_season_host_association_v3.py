"""Independent combinatorial checks; tiny data are tests, not field evidence."""
import importlib
import itertools
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'analysis'))
try:
    m = importlib.import_module('cirsium_within_season_host_association_v3')
except ModuleNotFoundError:
    m = None


def call(name, *a, **kw):
    assert m is not None, 'v3 implementation not written yet'
    return getattr(m, name)(*a, **kw)


def direct_u(a, b):
    within_a = sum(x == y for i,x in enumerate(a) for j,y in enumerate(a) if i != j)/(len(a)*(len(a)-1))
    within_b = sum(x == y for i,x in enumerate(b) for j,y in enumerate(b) if i != j)/(len(b)*(len(b)-1))
    cross = sum(x == y for x in a for y in b)/(len(a)*len(b))
    return (within_a+within_b)/2-cross


def allocations(a,b):
    pool = np.repeat(np.arange(len(a)), np.array(a)+np.array(b))
    for taken in itertools.combinations(range(len(pool)), int(sum(a))):
        mask = np.zeros(len(pool), dtype=bool);mask[list(taken)] = True
        yield np.bincount(pool[mask], minlength=len(a)),np.bincount(pool[~mask], minlength=len(a))


def example():
    rows=[]
    for c in ['c1','c2','c3']:
        for host,visitor in [('A','x'),('B','y')]:
            for i in range(4):
                rows.append(dict(event_id=f'{c}-{host}-{i}', focal_taxon=host,
                    partner_taxon=visitor,cell=c,month='07',year='2020',
                    year_month='2020-07',date=f'2020-07-{i+1:02d}',observer_key='o'+str(i)))
    return pd.DataFrame(rows)


def test_identical_single_category_zero():
    assert call('identity_contrast',[7,0],[3,0])==pytest.approx(0)


def test_disjoint_single_categories_one():
    assert call('identity_contrast',[6,0],[0,4])==pytest.approx(1)


def test_negative_estimate_not_truncated():
    assert call('identity_contrast',[1,1],[1,1])==pytest.approx(-.5)


def test_formula_matches_ordered_pair_enumeration():
    a=[0,0,0,1,2];b=[0,1,1,2]
    assert call('identity_contrast',[3,1,1],[1,2,1])==pytest.approx(direct_u(a,b))


def test_label_swap_and_zeros():
    a=[4,2,1];b=[1,1,3]
    assert call('identity_contrast',a,b)==pytest.approx(call('identity_contrast',b,a))
    assert call('identity_contrast',a+[0],b+[0])==pytest.approx(call('identity_contrast',a,b))


@pytest.mark.parametrize('a,b', [([1],[1]),([-1,3],[2,2]),([1.5,2],[2,2]),([2,2],[3]),([float('nan'),2],[2,2])])
def test_invalid_counts(a,b):
    assert m is not None, 'v3 implementation not written yet'
    with pytest.raises(ValueError):m.identity_contrast(a,b)


def test_unrestricted_exact_null_mean_zero():
    samples=list(allocations([2,1],[1,2]))
    brute=np.mean([call('identity_contrast',a,b) for a,b in samples])
    assert brute==pytest.approx(0,abs=1e-12)
    assert call('block_reference_expectation',[(np.array([2,1]),np.array([1,2]))])==pytest.approx(brute,abs=1e-12)


def test_observer_conditioned_mean_matches_all_block_allocations():
    blocks=[(np.array([2,0]),np.array([0,1])),(np.array([0,1]),np.array([1,1]))]
    values=[]
    for choices in itertools.product(*(list(allocations(a,b)) for a,b in blocks)):
        a=sum((c[0] for c in choices));b=sum((c[1] for c in choices))
        values.append(call('identity_contrast',a,b))
    assert call('block_reference_expectation',blocks)==pytest.approx(np.mean(values),abs=1e-12)


def test_single_host_observers_no_exchange_information():
    blocks=[(np.array([3,0]),np.array([0,0])),(np.array([0,0]),np.array([0,3]))]
    assert call('block_reference_expectation',blocks)==pytest.approx(1)
    draws=call('reference_draws',blocks,100,30)
    assert np.allclose(draws,1)


def test_reference_draws_reproducible_and_analytic_mean():
    blocks=[(np.array([4,1]),np.array([0,2])),(np.array([0,3]),np.array([3,0]))]
    d=call('reference_draws',blocks,12000,17)
    assert np.array_equal(d,call('reference_draws',blocks,12000,17))
    assert np.mean(d)==pytest.approx(call('block_reference_expectation',blocks),abs=.01)


def test_duplicate_ids_rejected():
    assert m is not None,'v3 implementation not written yet'
    d=example();d=pd.concat([d,d.iloc[:1]],ignore_index=True)
    with pytest.raises(ValueError,match='DUPLICATE'):m.analyse(d,'month',4,permutations=99)


def test_new_estimand_allows_different_matched_months_across_regions():
    d=example();d.loc[d.cell=='c3','month']='08'
    r=call('analyse',d,'month',4,permutations=99)
    assert r['n_replicated_host_pairs']==1
    q=r['host_pairs'][0]
    assert q['n_regions']==3 and q['n_strata']==3
    assert q['observed_u']==pytest.approx(1)


def test_insufficient_geographic_replication_not_zero():
    d=example();d=d[d.cell=='c1']
    r=call('analyse',d,'month',4,permutations=99)
    assert r['n_replicated_host_pairs']==0
    assert r['equal_replicated_pair_mean_u'] is None
    assert r['host_pairs'][0]['status']=='EXPLORATORY_FEWER_THAN_THREE_REGIONS'


def test_empty_support_not_zero():
    r=call('analyse',example(),'month',10,permutations=99)
    assert r['n_host_pairs']==0 and r['equal_replicated_pair_mean_u'] is None


def test_region_equal_weight_not_record_or_stratum_count():
    got=call('region_balanced_mean', [('a',1.),('a',1.),('b',0.)])
    assert got==pytest.approx(.5)


def test_missing_observers_are_not_a_shared_person():
    d=example();d['observer_key']=''
    r=call('analyse',d,'month',4,permutations=99)
    q=r['host_pairs'][0]
    assert q['observer_reference_mean_u']==pytest.approx(q['observed_u'])
    assert q['n_exchangeable_observer_block_records']==0


def test_observer_year_blocks_do_not_exchange_across_years():
    d=example();d.loc[d.focal_taxon=='B','year']='2021'
    q=call('analyse',d,'month',4,permutations=99)['host_pairs'][0]
    r=call('additional_reference_for_pair',d,q,['observer_key','year'],99)
    assert r['n_exchangeable_records']==0
    assert r['exact_reference_mean_u']==pytest.approx(q['observed_u'])
    assert r['observed_minus_reference_u']==pytest.approx(0)


def test_same_observer_same_day_blocks_cannot_exchange_different_dates():
    d=example();d.loc[d.focal_taxon=='B','date']='2020-08-01'
    q=call('analyse',d,'month',4,permutations=99)['host_pairs'][0]
    r=call('additional_reference_for_pair',d,q,['observer_key','date'],99)
    assert r['n_exchangeable_records']==0
    assert r['reference']['upper_tail_label_reference_probability']==pytest.approx(1)
