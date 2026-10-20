import itertools
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'analysis'))
from cirsium_flower_partner_regional_turnover_v1 import (
    rarefied_bray_curtis, rectangle_contrasts, prepare_events, complete_core,
)

def test_rarefaction_matches_exhaustive_two_record_draws():
    counts=np.array([[2,1,1],[1,2,1]])
    a=[0,0,1,2];b=[0,1,1,2]
    ds=[]
    for ix in itertools.combinations(range(4),2):
        for iy in itertools.combinations(range(4),2):
            ca=np.bincount([a[k] for k in ix],minlength=3)
            cb=np.bincount([b[k] for k in iy],minlength=3)
            ds.append(np.abs(ca-cb).sum()/4)
    got=rarefied_bray_curtis(counts,2)
    assert got is not None, 'exact rarefaction estimator missing'
    assert np.isclose(got[0,1],np.mean(ds))

def test_disjoint_partners_have_distance_one():
    got=rarefied_bray_curtis(np.array([[10,0],[0,10]]),5)
    assert got is not None
    assert got[0,1]==1

def test_no_false_population_self_distance_zero():
    got=rarefied_bray_curtis(np.array([[2,2],[2,2]]),2)
    assert got is not None
    assert got[0,1]>0 # independent finite draws from identical compositions

def test_reject_insufficient_record_count():
    with pytest.raises(ValueError):rarefied_bray_curtis(np.array([[1,1],[3,2]]),3)

def test_reject_fractional_counts():
    with pytest.raises(ValueError):rarefied_bray_curtis(np.array([[2.5,1.5],[3,1]]),2)

def test_region_only_rectangles_are_positive():
    units=[('A','c1'),('A','c2'),('A','c3'),('B','c1'),('B','c2'),('B','c3')]
    # Every host in one cell has the same exclusive visitor identity.
    d=np.array([[float(x[1]!=y[1]) for y in units] for x in units])
    got=rectangle_contrasts(d,units)
    assert got is not None
    assert len(got)==3
    assert np.allclose(got['regional_minus_host'],1)

def test_host_only_rectangles_are_negative():
    units=[('A','c1'),('A','c2'),('A','c3'),('B','c1'),('B','c2'),('B','c3')]
    d=np.array([[float(x[0]!=y[0]) for y in units] for x in units])
    got=rectangle_contrasts(d,units)
    assert got is not None
    assert np.allclose(got['regional_minus_host'],-1)

def test_two_shared_cells_do_not_pass_three_cell_gate():
    got=rectangle_contrasts(np.zeros((4,4)),[('A','1'),('A','2'),('B','1'),('B','2')])
    assert got is not None
    assert len(got)==0

def record(obs='1',host='Cirsium arvense',partner='Vanessa cardui',stage='adult'):
    return dict(focal_taxon=host,partner_taxon=partner,focal_position='target',
      interaction_type='visitsFlowersOf',
      relation_direction_status='DIRECTION_COMPATIBLE_NOT_FUNCTION_VERIFIED',
      study_title='https://www.inaturalist.org/observations/globi-observations-resource-relationships-dwca.zip',
      source_taxon_path='Animalia | Arthropoda | Insecta | Lepidoptera',
      source_specimen_occurrence_id='https://www.inaturalist.org/observations/'+obs,
      source_specimen_life_stage=stage,event_date='',latitude='45',longitude='10')

def test_duplicate_observation_cannot_inflate_sampling():
    out=prepare_events(pd.DataFrame([record(),record()]))
    assert out is not None
    events,audit=out
    assert len(events)==1

def test_ambiguous_host_event_is_excluded_not_arbitrarily_resolved():
    out=prepare_events(pd.DataFrame([record(),record(host='Cirsium vulgare')]))
    assert out is not None
    events,audit=out
    assert len(events)==0
    assert audit['ambiguous_observation_ids_excluded']==1

def test_reject_general_eating_and_nonadult_flower_rows():
    a=record();a['interaction_type']='eats'
    out=prepare_events(pd.DataFrame([a,record('2',stage='larva'),record('3')]))
    assert out is not None
    events,audit=out
    assert len(events)==1

def test_fish_and_unidentified_host_cannot_enter_insect_matrix():
    fish=record('1',partner='Platygobio gracilis');fish['source_taxon_path']='Animalia | Actinopterygii'
    out=prepare_events(pd.DataFrame([fish,record('2',host='Cirsium'),record('3')]))
    assert out is not None
    assert len(out[0])==1


def test_complete_core_uses_all_hosts_and_only_shared_regions():
    units=[(h,c) for h in ['A','B','C'] for c in ['1','2','3']]+[('A','4')]
    counts=pd.DataFrame([[10,0]]*len(units),index=pd.MultiIndex.from_tuples(units))
    rect=rectangle_contrasts(np.zeros((len(units),len(units))),units)
    got=complete_core(counts,rect,5)
    assert got is not None
    assert got['hosts']==['A','B','C']
    assert got['cells']==['1','2','3']
    assert got['n_recorded_visits']==90
    assert got['contrasts']['n_reused_rectangles']==9

def test_complete_core_does_not_select_a_favorable_host_subset():
    units=[('A','1'),('A','2'),('A','3'),('B','1'),('B','2'),('B','3'),
           ('B','4'),('B','5'),('B','6'),('C','4'),('C','5'),('C','6')]
    counts=pd.DataFrame([[10,0]]*len(units),index=pd.MultiIndex.from_tuples(units))
    rect=rectangle_contrasts(np.zeros((len(units),len(units))),units)
    got=complete_core(counts,rect,5)
    assert got is not None
    assert got['status']=='INSUFFICIENT_ALL_HOST_COMMON_REGIONS'
