#!/usr/bin/env python3
from __future__ import annotations
import csv, json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
CONTRACT=ROOT/'data/planning/chapter3_backbone_hole_fill_v6.json'
REPAIR=ROOT/'data/planning/chapter3_moreyra_anchor_repair_v6.csv'
DOC=ROOT/'docs/CHAPTER3_BACKBONE_HOLE_FILL_V6.md'

def main()->int:
    d=json.loads(CONTRACT.read_text(encoding='utf-8'))
    assert d['contract_version']=='chapter3_backbone_hole_fill_v6'
    assert d['authoritative_for_level1_acquisition_priority'] is True
    assert d['planning_species_count']==120
    assert d['public_skeleton_species']==38
    h={x['hole_id']:x for x in d['hole_classes']}
    assert h['H1_ABSENT_FROM_PUBLIC_SKELETON']['nominal_species']==82
    assert h['H1_ABSENT_FROM_PUBLIC_SKELETON']['nominal_new_individuals']==164
    assert h['H2_WEAK_OR_BROKEN_MOREYRA_TIP']['nominal_species']==9
    assert h['H2_WEAK_OR_BROKEN_MOREYRA_TIP']['nominal_new_individuals']==18
    assert h['H3_SINGLE_GOOD_PUBLIC_TIP_NEEDS_OWN_TRAIT_LINK']['nominal_species']==29
    assert h['H3_SINGLE_GOOD_PUBLIC_TIP_NEEDS_OWN_TRAIT_LINK']['nominal_new_individuals']==29
    counts=d['nominal_counts_at_120_species']
    assert counts['phase_A_empty_and_broken_holes_new_individuals']==182
    assert counts['phase_B_good_public_tip_trait_link_new_individuals']==29
    assert counts['new_individuals_to_reach_two_representatives_per_species_using_public_skeleton']==211
    assert counts['primary_two_representative_nationwide_matrix']==240
    assert d['formula']['at_N120_M38_R9']==211
    assert d['weak_tip_registry']['count']==9
    rows=list(csv.DictReader(REPAIR.open(encoding='utf-8',newline='')))
    assert len(rows)==9
    assert sum(int(r['new_primary_samples']) for r in rows)==18
    assert {r['member_id'] for r in rows}=={'JPN_29','JPN_31','JPN_33','JPN_32','JPN_34','JPN_35','JPN_36','JPN_37','JPN_38'}
    doc=DOC.read_text(encoding='utf-8')
    for term in ('164 new individuals','18 new individuals','29 new individuals','211 new individuals','240 primary representatives','Exact branch placement'):
        assert term in doc, term
    state=d['current_state']
    assert state['physical_samples']==0
    assert all(v is False for k,v in state.items() if k!='physical_samples')
    print('chapter3_backbone_hole_fill_v6_valid=true')
    print('phase_A_new=182')
    print('phase_B_new=29')
    print('two_representative_new_total=211')
    return 0
if __name__=='__main__':
    raise SystemExit(main())
