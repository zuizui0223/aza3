"""Head-only morphometric tests (synthetic correctness controls, no synthetic biological inference)."""
import sys
import unittest
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'analysis'))
from audit_cirsium_head_core_external_v1 import (  # noqa: E402
    read_pinned_csv, compute_metrics, bootstrap_metrics, analyze,
)


def fixture():
    rows=[]
    # The two stratum means are far apart: log centering must prevent
    # between-stratum differences from creating spurious covariance.
    for group,factor in [("A",1),("B",4)]:
        for i in range(8):
            r=float(i+1)
            width=factor*(8.0+r*.25)
            tip=width*(1.5 + 0.025*((3*i)%7))
            rows.append(dict(excel_row=len(rows)+20,taxon=group,sex="f",
                Cor_L=factor*(12+0.60*r+0.21*(i%3)),
                Inv_L=factor*(15+0.45*r+0.15*(i%2)),
                Inv_W=width,Ibr_W=tip,
                Cap_W=factor*(24+r*.5),Nec_P=factor*(6+r*.4),
                Nec_W=factor*(5+r*.2),Inv_Sp=tip/width))
    return rows

class HeadMorphometryTests(unittest.TestCase):
    def test_source_real_size_and_provenance(self):
        fn=Path(__file__).resolve().parents[1]/'data/evidence/cirsium_2023_same_plant_head_morphology_subset_v1.csv'
        rows=read_pinned_csv(fn)
        self.assertEqual(len(rows),114)
        self.assertEqual(len({x['excel_row'] for x in rows}),114)
        self.assertEqual(sorted(len([r for r in rows if (r['taxon'],r['sex'])==g])
                                for g in {(r['taxon'],r['sex']) for r in rows}),[10,10,10,24,30,30])
    def test_source_ratio_is_real_recorded_head_proxy(self):
        rows=fixture()
        rows[0]['Inv_Sp']+=0.08
        with self.assertRaisesRegex(ValueError,'BRAC_TIP_RATIO_MISMATCH'):
            compute_metrics(rows)
    def test_negative_and_zero_head_measurement_fail(self):
        for field in ['Cor_L','Inv_L','Inv_W','Ibr_W','Cap_W']:
            rows=fixture()
            rows[0][field]=0
            with self.assertRaisesRegex(ValueError,'NONPOSITIVE_HEAD_MEASURE'):
                compute_metrics(rows)
    def test_observational_core_correlation_is_bounded(self):
        m=compute_metrics(fixture())
        self.assertTrue(all(np.isfinite(m[k]) for k in ['r_core','r_spread','delta','adjusted_delta']))
        self.assertTrue(-1<=m['r_core']<=1)
    def test_group_mean_confounding_does_not_create_correlation(self):
        rows=fixture()
        m1=compute_metrics(rows)
        other=[dict(r,Cor_L=r['Cor_L']*10,Inv_L=r['Inv_L']*10,
                    Inv_W=r['Inv_W']*10,Ibr_W=r['Ibr_W']*10,Cap_W=r['Cap_W']*10)
               if r['taxon']=='B' else r.copy() for r in rows]
        m2=compute_metrics(other)
        for key in ['r_core','r_spread','delta','adjusted_delta']:
            self.assertAlmostEqual(m1[key],m2[key],places=12)
    def test_direct_tip_width_adjusted_equals_tip_spread_conditional_correlation(self):
        m=compute_metrics(fixture())
        self.assertAlmostEqual(m['adjusted_r_spread'],m['adjusted_r_tip_width'],places=11)
    def test_bootstrap_resamples_original_number_of_plants_per_stratum(self):
        out=bootstrap_metrics(fixture(),n_boot=100,seed=19)
        self.assertEqual(out['n_draws'],100)
        self.assertEqual(out['stratum_n'],{'A|f':8,'B|f':8})
        self.assertEqual(set(out['ci95']),{'r_core','r_spread','delta','adjusted_r_core','adjusted_r_spread','adjusted_delta'})
        self.assertTrue(all(hi>=lo for lo,hi in out['ci95'].values()))
    def test_invalid_duplicate_and_corrupted_source_fail_closed(self):
        fn=Path(__file__).resolve().parents[1]/'data/evidence/cirsium_2023_same_plant_head_morphology_subset_v1.csv'
        raw=fn.read_bytes()
        corrupted=fn.parent/'_negative_corrupted.csv'
        try:
            corrupted.write_bytes(raw.replace(b'Cirsium acaulon',b'Cirsium acaulon_typo',1))
            with self.assertRaisesRegex(ValueError,'FROZEN_SOURCE_CSV_SHA256_MISMATCH'):
                read_pinned_csv(corrupted)
        finally:
            corrupted.unlink(missing_ok=True)
    def test_real_cohort_metrics_are_source_based_not_simulated(self):
        fn=Path(__file__).resolve().parents[1]/'data/evidence/cirsium_2023_same_plant_head_morphology_subset_v1.csv'
        out=analyze(read_pinned_csv(fn),n_boot=100,seed=20261020)
        self.assertEqual(out['source']['n_real_nonhybrid_individuals'],114)
        self.assertAlmostEqual(out['primary_observed']['r_core'],.669896381550174,places=11)
        self.assertAlmostEqual(out['primary_observed']['r_spread'],-.05296790457029559,places=11)
        self.assertEqual(out['status'],'DESCRIPTIVE_WITHIN_HEAD_PHENOTYPIC_COVARIANCE_ONLY')
        self.assertFalse(out['claim_ceiling']['functional_spine_defence_identified'])

if __name__=='__main__':unittest.main()
