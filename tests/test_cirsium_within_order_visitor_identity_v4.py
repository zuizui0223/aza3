"""TDD: preserve the original 191-record source cohort; test observer×order null."""
from __future__ import annotations
import itertools
import math
import sys
from pathlib import Path
import unittest
import numpy as np
import pandas as pd

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'analysis'))
from cirsium_within_order_visitor_identity_v4 import (
    source_insect_order, exact_category_moments, category_identity_components,
    observer_fixed_order_difference, order_conditioned_units,
    frozen_cohort, summarize_audited_cohort)
from cirsium_within_season_host_association_v3 import identity_contrast,block_reference_expectation

class TestWithinOrderVisitorIdentityV4(unittest.TestCase):
    def test_unique_order_only(self):
        self.assertEqual(source_insect_order('Eukaryota | Insecta | Lepidoptera | Danaus plexippus'),'Lepidoptera')
        self.assertEqual(source_insect_order('Insecta | Hymenoptera | Bombus impatiens'),'Hymenoptera')
        with self.assertRaisesRegex(ValueError,'ORDER_AMBIGUOUS'):
            source_insect_order('Insecta | Lepidoptera | Hymenoptera')
        with self.assertRaisesRegex(ValueError,'ORDER_AMBIGUOUS'):
            source_insect_order('Insecta')

    def test_exact_category_moments_from_literal_observation_allocations(self):
        # Two observer/order blocks, probabilities follow hypergeometric counts.
        block_values=[(['A','A','B','C'],2),(['A','B','C'],1)]
        blocks=[(np.array([2,0,0]),np.array([0,1,1])),
                (np.array([0,1,0]),np.array([1,0,1]))]
        mu,v=exact_category_moments(blocks)
        self.assertAlmostEqual(mu[0],2*2/4+1/3)
        self.assertTrue(np.all(v>=0))
        ua,ub=sum((b[0] for b in blocks)),sum((b[1] for b in blocks))
        observed,ref,diff=category_identity_components(blocks)
        self.assertTrue(math.isclose(observed.sum(),identity_contrast(ua,ub),abs_tol=1e-12))
        self.assertTrue(math.isclose(ref.sum(),block_reference_expectation(blocks),abs_tol=1e-12))
        self.assertTrue(np.allclose(diff,observed-ref))
        # Exhaustive assignment of each observer's host labels to their fixed visitors.
        outcomes=[]
        for ia in itertools.combinations(range(4),2):
            for ib in itertools.combinations(range(3),1):
                x=np.zeros(3,dtype=int);y=np.zeros(3,dtype=int)
                for letters,keep in [(block_values[0][0],set(ia)),(block_values[1][0],set(ib))]:
                    for i,t in enumerate(letters):
                        idx='ABC'.index(t)
                        if i in keep:x[idx]+=1
                        else:y[idx]+=1
                outcomes.append(identity_contrast(x,y))
        self.assertEqual(len(outcomes),18)
        self.assertAlmostEqual(float(np.mean(outcomes)),float(ref.sum()),places=12)

    def test_no_false_positive_when_order_observer_labels_cannot_exchange(self):
        g=pd.DataFrame([
            ('A','Lepidoptera','s1','o1'),('A','Lepidoptera','s1','o1'),
            ('B','Hymenoptera','s2','o1'),('B','Hymenoptera','s2','o1'),
            ('A','Hymenoptera','s3','o2'),('A','Hymenoptera','s3','o2'),
            ('B','Lepidoptera','s4','o2'),('B','Lepidoptera','s4','o2')],
            columns=['focal_taxon','visitor_order','partner_taxon','observer_key'])
        s=order_conditioned_units(g,'A','B',('observer_key','visitor_order'))
        self.assertEqual(s['n_exchangeable_records'],0)
        self.assertEqual(s['n_exchangeable_blocks'],0)
        self.assertAlmostEqual(s['observer_order_excess_u'],0)

    def test_one_observer_two_hosts_within_order_is_exchangeable(self):
        g=pd.DataFrame([
            ('A','Lepidoptera','s1','o1'),('B','Lepidoptera','s2','o1'),
            ('A','Hymenoptera','s3','o1'),('B','Hymenoptera','s3','o2'),
            ('A','Lepidoptera','s1','o3'),('B','Lepidoptera','s2','o4')],
            columns=['focal_taxon','visitor_order','partner_taxon','observer_key'])
        s=order_conditioned_units(g,'A','B',('observer_key','visitor_order'))
        self.assertEqual(s['n_exchangeable_records'],2)
        self.assertEqual(s['n_exchangeable_blocks'],1)
        self.assertTrue(math.isfinite(s['observer_order_excess_u']))

    def test_observer_fixed_hym_difference_exact_enumeration(self):
        g=pd.DataFrame([
            ('A','Lepidoptera','o1'),('A','Hymenoptera','o1'),
            ('B','Hymenoptera','o1'),('B','Lepidoptera','o1'),
            ('A','Hymenoptera','o2'),('B','Lepidoptera','o2')],
            columns=['focal_taxon','visitor_order','observer_key'])
        s=observer_fixed_order_difference(g,'A','B','Hymenoptera')
        self.assertAlmostEqual(s['observed_host_b_minus_a'],-1/3)
        self.assertAlmostEqual(s['observer_reference_difference'],0.)
        self.assertEqual(s['n_exchangeable_records'],6)

    def test_frozen_cohort_must_match_predeclared_source_counts(self):
        # The earlier v3 cohort has been frozen; never reselect based on new outcome.
        rows=[]
        for cell,month,na,nb in [
          ('65_50','07',39,15),('65_52','08',12,33),('66_50','07',29,12),
          ('66_54','07',15,14),('66_55','08',12,10)]:
          rows.extend([{'event_id':f'{cell}_{month}_A_{i}','cell':cell,'month':month,
                        'focal_taxon':'Cirsium arvense'} for i in range(na)])
          rows.extend([{'event_id':f'{cell}_{month}_B_{i}','cell':cell,'month':month,
                        'focal_taxon':'Cirsium vulgare'} for i in range(nb)])
        source=pd.DataFrame(rows)
        expected={('65_50','07'):(39,15),('65_52','08'):(12,33),('66_50','07'):(29,12),
                  ('66_54','07'):(15,14),('66_55','08'):(12,10)}
        sub=frozen_cohort(source,expected)
        self.assertEqual(len(sub),191)
        with self.assertRaisesRegex(ValueError,'SOURCE_PRIMARY_COUNT_DRIFT'):
            frozen_cohort(source.iloc[1:],expected)
        with self.assertRaisesRegex(ValueError,'SOURCE_PRIMARY_DUPLICATE_ID'):
            frozen_cohort(pd.concat([source,source.iloc[:1]]),expected)

    def test_result_does_not_store_event_location_or_observer(self):
        # Public report output is aggregate, source IDs/coordinates stay in transient source.
        x={'result':'a', 'nested':{'count':3}}
        with self.assertRaisesRegex(ValueError,'REDACTED_FIELD_PRESENT'):
            summarize_audited_cohort({**x,'observer_key':'secret'})
        self.assertEqual(summarize_audited_cohort(x),x)

if __name__=='__main__':unittest.main()