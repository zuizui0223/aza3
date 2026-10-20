import importlib.util
from pathlib import Path
import sys
import unittest
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'analysis'))
MODULE='cirsium_visit_dates_v2'

class DatesTest(unittest.TestCase):
    def module(self):
        self.assertIsNotNone(importlib.util.find_spec(MODULE),'Date recovery and same-season analysis module must exist')
        return __import__(MODULE)
    def test_metadata_is_minimal_and_observer_pseudonymous(self):
        m=self.module(); r=m.compact_observation({'id':1,'observed_on':'2020-07-14','user':{'id':9,'login':'do-not-store'},'geojson':{'coordinates':[0,0]},'taxon':{'id':44,'name':'Bombus test','ancestors':[{'name':'Bombus'}]}})
        self.assertEqual(r['event_id'],'1');self.assertEqual(r['observed_on'],'2020-07-14')
        self.assertEqual(len(r['observer_key']),64);self.assertNotIn('user',r);self.assertNotIn('geojson',r)
    def test_created_date_is_never_observation_date(self):
        m=self.module(); r=m.compact_observation({'id':1,'created_at':'2020-07-14'})
        self.assertEqual(r['observed_on'],'')
    def test_join_dates_rejects_future_partial_and_invalid(self):
        m=self.module(); e=pd.DataFrame({'event_id':['1','2','3','4','5']})
        meta=[{'event_id':str(i+1),'observed_on':d,'observer_key':''} for i,d in enumerate(['2020-07-14','2026-10-09','2020-07','2020-02-30',''])]
        x,a=m.join_dates(e,meta)
        self.assertEqual(len(x),1);self.assertEqual(x.iloc[0]['month'],'07');self.assertEqual(a['dated_valid'],1)
        self.assertEqual(a['after_source_snapshot'],1);self.assertEqual(a['invalid_or_partial_date'],2)
    def test_duplicate_or_unrequested_metadata_fails(self):
        m=self.module();e=pd.DataFrame({'event_id':['1']})
        for records in [[{'event_id':'1'},{'event_id':'1'}],[{'event_id':'2'}]]:
            with self.assertRaises(ValueError):m.join_dates(e,records)
    def test_missing_observers_are_not_collapsed(self):
        m=self.module();e=pd.DataFrame({'event_id':['1','2','3','4'], 'focal_taxon':['H']*4,'partner_taxon':['V']*4,'cell':['c']*4,'date':['2020-07-01']*4,'month':['07']*4,'observer_key':['x','x','','']})
        out=m.reduce_observer_days(e);self.assertEqual(set(out.event_id),{'1','3','4'})
    def fixture(self,split_month=False):
        rows=[];i=0
        for h in ['Cirsium a','Cirsium b']:
            for cell in ['c1','c2','c3']:
                for k in range(3):
                    i+=1;mo='07' if (not split_month or h=='Cirsium a') else '08'
                    rows.append(dict(event_id=str(i),focal_taxon=h,partner_taxon='Visitor '+h[-1],cell=cell,month=mo,year_month='2020-'+mo,date='2020-'+mo+'-01',observer_key=str(i)))
        return pd.DataFrame(rows)
    def test_season_never_mixes_months(self):
        m=self.module();x=self.fixture(True)
        out=m.season_summary(x,'month',2)
        self.assertEqual(out['n_host_pairs'],0);self.assertIsNone(out['equal_host_pair_mean_delta'])
    def test_three_region_same_month_species_signal(self):
        m=self.module();out=m.season_summary(self.fixture(),'month',2)
        self.assertEqual(out['n_host_pairs'],1);self.assertEqual(out['n_reused_rectangles'],3)
        self.assertAlmostEqual(out['equal_host_pair_mean_delta'],-1.)
    def test_year_month_is_not_month_pool(self):
        m=self.module();x=self.fixture();x.loc[x.focal_taxon=='Cirsium b','year_month']='2021-07'
        self.assertEqual(m.season_summary(x,'year_month',2)['n_host_pairs'],0)
        self.assertEqual(m.season_summary(x,'month',2)['n_host_pairs'],1)
    def test_coordinate_boundaries_match_previous_grid(self):
        m=self.module();x=pd.DataFrame({'lat':[90,0,-90],'lon':[180,0,-180]})
        self.assertEqual(m.add_cells(x)['cell'].tolist(),['89_0','45_90','0_0'])

if __name__=='__main__':unittest.main()
