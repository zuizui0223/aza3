"""Negative and positive unit checks for source-authenticated Cirsium four-arm audit."""
from __future__ import annotations
import unittest
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"analysis"))
from audit_cirsium_2022_four_arm_head_outcomes_v1 import (
    parse_rows,group_estimate,prep,sign_counts,summarize_contrast,
    synthetic_controls,complete_block_summary,check,
)

def row(site="campus",block="B1",arm="Control",flowers=10,viable=7,
        nonviable=2,damaged=1,missing=0,insects=0):
    return {"Site_Code":site,"Block_Code":block,"Treatment_Code":arm,
            "Vial_Code":"1","Flowers":str(flowers),
            "Viable_Seeds":str(viable),"NViable_Seeds":str(nonviable),
            "Damaged_Seeds":str(damaged),"Missing_Seeds":str(missing),
            "Insects":str(insects)}

class FieldAudit(unittest.TestCase):
    def test_synthetic_gate(self):
        self.assertTrue(synthetic_controls().startswith("PASS_"))

    def test_bad_git_blob_rejected_before_biology(self):
        with self.assertRaisesRegex(ValueError,"GIT_BLOB_CHANGED"):
            parse_rows(b"Site_Code,Block_Code,Treatment_Code\n",strict=True)

    def test_missing_required_column(self):
        with self.assertRaisesRegex(ValueError,"ORIGINAL_SCHEMA_MISSING"):
            parse_rows(b"Site_Code,Block_Code\ncampus,B1\n",strict=False)

    def test_negative_accounting_residual_kept(self):
        r=row(viable=8,nonviable=1,damaged=2,missing=-1)
        b=( ",".join(r.keys())+"\n"+",".join(r.values())+"\n").encode()
        x,h=parse_rows(b,strict=False)
        s=group_estimate(x)
        self.assertAlmostEqual(s["signed_missing_fraction"],-.1)
        self.assertAlmostEqual(s["inferred_eaten_signed_fraction"],.1)

    def test_bad_mass_balance_rejected(self):
        z=row(viable=8,nonviable=1,damaged=2,missing=0)
        b=( ",".join(z.keys())+"\n"+",".join(z.values())+"\n").encode()
        with self.assertRaisesRegex(ValueError,"SEED_ACCOUNTING_IDENTITY_BROKEN"):
            parse_rows(b,strict=False)

    def test_zero_flowers_rejected(self):
        z=row(flowers=0,viable=0,nonviable=0,damaged=0)
        b=( ",".join(z.keys())+"\n"+",".join(z.values())+"\n").encode()
        with self.assertRaisesRegex(ValueError,"FLOWERS_DENOMINATOR_ZERO"):
            parse_rows(b,strict=False)

    def test_negative_viable_rejected(self):
        z=row(viable=-1,nonviable=9,damaged=2,missing=0)
        b=( ",".join(z.keys())+"\n"+",".join(z.values())+"\n").encode()
        with self.assertRaisesRegex(ValueError,"INVALID_HEAD_COUNT"):
            parse_rows(b,strict=False)

    def test_group_uses_sum_flowers_not_mean_head_fraction(self):
        a=row(flowers=10,viable=10,nonviable=0,damaged=0,missing=0)
        b=row(flowers=100,viable=0,nonviable=100,damaged=0,missing=0)
        s=group_estimate([a,b])
        self.assertAlmostEqual(s["viable_fraction"],10/110)
        self.assertNotAlmostEqual(s["viable_fraction"],.5)
        self.assertEqual(s["n_sampled_heads"],2)

    def test_comparison_is_within_block_and_arm(self):
        z=[
            row("campus","A","Control",viable=9,nonviable=1,damaged=0),
            row("campus","A","Pollinator Exclusion",viable=1,nonviable=9,damaged=0),
            row("campus","A","Predator Exclusion",viable=8,nonviable=2,damaged=0),
            row("campus","A","Predator Model",viable=8,nonviable=2,damaged=0),
            row("ditch","A","Control",viable=7,nonviable=3,damaged=0),
            row("ditch","B","Pollinator Exclusion",viable=0,nonviable=10,damaged=0),
        ]
        _,proxy,blocks,complete=prep(z,False)
        self.assertEqual(len(complete),1)
        o=summarize_contrast(blocks,"Pollinator Exclusion","Control")
        self.assertEqual(o["n_matched_blocks"],1)
        self.assertAlmostEqual(o["outcomes"]["viable_fraction"]["mean_a_minus_b"],-.8)

    def test_complete_block_requires_all_four_sites(self):
        rows=[row("campus","A",tr,viable=7,nonviable=3,damaged=0)
              for tr in ["Control","Pollinator Exclusion","Predator Exclusion","Predator Model"]]
        complete=prep(rows,False)[3]
        with self.assertRaisesRegex(ValueError,"REQUIRE_ALL_FOUR_ORIGINAL_STUDY_SITES"):
            complete_block_summary(complete,bootstrap=100)

    def test_signed_missing_inferred_identity(self):
        z=[row(viable=7,nonviable=1,damaged=0,missing=2),
           row(viable=3,nonviable=1,damaged=2,missing=4)]
        s=group_estimate(z)
        self.assertAlmostEqual(s["inferred_eaten_signed_fraction"],
            s["observed_damaged_fraction"]+s["signed_missing_fraction"])
        self.assertEqual(sign_counts([1,-1,0]),{"positive":1,"negative":1,"zero":1})

if __name__=="__main__":
    unittest.main()
