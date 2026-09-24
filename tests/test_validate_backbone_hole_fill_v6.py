#!/usr/bin/env python3
import subprocess, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

def main()->int:
    proc=subprocess.run([sys.executable,str(ROOT/'analysis/validate_backbone_hole_fill_v6.py')],cwd=ROOT,text=True,capture_output=True)
    if proc.returncode!=0:
        print(proc.stdout)
        print(proc.stderr,file=sys.stderr)
        return proc.returncode
    assert 'chapter3_backbone_hole_fill_v6_valid=true' in proc.stdout
    assert 'two_representative_new_total=211' in proc.stdout
    print('backbone_hole_fill_v6_test=true')
    return 0
if __name__=='__main__':
    raise SystemExit(main())
