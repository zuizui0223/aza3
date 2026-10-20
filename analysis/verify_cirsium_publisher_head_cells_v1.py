#!/usr/bin/env python3
"""Check a published morphometry XLSX against the reduced aza3 head-only CSV.

Uses artifact_tool to inspect the exact publisher workbook and checks each
selected original Excel row, sex/taxon, and all eight head measurements.
This is a reproduction gate, not additional field observations.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

from artifact_tool import Blob, SpreadsheetFile

EXPECTED_XLSX = 'c06423cc6f9800ba4775a4428afa072ba1b3dc6fb48ccd05fdb7e42a2f03ee79'
EXPECTED_CSV = 'e69b20badd340e08b366ac51596661a2f8951b020994bc2f3154d838f1c2c9dc'
BOTANICAL_COLUMNS = ('Cor_L', 'Inv_L', 'Inv_W', 'Ibr_W', 'Cap_W', 'Nec_P', 'Nec_W', 'Inv_Sp')

def verify(xlsx_path: Path, csv_path: Path):
    xsha = hashlib.sha256(xlsx_path.read_bytes()).hexdigest()
    csha = hashlib.sha256(csv_path.read_bytes()).hexdigest()
    if xsha != EXPECTED_XLSX or csha != EXPECTED_CSV:
        raise ValueError('PINNED_SOURCE_BYTES_CHANGED')
    wb = SpreadsheetFile.import_xlsx(Blob.load(str(xlsx_path)))
    sheet = wb.worksheets.get_item_at(0)
    if sheet.name != 'S5 Morphometry':
        raise ValueError('WRONG_AUTHOR_WORKSHEET')
    grid = sheet.get_range('A9:AE181').values
    codes = {str(c):i for i,c in enumerate(grid[1]) if c is not None}
    if not set(BOTANICAL_COLUMNS).issubset(codes):
        raise ValueError('SOURCE_COLUMN_CODE_MISMATCH')
    with csv_path.open(encoding='utf-8', newline='') as f:
        rows = list(csv.DictReader(f))
    if len(rows) != 114 or len({int(r['excel_row']) for r in rows}) != 114:
        raise ValueError('SOURCE_COHORT_N_DRIFT')
    checks = 0
    for row in rows:
        idx = int(row['excel_row'])
        if not 12 <= idx <= 181:
            raise ValueError('SOURCE_EXCEL_ROW_OUT_OF_RANGE')
        original = grid[idx - 9]
        if str(original[1]) != row['taxon'] or str(original[3]) != row['sex']:
            raise ValueError(f'SOURCE_TAXON_OR_SEX_MISMATCH_ROW_{idx}')
        for field in BOTANICAL_COLUMNS:
            x = original[codes[field]]
            y = float(row[field])
            if not isinstance(x,(int,float)) or not math.isclose(float(x), y, rel_tol=1e-12, abs_tol=1e-10):
                raise ValueError(f'SOURCE_NUMERIC_MISMATCH_ROW_{idx}_{field}')
            checks += 1
    return {
        'version':'cirsium_michalkova_2023_original_publisher_cell_audit_v1',
        'status':'PASS_ALL_SELECTED_AUTHOR_CELLS',
        'xlsx_sha256':xsha, 'head_subset_csv_sha256':csha,
        'n_nonhybrid_source_plants':len(rows),
        'head_numeric_columns':list(BOTANICAL_COLUMNS),
        'source_numeric_cells_confirmed':checks,
        'source_sheet':sheet.name,
        'provenance_only_not_an_independent_experiment':True,
    }

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--xlsx',required=True,type=Path)
    p.add_argument('--subset-csv',required=True,type=Path)
    p.add_argument('--out',type=Path,required=True)
    args = p.parse_args()
    receipt = verify(args.xlsx,args.subset_csv)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps(receipt,ensure_ascii=False))

if __name__ == '__main__': main()
