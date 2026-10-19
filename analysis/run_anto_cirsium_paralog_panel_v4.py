"""Run multi-protein control panel on real, source-verified Cirsium transcripts.

No pollen experiments, Ryukyu 2026 source alignment, functional DFR/ANS
enzymology or floral re-expression data are implied.
"""
import argparse
import json
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"analysis/anto_empirical"))
from cirsium_paralog_specificity_v4 import run_panel

DATA=ROOT/"data/anto_source_handoff/molecular"


def main(output=None):
    report=run_panel(
        DATA/"cirsium_japonicum_2020_26_source_sequences_v6.fasta",
        DATA/"cirsium_japonicum_2020_26_sequence_provenance_v6.json",
        DATA/"extended_DFR_ANS_ANR_FLS_reference_proteins_v4.fasta",
        DATA/"extended_reference_protein_provenance_v4.json",
    )
    report["repository_owner"]="zuizui0223/aza3"
    report["upstream_theoretical_programme"]="zuizui0223/anto"
    text=json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True)+"\n"
    if output is None:print(text)
    else:
        output=Path(output)
        output.parent.mkdir(parents=True,exist_ok=True)
        output.write_text(text,encoding="utf-8")
    return report


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path)
    args=ap.parse_args()
    main(args.out)
