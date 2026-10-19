# Source migrated from zuizui0223/anto/src/anto_cirsium_organ_evidence_v6.py (anto v4-v6); empirical ownership now aza3.
"""Author-reported *Cirsium japonicum* organ contrasts, never Ryukyu observations.

The 2020 study (DOI 10.1021/acsomega.0c04001) measured expression
and phenylpropanoid metabolites in flowers, leaves, stems and roots.
A reported flower:leaf CHS ratio is NOT an anthocyanin-specific assay.
A reported root:leaf DFR ratio is NOT a root:flower comparison.
Do not pool fold ratios of different genes, experimental outcomes or
different denominators as a single effect size.
"""
from collections import Counter
import json
from math import log2, isfinite
from pathlib import Path

DOI = "10.1021/acsomega.0c04001"
EXPRESSION = {
    "Park2020_CjCHS1": ("CHS1", "flower", "leaf", 500.0),
    "Park2020_CjCHS2": ("CHS2", "flower", "leaf", 112.0),
    "Park2020_CjF3pH": ("F3pH", "flower", "leaf", 16.6),
    "Park2020_CjCHI": ("CHI", "flower", "leaf", 13.0),
    "Park2020_CjFLS": ("FLS", "flower", "leaf", 5.6),
    "Park2020_CjFSII": ("FSII", "flower", "leaf", 3.3),
    "Park2020_CjDFR_root_isoform_A": ("DFR_isoform_unresolved_A", "root", "leaf", 22.0),
    "Park2020_CjDFR_root_isoform_B": ("DFR_isoform_unresolved_B", "root", "leaf", 2.6),
}
METABOLITES = {
    "Park2020_total_flavonols": ("total_flavonols", "flower", "root", 21.0),
    "Park2020_total_flavones": ("total_flavones", "leaf", "flower", 58.0),
    "Park2020_total_flavanols": ("total_flavanols", "leaf", "root", 19.0),
}


def validate_registry(path):
    doc = json.loads(Path(path).read_text(encoding="utf-8"))
    if doc.get("version") != "anto_cirsium_japonicum_organ_evidence_v6":
        raise ValueError("Unknown evidence registry")
    if doc.get("floral_anthocyanin_directly_measured") is not False:
        raise ValueError("This author study has no admitted direct anthocyanin chemistry")
    if doc.get("six_ryukyu_specimens_linked") is not False:
        raise ValueError("Cannot splice this study onto Ryukyu specimens")
    if doc.get("trait_regain_event_observed") is not False:
        raise ValueError("No trait regain event was measured")
    if doc.get("within_system_independent_origin_replicates") != 0:
        raise ValueError("Tissue samples are not independent evolutionary reversals")
    if doc.get("published_transcriptome_unigenes") != 104890:
        raise ValueError("Published source transcriptome denominator altered")
    recs = doc.get("records")
    if not isinstance(recs, list) or len(recs) != 11:
        raise ValueError("Expect exactly eight mRNA and three metabolite contrasts")
    seen=set()
    for r in recs:
        name=r.get("record_id")
        if not name or name in seen:
            raise ValueError("Duplicated/empty contrast")
        seen.add(name)
        kind = r["data_kind"]
        ref = (EXPRESSION if kind == "qRT-PCR mRNA" else
               METABOLITES if kind == "reported metabolite total" else None)
        if ref is None or name not in ref:
            raise ValueError("Unverified molecular or metabolite record")
        gene, numerator, denominator, value=ref[name]
        if (r["molecule"],r["numerator_organ"],r["denominator_organ"]) != (
                gene,numerator,denominator):
            raise ValueError("A quoted fold ratio had its organ or identity changed")
        observed=r["fold_numerator_over_denominator"]
        if not isinstance(observed,(int,float)) or not isfinite(observed) or abs(observed-value)>1e-9:
            raise ValueError("Fold ratio changed from source evidence")
        if numerator == denominator:
            raise ValueError("Expected a comparison of different organs")
        if r.get("citation_doi") != DOI:
            raise ValueError("Missing original evidence provenance")
        if r.get("independent_samples") != "not extractable from text":
            raise ValueError("Never invent biological replication from author summary")
    if seen != (set(EXPRESSION)|set(METABOLITES)):
        raise ValueError("Incomplete or substituted contrast IDs")
    return doc


def report(path):
    doc=validate_registry(path)
    recs=doc["records"]
    groups=Counter((r["data_kind"],r["numerator_organ"],r["denominator_organ"]) for r in recs)
    return {
        "source":doc["source"],
        "real_empirical_source":True,
        "n_author_reported_contrasts":len(recs),
        "grouped_by_kind_and_denominator":[
            {"assay":k[0],"numerator":k[1],"denominator":k[2],"n_distinct_contrasts":v}
            for k,v in sorted(groups.items())],
        "evidence_entries":[
            {"id":r["record_id"],"target":r["molecule"],
             "numerator_organ":r["numerator_organ"],
             "denominator_organ":r["denominator_organ"],
             "published_fold":r["fold_numerator_over_denominator"],
             "log2_fold_descriptive_only":log2(r["fold_numerator_over_denominator"])}
            for r in recs],
        "flower_vs_leaf_flavonoid_mrna_n":sum(
            r["data_kind"] == "qRT-PCR mRNA" and r["numerator_organ"] == "flower"
            and r["denominator_organ"] == "leaf" for r in recs),
        "root_vs_leaf_dfr_mrna_n":sum(
            r["data_kind"] == "qRT-PCR mRNA" and r["molecule"].startswith("DFR_")
            and r["numerator_organ"] == "root" and r["denominator_organ"] == "leaf"
            for r in recs),
        "ryukyu_competence_assessed":False,
        "ancestral_reexpression_demonstrated":False,
        "independent_colonization_history_assessed":False,
        "claim_ceiling":"Real, published organ-specific gene/metabolite allocation in C. japonicum sensu source; no recovered ancestry/pathway-reuse outcome for Ryukyu Cirsium.",
        "comparability_note":"Never compute across-gene mean fold ratios or interpret root:leaf as root:flower.",
    }
