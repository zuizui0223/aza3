#!/usr/bin/env python3
"""Evidence-and-join audit for three linked Cirsium evolutionary questions.

Audit REAL existing public accession metadata and the AUTHORITATIVE aza3
sampling targets, without importing a single RNA read, phenotype or model
simulation as if it were a proven floral molecular resource.

Q1: Phylogenetic memory has history data (EAzami), but independent
    cross-lineage prediction on linked individual head phenotypes is open.
Q2: Partner-specific physical access -> seed benefit has no measured
    head-level observations in the pinned native EAzami intake.
Q3: Ryukyu/takaoense LEAF RNA metadata are real, but no linked floral
    enzyme competence, induced re-expression or heritable regain.

Source accessions and planning denominators must not be interchanged.
"""
from __future__ import annotations
import argparse
import collections
import copy
import csv
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
RUNS=ROOT/"data/anto_source_handoff/public_sra/prjna1311153_run_manifest_20261017.csv"
MORPHS=ROOT/"data/anto_source_handoff/public_sra/takaoense_aza3_morph_crosswalk_20261017.csv"
GENES=ROOT/"data/anto_source_handoff/molecular/ryukyu_gene_screen_register_v5.csv"
PANEL=ROOT/"data/planning/aza3_nature_genomic_sampling_panel_v1.csv"
BRIDGE=ROOT/"data/contracts/aza3_anto_empirical_lego_bridge_v1.json"
LEGO=ROOT/"data/contracts/aza3_capitulum_higher_order_lego_v2.json"

FOCAL=("Cirsium brevicaule","Cirsium irumtiense",
       "Cirsium japonicum var. takaoense")
RUN_COLS={"study_accession","run_accession","experiment_accession",
          "biosample_accession","taxon_as_experiment_titled",
          "specimen_code_from_experiment_title","tissue_from_project_description",
          "individual_flower_colour","library_strategy"}
MORPH_COLS={"run_accession","biosample_accession",
            "specimen_code_from_experiment_title","morph_code_from_aza3",
            "laboratory_tissue","floral_rna_measured",
            "floral_pigment_chemistry_measured",
            "independent_population_replication"}
GENE_COLS={"run_accession","biosample_accession","source_specimen_code",
           "taxon_title","screen_locus","candidate_accession","status",
           "floral_expression_measured","floral_pigment_measured",
           "orthology_verified","genomic_coding_variant_verified",
           "independent_event_origin_inferred"}
PANEL_COLS={"panel_id","phase","taxon","recommended_n"}

def read_csv(path):
    with Path(path).open(encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f)
        return set(reader.fieldnames or ()),list(reader)

def require(condition,code):
    if not condition:
        raise ValueError(code)

def index_unique(rows,col,label):
    result={}
    for row in rows:
        val=(row.get(col) or "").strip()
        require(bool(val),label+"_EMPTY_KEY")
        require(val not in result,label+"_DUPLICATE_KEY")
        result[val]=row
    return result

def checked_bout_grain(heads,bouts):
    if heads is None or bouts is None:
        return {"status":"ORIGINAL_EAZAMI_INTAKES_NOT_EXAMINED",
                "measured_heads":None,"observed_bouts":None}
    # The grain is exactly population × individual × capitulum
    # (not globally unique local head labels); bouts additionally need ID.
    # Otherwise unrelated plants with a shared H1 label would be conflated.
    head_keys=set()
    for head in heads:
        key=tuple((head.get(k) or "").strip() for k in
                  ("population_id","individual_id","capitulum_id"))
        require(all(key),"HEAD_COMPOSITE_ID_EMPTY")
        require(key not in head_keys,"HEAD_COMPOSITE_ID_DUPLICATE")
        head_keys.add(key)
    bout_keys=set()
    for bout in bouts:
        key=tuple((bout.get(k) or "").strip() for k in
                  ("population_id","individual_id","capitulum_id",
                   "observation_bout_id"))
        require(all(key),"BOUT_COMPOSITE_ID_EMPTY")
        require(key not in bout_keys,"BOUT_COMPOSITE_ID_DUPLICATE")
        require(key[:3] in head_keys,
                "OBSERVATION_BOUT_WITHOUT_ENROLLED_CAPITULUM")
        bout_keys.add(key)
    return {
        "status":("PINNED_EAZAMI_ORIGINAL_FIELD_INTAKE_EMPTY" if
                  not heads and not bouts else
                  "SOURCE_FIELD_ROWS_PRESENT_REQUIRE_LINKED_OUTCOME_AUDIT"),
        "measured_heads":len(heads),
        "observed_bouts":len(bouts),
        "independent_functional_seed_or_pollen_receipts":0
            if not heads and not bouts else None,
    }

def audit(runs,morphs,genes,panel,bridge,lego,original_heads=None,
          original_bouts=None):
    run_idx=index_unique(runs,"run_accession","ENA_RUN")
    bio_idx=index_unique(runs,"biosample_accession","ENA_BIOSAMPLE")
    code_idx=index_unique(runs,"specimen_code_from_experiment_title",
                          "ENA_SOURCE_SPECIMEN_CODE")
    require(len(runs)==25,"FROZEN_ENA_MANIFEST_NOT_25")
    require(len(bio_idx)==len(code_idx)==25,"SOURCE_SPECIMEN_GRAIN_DRIFT")
    require(all(row.get("study_accession")=="PRJNA1311153"
                for row in runs),"MULTIPLE_UNRELATED_BIOPROJECTS")
    require(all(row.get("tissue_from_project_description")=="young leaf"
                and row.get("individual_flower_colour")=="NOT_MEASURED_IN_THIS_DATASET"
                and row.get("library_strategy")=="RNA-Seq"
                for row in runs),"LEAF_RNA_WAS_FALSELY_DECLARED_FLORAL_ASSAY")
    taxa=collections.Counter(row["taxon_as_experiment_titled"] for row in runs)
    require(len(taxa)==8,"NOT_THE_8_DECLARED_TAXON_LABELS")
    require([taxa[x] for x in FOCAL]==[3,3,6],
            "RYUKYU_OR_TAKAOENSE_PUBLIC_RNA_SAMPLE_COUNT_DRIFT")

    morph_idx=index_unique(morphs,"run_accession","MORPH_CROSSWALK")
    require(len(morph_idx)==6 and
            set(morph_idx)=={
                row["run_accession"] for row in runs
                if row["taxon_as_experiment_titled"]==FOCAL[2]},
            "TAKAOENSE_RNA_CROSSWALK_DROPS_OR_SWAPS_SOURCE_RUNS")
    morph_counts=collections.Counter()
    for acc,row in morph_idx.items():
        parent=run_idx[acc]
        require(row["biosample_accession"]==parent["biosample_accession"]
                and row["specimen_code_from_experiment_title"]==
                parent["specimen_code_from_experiment_title"],
                "SAMPLE_BIOSAMPLE_OR_SPECIMEN_CROSSWALK_MISMATCH")
        require(row["laboratory_tissue"]=="young leaf"
                and row["floral_rna_measured"]=="false"
                and row["floral_pigment_chemistry_measured"]=="false"
                and row["independent_population_replication"]=="false",
                "TAKAOENSE_FLORAL_FUNCTION_OR_REPLICATION_FALSE_PROMOTION")
        morph=row["morph_code_from_aza3"]
        require(morph in {"W","BP"},"UNREGISTERED_COLOUR_MORPH_CODE")
        morph_counts[morph]+=1
    require(morph_counts=={"W":3,"BP":3},
            "TAKAOENSE_PUBLIC_THREE_THREE_MORPH_BALANCE_DRIFT")

    gene_index=set()
    for g in genes:
        k=(g.get("run_accession"),g.get("screen_locus"))
        require(k not in gene_index,"DUPLICATE_RECONSTRUCTED_MOLECULAR_LOCUS")
        gene_index.add(k)
        accession=g["run_accession"]
        require(accession in run_idx and
                run_idx[accession]["taxon_as_experiment_titled"] in FOCAL[:2],
                "RYUKYU_GENE_GATEMATRIX_FOREIGN_RUN")
        parent=run_idx[accession]
        require(g["biosample_accession"]==parent["biosample_accession"]
                and g["source_specimen_code"]==
                parent["specimen_code_from_experiment_title"]
                and g["taxon_title"]==parent["taxon_as_experiment_titled"],
                "GENE_GATEMATRIX_BIOSAMPLE_TAXON_MISMATCH")
        require(g["screen_locus"] in {"ANS","DFR","GMYC1"},
                "UNEXPECTED_SCREEN_LOCI")
        require(g["candidate_accession"]=="NOT_RETRIEVED"
                and g["status"]=="NOT_TESTED",
                "NOT_TESTED_GENE_WAS_PROMOTED_TO_MEASURED")
        for field in ("floral_expression_measured","floral_pigment_measured",
                      "orthology_verified","genomic_coding_variant_verified",
                      "independent_event_origin_inferred"):
            require(g[field]=="false",
                    "ABSENT_FUNCTIONAL_OR_ANCESTRAL_ASSAY_WAS_FABRICATED")
    ryukyu={row["run_accession"] for row in runs if
            row["taxon_as_experiment_titled"] in FOCAL[:2]}
    require(gene_index=={(run,locus) for run in ryukyu
                         for locus in ("ANS","DFR","GMYC1")},
            "RYUKYU_GENE_LABELLING_HAS_MISSING_OR_EXTRA_SOURCE_CELLS")
    require(len(genes)==18,"NOT_EIGHTEEN_UNTESTED_SOURCE_CELLS")

    phase_a={}
    for row in panel:
        if row.get("phase")=="A":
            taxon=row["taxon"]
            require(taxon not in phase_a,"PANEL_A_DUPLICATE_TAXON")
            try:n=int(row["recommended_n"])
            except ValueError:
                raise ValueError("PANEL_A_N_NOT_INTEGER") from None
            require(n>0,"PANEL_A_NONPOSITIVE_TARGET")
            phase_a[taxon]=n
    expected={
        x["taxon"]:x["planned_n"] for x in
        bridge["samplings"]["by_taxon_target"]
    }
    require(phase_a==expected and len(phase_a)==6,
            "PHASE_A_SOURCE_AND_MASTER_BRIDGE_DISAGREE")
    require(sum(phase_a.values())==bridge["samplings"]["phase_A_original_total"]==298,
            "PHASE_A_298_TARGET_NOT_CONSISTENT")
    require(phase_a["Cirsium sieboldii"]==
            lego["study_design"]["planned_phase_A_total"]==40,
            "SIEBOLDII_40_TARGET_VS_COHORT_GATE_DRIFT")
    require(bridge["samplings"]["phase_A_collected_n_currently_admitted"]==0
            and lego["explicit_claim_ceilings"]["actual_collected_individuals_in_this_new_protocol"]==0,
            "PROSPECTIVE_SAMPLE_TARGET_FALSELY_PROMOTED_TO_COLLECTED")

    field=checked_bout_grain(original_heads,original_bouts)
    return {
        "status":"SOURCE_AUDITED_THREE_QUESTION_EMPIRICAL_CHAIN_INCOMPLETE",
        "sources":{
            "ENA_RNA_run_metadata_rows":len(runs),
            "ENA_distinct_biosamples":len(bio_idx),
            "ENA_taxon_title_count":len(taxa),
            "ENA_taxon_title_counts":dict(sorted(taxa.items())),
            "ryukyu_leaf_runs_by_species":{x:taxa[x] for x in FOCAL[:2]},
            "takaoense_leaf_RNA_w":morph_counts["W"],
            "takaoense_leaf_RNA_bp":morph_counts["BP"],
            "ryukyu_3_locus_unretrieved_orthology_screen_cells":len(genes),
            "validated_functional_Ryukyu_ANS_DFR_GMYC1_screen_cells":0,
            "phase_A_planned_individuals_not_collected":sum(phase_a.values()),
            "phase_A_target_by_taxon":dict(sorted(phase_a.items())),
            "planned_C_sieboldii_joint_module_individuals_not_collected":
                phase_a["Cirsium sieboldii"],
        },
        "field_source":field,
        "claim_gates":{
            "Q1_lineage_trait_memory":{
                "source_history":"EAzami repeated phenotypes and turnover/chun original results belong to their independent repositories",
                "predictive_Cirsium_lineage_holdout_completed":False,
                "reason":"RNA metadata and six Phase-A targets are not homologous trait-state matrices scored across independent Cirsium radiations",
                "maximum_claim":"EXTERNAL_EVOLUTIONARY_HISTORY_STARTING_POINT_ONLY",
            },
            "Q2_partner_conditional_payoff":{
                "actual_verified_Cirsium_anatomical_contact_and_functional_seed_joins":0
                    if field["measured_heads"]==0 else None,
                "status":"NO_NEW_CIRSIUM_PARTNER_CAUSAL_TRAIT_FUNCTION_REPORTED",
                "reason":"Source RNA is nonfloral and original EAzami native head/bout intake has no real rows when pinned and inspected",
                "maximum_claim":"PROSPECTIVE_EXPERIMENT_NOT_NATURAL_SELECTION_RESULT",
            },
            "Q3_after_phenotype_loss":{
                "public_takaoense_W_BP_metadata_pair":True,
                "public_young_leaf_RNA_source_available_for_genealogical_study":True,\n                "RNA_metadata_alone_proves_genealogy":False,
                "corolla_flower_RNA_linked_to_same_source_individual":False,
                "validated_enzymatic_toolkit_retained_after_loss":False,
                "verified_ancestral_loss_then_same_lineage_reactivation":False,
                "heritable_multigeneration_return_assayed":False,
                "maximum_claim":"ANATOMICAL_AND_LEAF_RNA_MOTIVATION_NOT_REEXPRESSION",
            },
            "Q1_x_Q3_history_mediated_by_latent_state":"UNIDENTIFIABLE_NO_LATER_EXPRESSION_OR_FUNCTIONAL_RESOURCE_STATE",
            "Q2_x_Q3_selection_value_vs_resource_accessibility":"UNIDENTIFIABLE_NO_CO_MEASURED_INTERVENTION_OR_SEGREGATION",
        },
        "one_minimum_progress_package":{
            "unit":"same pre-enrolled plant and its capitulum with nested independent approach events",
            "first_needed":"sample registry + voucher/authorization + true head morphology + flower-organ pigment and flower transcript + cytotype/ancestry crosswalk",
            "functional_followup":"concurrent independently identified pollinator/seed feeder approach, true barrier contact, independently assayed oviposition/pollen and per-registered-head filled achenes",
            "history_followup":"independent gene genealogies, alternative introgression/de_novo histories, held-out populations/lineages and heritable expression assays",
            "not_a_claim_of_power":"No sample size guarantees identification of triple functional resources, selection or independent evolutionary reversals",
        },
        "inferential_limit":"PUBLIC_LEAF_LIBRARY_COUNTS_ARE_NOT_TISSUE_SPECIFIC_FLORAL_CAPACITY_OR_NEW_FIELD_FITNESS",
    }

def adversarial_tests(runs,morphs,genes,panel,bridge,lego):
    def will_fail(a,b,c,d,needle):
        try:audit(a,b,c,d,bridge,lego,[],[])
        except ValueError as err:
            require(needle in str(err),
                    "UNEXPECTED_NEGATIVE_CONTROL_ERROR:"+str(err))
            return
        raise AssertionError("FALSE_BIOLOGICAL_EVIDENCE_ACCEPTED:"+needle)
    base=audit(runs,morphs,genes,panel,bridge,lego,[],[])
    require(base["field_source"]["status"]=="PINNED_EAZAMI_ORIGINAL_FIELD_INTAKE_EMPTY",
            "SYNTHETIC_EMPTY_FIELD_NOT_EXPLICIT")
    x=copy.deepcopy(runs);x.append(copy.deepcopy(runs[0]))
    will_fail(x,morphs,genes,panel,"ENA_RUN_DUPLICATE_KEY")
    x=copy.deepcopy(runs);x[0]["tissue_from_project_description"]="corolla"
    will_fail(x,morphs,genes,panel,"LEAF_RNA_WAS_FALSELY_DECLARED_FLORAL_ASSAY")
    x=copy.deepcopy(morphs);x[0]["biosample_accession"]=runs[0]["biosample_accession"]
    will_fail(runs,x,genes,panel,"SAMPLE_BIOSAMPLE_OR_SPECIMEN_CROSSWALK_MISMATCH")
    x=copy.deepcopy(morphs);x[0]["floral_rna_measured"]="true"
    will_fail(runs,x,genes,panel,"TAKAOENSE_FLORAL_FUNCTION_OR_REPLICATION_FALSE_PROMOTION")
    x=copy.deepcopy(genes);x[0]["status"]="CONFIRMED_ENZYME"
    will_fail(runs,morphs,x,panel,"NOT_TESTED_GENE_WAS_PROMOTED_TO_MEASURED")
    x=copy.deepcopy(genes);x[0]["taxon_title"]="Cirsium brevicaule"
    if genes[0]["taxon_title"]=="Cirsium brevicaule":
        x[0]["taxon_title"]="Cirsium irumtiense"
    will_fail(runs,morphs,x,panel,"GENE_GATEMATRIX_BIOSAMPLE_TAXON_MISMATCH")
    x=copy.deepcopy(panel)
    next(a for a in x if a["panel_id"]=="SIE")["recommended_n"]="400"
    will_fail(runs,morphs,genes,x,"PHASE_A_SOURCE_AND_MASTER_BRIDGE_DISAGREE")
    return {
        "status":"PASS_FROZEN_CIRSIUM_LEAF_SOURCE_AND_NO_FALSE_REEXPRESSION",
        "source_backed_real_RNA_run_rows_tested":25,
        "synthetic_corruptions_rejected":7,
        "actual_new_flower_or_fitness_measurements_added":0,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-eazami-heads",type=Path)
    p.add_argument("--original-eazami-bouts",type=Path)
    p.add_argument("--out",type=Path)
    x=p.parse_args()
    require((x.original_eazami_heads is None)==
            (x.original_eazami_bouts is None),
            "ORIGINAL_EAZAMI_INTAKE_TWO_SOURCE_FILES_REQUIRED_TOGETHER")
    manifests=[]
    for source,columns in [(RUNS,RUN_COLS),(MORPHS,MORPH_COLS),
                           (GENES,GENE_COLS),(PANEL,PANEL_COLS)]:
        header,values=read_csv(source)
        require(columns.issubset(header),
                "ORIGINAL_CIRSIUM_SOURCE_SCHEMA_CHANGED:"+source.name)
        manifests.append(values)
    runs,morphs,genes,panel=manifests
    bridge=json.loads(BRIDGE.read_text(encoding="utf-8"))
    lego=json.loads(LEGO.read_text(encoding="utf-8"))
    if x.original_eazami_heads:
        hhead,heads=read_csv(x.original_eazami_heads)
        bhead,bouts=read_csv(x.original_eazami_bouts)
        require({"individual_id","population_id","capitulum_id"}.issubset(hhead)
                and {"individual_id","population_id","capitulum_id",
                     "observation_bout_id"}.issubset(bhead),
                "ORIGINAL_EAZAMI_FIELD_HEADER_SOURCE_DRIFT")
    else:
        heads=bouts=None
    result={
        "version":"aza3_three_question_source_coverage_v1",
        "source": {
            "ENA_project":"PRJNA1311153",
            "repository_base":"zuizui0223/aza3",
            "sample_panel":"Phase-A n=298 prospective, not observed",
            "original_EAzami":"Frozen main c9870e91604ed1514344ddc8f765eae5237dbf27 when supplied",
        },
        "negative_controls":adversarial_tests(
            runs,morphs,genes,panel,bridge,lego),
        "actual_evidence":audit(runs,morphs,genes,panel,bridge,lego,heads,bouts),
    }
    if x.out:
        x.out.parent.mkdir(parents=True,exist_ok=True)
        x.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",
                         encoding="utf-8")
    print(json.dumps({
        "source_runs":result["actual_evidence"]["sources"]["ENA_RNA_run_metadata_rows"],
        "morph_libraries":result["actual_evidence"]["sources"]["takaoense_leaf_RNA_w"]+
                          result["actual_evidence"]["sources"]["takaoense_leaf_RNA_bp"],
        "Q3_functional_loci":result["actual_evidence"]["sources"]["validated_functional_Ryukyu_ANS_DFR_GMYC1_screen_cells"],
        "original_EAzami_head_rows":result["actual_evidence"]["field_source"]["measured_heads"],
        "design_safety":result["negative_controls"]["status"],
    },ensure_ascii=False,indent=2))

if __name__=="__main__":
    main()
