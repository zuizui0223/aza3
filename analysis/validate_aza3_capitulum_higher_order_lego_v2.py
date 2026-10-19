"""Prospective *Cirsium* LEGO higher-order measurement and precision gates.

The cross-plant theory is anto v7. All Cirsium biological measurements,
same-individual inheritance tests and trait-fitness conclusions belong to
aza3. This is an outcome-blind DESIGN audit only: there are no new genomic
or floral functional observations.

Read:
    data/contracts/aza3_capitulum_higher_order_lego_v2.json
    data/contracts/aza3_anto_empirical_lego_bridge_v1.json

The 40 target C. sieboldii plants are NOT 40 independent trait re-evolutions,
and are not yet genotyped for functional module inventories. Phenotypic
colour/angle/stickiness measurements cannot substitute for those inventories.
"""
from fractions import Fraction
from pathlib import Path
from math import exp,log,sqrt
import json

ROOT=Path(__file__).resolve().parents[1]
CONTRACT=ROOT/"data/contracts/aza3_capitulum_higher_order_lego_v2.json"
BRIDGE=ROOT/"data/contracts/aza3_anto_empirical_lego_bridge_v1.json"
MOMENTS={
    "marginal":(Fraction(1,2),)*3,
    "pairs":{(0,1):Fraction(1,4),(0,2):Fraction(1,4),(1,2):Fraction(1,4)}
}


def feasible_triple_bounds(m,p):
    """Elementary, exact three-event Fréchet moment constraints.

    The code is deliberately independent of anto's implementation so that
    aza3 validates the planned assay target without cross-repo imports.
    """
    if len(m)!=3 or set(p)!={(0,1),(0,2),(1,2)}:
        raise ValueError("Three marginal and three pairwise joint moments required")
    m1,m2,m3=(Fraction(x) for x in m)
    p12,p13,p23=(Fraction(p[x]) for x in ((0,1),(0,2),(1,2)))
    if any(x<0 or x>1 for x in (m1,m2,m3,p12,p13,p23)):
        raise ValueError("Probabilities must be in [0,1]")
    low=max(Fraction(0),p12+p13-m1,p12+p23-m2,p13+p23-m3)
    high=min(p12,p13,p23,1-m1-m2-m3+p12+p13+p23)
    if low>high:raise ValueError("Inconsistent first and second moments")
    return low,high


def no_complete_detection_probability(p,n):
    """Exact binomial null conditional on INDEPENDENT founders/individuals."""
    if not 0<=p<=1 or int(n)!=n or n<1:
        raise ValueError("Valid positive sample count and event probability required")
    return (1-p)**n


def design_effect(n,individuals_per_population,intraclass_correlation):
    """Kish-style sensitivity approximation under equal cluster sizes.

    THIS IS NOT a general estimator for genotype LD or higher-order
    interaction power, and unmeasured ICC must never be copied as a fact.
    """
    if n<1 or individuals_per_population<1 or n%individuals_per_population:
        raise ValueError("Equal-size illustrative cluster design required")
    if not 0<=intraclass_correlation<=1:
        raise ValueError("Intraclass correlation must be between zero and one")
    deff=1+(individuals_per_population-1)*intraclass_correlation
    return {"design_effect":deff,"illustrative_effective_n":n/deff}


def exact_zero_success_two_sided_upper(n,alpha=.05):
    """Clopper-Pearson upper confidence bound when 0 of n Bernoulli successes.

    All observations must be independent, otherwise this upper confidence
    limit is overprecise and cannot be interpreted as effective n.
    """
    if n<1 or not 0<alpha<1:raise ValueError("Invalid count or error rate")
    return 1-(alpha/2)**(1/n)


def run():
    contract=json.loads(CONTRACT.read_text(encoding="utf-8"))
    bridge=json.loads(BRIDGE.read_text(encoding="utf-8"))
    assert contract["version"]=="aza3_capitulum_higher_order_lego_v2"
    assert contract["repositories"]=={
        "empirical_owner":"zuizui0223/aza3",
        "upstream_general_theory":"zuizui0223/anto",
        "theoretical_result":"anto/docs/JOINT_MODULE_ACCESSIBILITY_AND_DISPERSAL_V7.md"
    }
    assert bridge["canonical_empirical_repository"]=="zuizui0223/aza3"
    assert bridge["samplings"]["phase_A_original_total"]==298
    assert bridge["samplings"]["phase_A_collected_n_currently_admitted"]==0
    planned=next(x["planned_n"] for x in bridge["samplings"]["by_taxon_target"]
                 if x["taxon"]=="Cirsium sieboldii")
    assert planned==contract["study_design"]["planned_phase_A_total"]==40
    assert len(contract["study_design"]["three_module_theoretical_mapping"])==3
    assert all(item["functional_inventory_status"]=="UNMEASURED"
               and item["not_equal_to_latent_inventory"] is True
               for item in contract["study_design"]["three_module_theoretical_mapping"])
    assert len(contract["model_comparison"])==5
    assert [m["id"] for m in contract["model_comparison"]]==["M0","M1","M2","M3","M4"]
    assert contract["data_grain"]["joint_assay_required"]=="SAME_GENOTYPED_BIOLOGICAL_INDIVIDUAL"
    assert contract["power_precision_precondition"]["n_eff_UNKNOWN"] is True
    assert contract["explicit_claim_ceilings"]["actual_collected_individuals_in_this_new_protocol"]==0

    low,high=feasible_triple_bounds(MOMENTS["marginal"],MOMENTS["pairs"])
    assert (low,high)==(Fraction(0),Fraction(1,4))
    n=planned
    scenarios=[]
    for p in contract["power_precision_precondition"]["prospective_sensitivity_scenarios"]["p_complete"]:
        scenarios.append({
            "hypothetical_p_complete":p,
            "expected_complete_inventory_count_if_n_independent":n*p,
            "probability_zero_complete_in_40_independent_draws":
                no_complete_detection_probability(p,n),
            "not_fitted":True,
        })
    cluster_size=contract["power_precision_precondition"]["prospective_sensitivity_scenarios"]["hypothetical_cluster_sizes"][0]
    deff=[]
    for rho in contract["power_precision_precondition"]["prospective_sensitivity_scenarios"]["illustrative_ICC"]:
        x=design_effect(n,cluster_size,rho)
        deff.append({
            "assumed_ICC":rho,
            "illustrative_clusters":n//cluster_size,
            "illustrative_mean_cluster_size":cluster_size,
            **x,
            "warning":"Effective-n is a heuristic only, not a power calculation or observed population ICC"
        })
    return {
        "version":contract["version"],
        "study_population":"Cirsium sieboldii; planned Phase-A n=40, not sampled for this analysis",
        "registered_individuals_qualified_for_joint_resource_inference":0,
        "all_three_functional_resources_measured":False,
        "formal_anto_moment_bound":[str(low),str(high)],
        "identical_marginal_and_pairwise_example":{
            "M1_retention":str(MOMENTS["marginal"][0]),
            "each_pair_co_retention":str(MOMENTS["pairs"][(0,1)]),
            "triple_retention_nonidentifiable_interval":[str(low),str(high)]
        },
        "hypothetical_n40_draw_scenarios":scenarios,
        "zero_successes_upper95_if_all_40_independent":
            exact_zero_success_two_sided_upper(n),
        "cluster_assumption_sensitivity":deff,
        "real_empirical_power_study_completed":False,
        "new_biological_data_admitted":False,
        "proven_genetic_LEGO_reuse":False,
        "validated_Ryukyu_flower_regain":False,
        "fitness_adaptation_demonstrated":False,
        "claim_ceiling":"Prospective multi-module within-species design only; measure provenance, functional resources, whole ancestry and fitness before inference."
    }


if __name__=="__main__":
    print(json.dumps(run(),indent=2,ensure_ascii=False))
