# aza3 v2 — Testing higher-order genetic LEGO in real capitula

**2026-10-19 — Prospective protocol, not completed sampling or adaptation evidence.**

## Why the anto v7 result matters

The [anto theory](https://github.com/zuizui0223/anto/blob/main/docs/JOINT_MODULE_ACCESSIBILITY_AND_DISPERSAL_V7.md) proves a sharp conditional probability bound: for three necessary **functional** resources with each marginal retention 1/2 and each pairwise co-retention 1/4, the probability that all three resources occur **together in one plant** can be anywhere from 0 to 1/4. The first- and second-order moments cannot determine the third-order genotype inventory.

Thus the existing aza3 "LEGO" hypothesis cannot be validated just by repeated head phenotypes, pairwise trait correlations or separate association peaks. We need a **joint genomic/developmental inventory on the same individual** and evidence of compatibility between resources.

## Study domains — capitulum rather than leaves

The primary C. sieboldii system should record the following SAME plants:

| Head module | Direct phenotype measurements | Functional/genomic resource still needed |
|---|---|---|
| Colour | floret visible/UV calibrated colour, developmental onset and pigment chemistry | structural anthocyanin gene/pathway, MYB/bHLH/WD40 floral activation and validated retained capacity |
| Presentation | gravity-referenced capitulum orientation and growth trajectory | independent heritable developmental regulator or resolved causal module |
| Protective envelope | phyllary geometry, posture, involucre stickiness/gland secretion | independent heritable regulatory/biochemical apparatus versus global ancestry |

Head shape and elongation remain cross-module phenotype variables. Leaves may be sampled for pigment pathway pleiotropic retention, **not** treated as the primary capitulum LEGO morphology.

The term "module" has two noninterchangeable meanings. Observed phenotype modules Z can be correlated even if the genetic parts have a common origin. Functional retained resources G require verified biological mechanism. No current 40 C. sieboldii specimens are molecularly admitted to G; all three domains have state UNMEASURED until validated.

## Five nested models — all under aza3 empirical control

M0. Environment, stage, population, cytotype, whole-genome ancestry and per-module marginals.

M1. M0 plus pairwise associations, local ancestry and candidate locus genealogy.

M2. M1 plus independently measured triple co-retention with whole-population or ancestry-held-out evaluation.

M3. M2 plus flower-head developmental expressibility and experimental mediation: can the parts actually build the predicted head?

M4. M3 plus independently controlled ecological function and reproductive fitness, using filled viable achenes as primary outcome.

Three-way predictive improvement is not by itself evidence of ancestral allele reuse; verify local gene trees, standing variants, introgression, homeologous expression and de novo alternatives. Same-colour phenotypes may be retained ancestral states rather than re-expression. Unknown history must stay UNKNOWN.

**Stop conditions:** if only phenotypes exist -> LEGO_ONLY_PHENOTYPE; if G resources not assayed -> UNMEASURED_FUNCTIONAL_MODULE; if one ancestry axis explains co-occurrence -> ONE_ANCESTRY_GRADIENT; if no dated/functional genetic genealogy -> HISTORY_NOT_POLARIZED; if no causal viable-achene outcome -> FITNESS_NOT_DEMONSTRATED.

## Feasibility and sample size

The pre-existing aza3 Phase-A target includes **40 prospective C. sieboldii**, not 40 independent evolutionary origin events. The 298-individual overall Phase-A target, legal gates and prefield plan remain unchanged.

In the idealized case of 40 independent plants and theoretical triple-resource frequency 1/8, the expected count is 5, and P(no triple sample) is (7/8)^40 (~0.48%). With **zero** among 40 independent plants, the two-sided exact 95% binomial upper bound is 1-(0.025)^(1/40), ~8.8%. This cannot be claimed when within-population kinship, cryptic cytotype, linked assay missingness or ascertainment reduce independent effective sample size.

For illustration ONLY, four equal clusters of 10 and an unknown ICC rho give an approximate design effect 1+9 rho and an effective n of 40/(1+9 rho). If rho were 0.2, n_eff would be ~14.3. This is a methodological sensitivity, NOT a measured plant ICC or a formal high-order model power analysis.

A prespecified outcome-blind power/precision gate is mandatory before genomic LEGO, recurrent origin or fitness claims. Whole-population held-out prediction is required where population structure permits; random photographs from the same plant are not independent test records.

## Link to the general theory without duplicating ownership

1. **anto** predicts latent joint resource availability under dispersal, historical degradation, genotype correlation, recombination and regeneration.
2. **aza3** determines whether the resources have distinct origins and whether recombination/genetic sorting creates new head assemblages.
3. **aza3** evaluates environment-specific functional consequences and reproductive success. anto imports a versioned, non-sensitive summary only after validation, not raw empirical Cirsium source.
4. All *Cirsium* source records, pre-field RNA/GWAS plans, optical chemistry, organ-specific reactivation and future sampling are in aza3 (see [empirical master](ANTO_ACCESSIBILITY_TO_AZA3_HEAD_LEGO_EMPIRICAL_MASTER_V1.md)).

A strong negative result is still informative: the visible phenotypic LEGO pattern might not reflect separable genetic resource modules, or actual resources might fail to reassemble developmentally even if they are present.

## Machine and execution contracts

- [Module experiment and missingness requirements](../data/contracts/aza3_capitulum_higher_order_lego_v2.json)
- [Outcome-blind precision and coherence validation](../analysis/validate_aza3_capitulum_higher_order_lego_v2.py)

Reproduce using Python standard library:

    python analysis/validate_aza3_capitulum_higher_order_lego_v2.py

This script uses no actual plant phenotype/genotype outcomes and cannot confer field authorization.
