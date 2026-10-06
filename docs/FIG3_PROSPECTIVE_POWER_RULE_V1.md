# aza3 Fig.3 prospective power rule v1

Status: frozen before focal phenotype-linked genomic discovery.

## What is powered

The Nature-scale discriminator is not a single-SNP association.

The powered event is the frozen E1 result:

R (reconfigurable inheritance) beats M (frozen phenotypic-module inheritance) in untouched validation individuals under the paired one-sided sign-flip test.

For validation individual i:

d_RM_i = log p_R(y_i | genotype_i, frozen TRAIN fit) - log p_M(y_i | genotype_i, frozen TRAIN fit)

The power problem is therefore defined on the distribution of d_RM, not on an assumed GWAS odds ratio.

## Quantities that may come from the 8-individual pilot

Only phenotype-blind / genotype-only quantities may inform the final sample-size simulation:

- usable-read / callable fraction;
- heterozygosity / MAF spectrum;
- LD decay and effective genomic-block count;
- genome-wide relatedness;
- population structure;
- expected post-QC retention fraction;
- reference-dependent missingness measured before phenotype labels open.

The 8 pilot individuals cannot contribute phenotype-linked effect estimates.

## Frozen predictive-separation grid

Power is reported over a grid of standardized held-out predictive separation:

delta = E[d_RM] / SD(d_RM)

Use the frozen grid:

0.20, 0.30, 0.40, 0.50, 0.70

Interpretation:

- 0.20: weak but systematic predictive advantage;
- 0.30: modest;
- 0.40: moderate;
- 0.50: substantial;
- 0.70: large.

These labels are descriptive only. No one grid point is declared the true effect before data exist.

Primary design requirement:

- at least 80% power at delta = 0.40;
- report power at every frozen grid point;
- target 90% power at delta = 0.50 as a secondary design benchmark.

If feasible sampling cannot meet the primary design requirement after accounting for expected QC retention, the strongest Nature E1 route is prospectively classified as underpowered / NOT_IDENTIFIABLE rather than relaxed after seeing outcomes.

## Validation n grid

Evaluate untouched validation n:

20, 24, 30, 36, 40, 50, 60

The minimum remains 20, but 20 is a design floor, not a promise of adequate power.

## Training n

TRAIN must contain at least 40 usable high-confidence focal individuals.

After OWN_WGS_GREEN, TRAIN n may be increased using genotype-only estimates of effective block count, LD and expected model complexity.

No phenotype-linked genomic result may be used to decide the increase.

The power report therefore has two distinct checks:

1. validation discrimination power on the frozen d_RM grid;
2. genotype-only adequacy of TRAIN relative to the effective genomic dimension and frozen F/M/R complexity budget.

Failure of either check closes the strongest route.

## Simulation model

For the validation-discrimination calculation, simulate standardized paired score differences:

d_i ~ Normal(delta, 1)

and apply the exact same one-sided sign-flip decision rule used by the frozen scorer.

This is deliberately an estimand-level design calculation. It does not pretend that the full genomic training pipeline is known before focal data exist.

A second-stage end-to-end simulation may be added after the genotype-only pilot, but it must:

- use only genotype-only pilot quantities;
- preserve the frozen F/M/R fitting procedure;
- preserve TRAIN/VALIDATION separation;
- never use validation phenotypes for tuning.

## QC inflation

Let r be the expected usable-individual fraction estimated without phenotype labels.

Required collected n for a desired usable n is:

ceil(usable_n / r)

Report r and the resulting inflation explicitly.

Do not silently replace failed or excluded validation individuals after validation phenotype labels are opened.

## Decision output

POWER_READY:
- validation grid includes a feasible n with >=80% power at delta=0.40;
- TRAIN adequacy can be met;
- expected QC-inflated collection target is feasible.

POWER_LIMITED:
- only larger effects such as delta>=0.50 are adequately powered at feasible n.

POWER_NOT_IDENTIFIABLE:
- pilot quantities needed for QC inflation / effective genomic dimension are not yet estimable.

POWER_INFEASIBLE:
- the frozen primary power requirement cannot be met within the feasible sampling ceiling.

No power class is evidence for F, M or R.
