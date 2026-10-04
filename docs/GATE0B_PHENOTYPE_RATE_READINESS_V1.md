# aza3 Gate 0B — phenotype-rate readiness and primary estimand v1

**Status:** ACTIVE COVERAGE GATE — PRIMARY RATE TEST NOT YET AUTHORIZED  
**Date:** 2026-10-04

## Why Gate 0B exists

Gate 0A supplies the evolutionary-time denominator. Gate 0B supplies the phenotype numerator.

Azami provides broad continuous capitulum phenotypes globally, but exact Japan38 overlap is currently incomplete:

- 38 Japan paper concepts;
- 14 exact Japan38 concepts with trait matches in the strict-spatial Azami cohort;
- 17 paper concepts represented only at binomial level;
- only 6 distinct Japan38 trait taxa with at least 10 strict-spatial observations;
- 18 of 36 Japan38 species binomials are absent from the current exhaustive Azami source pool;
- three additional binomials were detected in the exhaustive source but removed before the strict-spatial cohort.

EAzami provides stronger Japanese coverage for three discrete traits (orientation 20, phyllary 10, stickiness 13), but there is no equally harmonized global ontology for all three traits.

Therefore no Nature-scale comparative phenotype-rate claim is authorized from the current continuous or discrete coverage alone.

## Primary future Gate 0B estimand

The primary comparative rate test will use the high-coverage Azami **core continuous phenotype axes**:

- presentation angle;
- floral lightness;
- floral chroma;
- head elongation;
- head compactness.

Floral hue and candidate involucre/armature constructs are secondary sensitivities because hue is circular and the architectural proxies require stronger botanical calibration.

### Rate statistic

Let sigma2_J,m be the dated phylogenetic evolutionary-rate estimate for core axis m inside the dominant Japanese radiation and sigma2_BG,m the corresponding background estimate from non-Japanese *Cirsium* on the same dated global tree and measurement scale.

The primary total-rate ratio is:

R_total = sum_m sigma2_J,m / sum_m sigma2_BG,m.

Each axis is standardized using the non-Japanese background distribution before rate estimation so that one measurement scale cannot dominate the sum.

### Null and uncertainty

The primary null is R_total = 1.

Significance / interval calibration will use parametric simulation on the **actual dated Japanese subtree**, preserving:

- branch times;
- admitted tip coverage;
- observed missingness pattern;
- measurement-error structure;
- background trait covariance in the declared sensitivity analysis.

The main result is not allowed to depend on selecting the fastest single trait after outcome inspection.

## Measurement uncertainty

Taxon phenotype estimates will be propagated from observation-level data, not treated as error-free species constants.

Preferred implementation:

- retain all admissible strict-spatial observations;
- bootstrap observations within taxon/source strata;
- propagate the resulting taxon-level uncertainty jointly with dated-tree uncertainty;
- preserve exact Japan38 concept identity and never assign broad binomial images to named varieties/subspecies without source-backed reconciliation.

## Coverage admission rule

No fixed tip-count threshold is chosen from the phenotype outcome.

After Gate 0A provides the dated tree but **before trait-rate outcomes are opened**, an outcome-blind simulation will test whether the observed/recovered Japan coverage can estimate R_total with useful precision under R_total = 1 and under a predeclared acceleration alternative.

If the rate multiplier is not identifiable at the available coverage, the analysis remains `COVERAGE_NOT_IDENTIFIABLE` and more trait recovery is required.

## Role of the EAzami discrete traits

Orientation, phyllary posture and stickiness remain an independent historical layer.

After a dated ensemble exists, they can supply secondary summaries such as:

- state-change accumulation per unit branch time;
- trait-specific CTMC rate estimates;
- timing of reconstructed changes;
- consistency with the continuous-rate result where phenotypes are homologous.

They do not substitute for the global continuous comparator because homologous global state coding is not currently available for all three modules.

## Current first recovery actions

1. Recover phenotype data for the 18 source-pool gaps using a frozen synonym-aware public-image/herbarium protocol.
2. Recover the three taxa present in the exhaustive source but removed only by strict spatial filtering for non-spatial trait-rate summaries under a separate sensitivity.
3. Prioritize exact nuclear tips and internal-branch coverage rather than maximizing raw photo count.
4. Resolve identity-conflicted taxa before any trait reuse.
5. Keep field collection as a later option if public/herbarium recovery cannot make the rate estimand identifiable.

## Stop rules

- Do not use broad species images as exact infraspecific phenotypes.
- Do not lower observation/identity standards after seeing rate results.
- Do not claim acceleration from the six >=10-observation Japanese taxa alone.
- Do not use candidate armature/image proxies as equivalent to botanical phyllary/spine states.
- Do not promote a single fast axis into a whole-capitulum evolvability claim.
- Do not call Gate 0B passed unless the outcome-blind precision simulation and dated-tree gate both pass.
