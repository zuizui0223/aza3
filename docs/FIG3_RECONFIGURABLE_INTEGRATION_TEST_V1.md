# aza3 Fig.3–4 contract v1 — fixed versus reconfigurable integration

Status: prospective contract; unopened until Gate CS and own-WGS reference transferability pass
Question: Must an integrated complex organ evolve as an integrated unit?

## Biological models

### Model F — fixed evolutionary integration

The focal capitulum components may vary phenotypically, but their heritable variation is carried by one shared genomic state.

Operational prediction:
- the same local genomic blocks / haplotype state explain the primary component traits;
- component effects are proportional to a common genomic axis (rank-1/shared-state architecture);
- a whole-head genomic score predicts the joint phenotype as well as or better than trait-specific scores;
- apparent recombinant phenotype combinations do not require discordant component-specific genomic states.

### Model R — reconfigurable integration

Present organ-level integration is reconstructed from component-specific genomic states that can assort independently.

Operational prediction:
- different component traits have reference-stable genomic predictors that need not occupy the same blocks or local histories;
- a componentwise genomic model predicts held-out joint phenotypes better than a forced shared-state model;
- natural individuals carry discordant component-specific genomic states (genomic mosaics);
- those mosaics prospectively predict the component combination carried by the individual.

The decisive test is prediction, not merely non-overlapping GWAS peaks.

## Gate prerequisites

The test is unopened unless:

1. Gate CS is classified.
   - CS_GREEN: strongest within-population route.
   - CS_GREEN_LOCAL_SET: local-set route allowed, but within-population natural-reassembly wording remains closed.
   - CS_AMBER / CS_RED_* / NOT_IDENTIFIABLE: do not force C. sieboldii into Fig.3–4.
2. The own-WGS reference-transferability pilot is OWN_WGS_GREEN.
3. Every included individual has same-individual phenotype, taxonomic-confidence, 2C/cytotype-QC and genomic records.
4. The primary trait pair is frozen before genomic discovery.

## Trait hierarchy

Primary:
- floral colour;
- A1–A2 anthesis orientation.

Secondary:
- phyllary/involucre state only if Gate CS prospectively opens it with sufficient individual variation.

Exploratory:
- stickiness.

Stickiness cannot rescue a failed primary test.

## Sample-size boundary

The 8-individual own-WGS pilot is technical only and cannot be used to discover phenotype-linked genomic regions.

The legacy 40-individual C. sieboldii allocation is not automatically Nature-ready.

For the prospective prediction route:
- training set: target at least 40 high-confidence focal individuals;
- untouched validation set: at least 20 individuals;
- preferred validation: a second independently segregating population;
- if no second population exists, a prospectively frozen spatially blocked holdout within the focal population is allowed, but the claim ceiling becomes within-population prediction only.

Thus the strongest route requires at least 60 usable individuals after QC.

If fewer than 20 untouched validation individuals remain, Fig.3 discovery may proceed but Fig.4 prospective-combination prediction is closed.

The final total n may be increased after the phenotype-blind WGS pilot using genotype-only callability, LD and allele-frequency information. No phenotype-linked genomic result may be used to justify that increase.

## Frozen data split

Before phenotype-linked genomic analysis:
1. freeze population/local-set membership;
2. designate TRAIN and VALIDATION individuals;
3. write the split and individual IDs to the registry;
4. hide validation phenotype labels from model fitting;
5. do not change the split after candidate genomic blocks are seen.

Priority: independent segregating population > predeclared local population > spatially blocked within-population holdout.

## Common covariates

Both models receive the same non-focal covariates:
- global ancestry / genome-wide structure;
- population or spatial block as appropriate;
- 2C/genome-size QC;
- cytotype where informative;
- sequencing/callability metrics required by R1B;
- batch variables frozen before phenotype-linked analysis.

Model R is not allowed extra demographic covariates unavailable to Model F.

## Training-stage genomic discovery

Candidate blocks are discovered in TRAIN only.

A block can enter the frozen candidate set only if:
- phenotype association / predictive contribution survives the predeclared population/relatedness control;
- the direction of contribution is stable under the required public-reference sensitivity analysis;
- the block is not driven by phenotype-linked missingness/callability;
- the same selection procedure is available to both competing models.

Matched genomic-window and within-population phenotype-label permutations are retained as negative controls.

## Model F implementation constraint

Model F may use multiple genomic blocks, but their joint contribution to the focal trait vector must pass through one shared latent genomic score or equivalent rank-1 block-by-trait effect structure.

Trait-specific intercepts/scales are allowed.

Model F therefore tests integration of inheritance, not the unrealistic hypothesis of one causal gene.

## Model R implementation constraint

Model R may estimate component-specific genomic scores using the same training data and complexity budget.

The genomic effect structure may have rank greater than one and blocks may be component-specific.

Residual phenotypic covariance is allowed. Model R does not assume phenotypic independence.

## Primary estimand — untouched joint prediction

For every VALIDATION individual, record without refitting:

d_i = log p_R(y_i | genotype_i, frozen training fit) - log p_F(y_i | genotype_i, frozen training fit)

Primary summary:

Delta_joint = sum_i d_i

Report:
- Delta_joint;
- mean and median individual log-score difference;
- exact/sign-flip paired one-sided test over validation individuals;
- componentwise log-score differences.

Primary R-support criterion:
- Delta_joint > 0;
- paired one-sided sign-flip P < 0.05;
- neither primary component has a negative total held-out log-score difference large enough that the joint result is driven solely by the other component.

No training fit statistic can substitute for this result.

## Secondary estimand — genomic-state separability

After TRAIN candidate blocks are frozen:
- derive the local genomic state/haplotype representation for each primary component;
- quantify cross-component concordance of individual genomic states;
- compare that concordance with matched genomic blocks of similar callability, MAF/heterozygosity and LD/recombination context.

Low concordance alone is not evidence for reconfigurability because recombination makes arbitrary genomic windows differ.

It becomes mechanistic evidence only when the discordant states are the states that independently predict the focal component phenotypes.

## Nature kill-shot estimand — prospective mosaics

Before opening VALIDATION phenotype labels, register individuals whose frozen component-specific genomic scores predict discordant component states.

A prospective mosaic is a validation individual whose colour-predictive genomic state and orientation-predictive genomic state correspond to a combination not forced by the shared Model-F state.

Strong natural-reassembly support requires:
- at least 3 high-confidence prospective mosaic individuals in VALIDATION;
- the componentwise model correctly predicts both primary component states/continuous directions for those mosaics under the frozen scoring rule;
- the prediction is reference-stable;
- the individuals are not explained by taxonomic mixture, post-anthesis reorientation or phenotype-linked missingness.

A single striking recombinant is illustrative, not decisive.

## Decision classes

### FIG3_R_STRONG

All are required:
- primary untouched prediction favors Model R under the frozen criterion;
- both primary components have stable genomic predictive contribution;
- component genomic states are not one shared state;
- at least 3 prospective validation mosaics are correctly predicted;
- inference is reference-stable.

This opens the Nature Fig.4 natural-reassembly claim.

### FIG3_R_PARTIAL

Model R improves held-out prediction but one primary component, the validation route or the mosaic count is insufficient for the full kill shot.

Interpretation: evidence for separable inheritance, but not the full natural-reassembly kill shot.

### FIG3_F_SUPPORTIVE

The shared-state Model F predicts held-out joint phenotypes as well as or better than Model R, and component-specific genomic scores collapse onto the same inherited state.

Interpretation: the focal organ behaves as a persistent integrated evolutionary unit in this system.

### NOT_IDENTIFIABLE

Use if candidate blocks are not reference-stable, genomic signal is too weak, validation n is below 20 for the prospective route, or population/taxonomic/cytotype confounding cannot be separated.

Do not reinterpret NOT_IDENTIFIABLE as support for Model F.

## Claim boundary

Even FIG3_R_STRONG does not establish adaptation. Fig.5 remains a separate functional/fitness test.

The Nature claim becomes:

Phenotypic integration need not be the unit of inheritance; an integrated reproductive organ can be reconstructed from independently assorting genomic components.

It must not be reduced to:
- different traits have different genes;
- GWAS peaks do not overlap;
- local genealogies differ somewhere in the genome.