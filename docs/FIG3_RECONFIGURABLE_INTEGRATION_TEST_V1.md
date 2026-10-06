# aza3 Fig.3–4 contract v2 — whole-organ, module-level or reconfigurable inheritance?

Status: prospective contract; unopened until Gate CS and own-WGS reference transferability pass
Question: Must an integrated complex organ evolve as an integrated unit?

## Why three models are required

A two-model comparison between one whole-head genomic state and separate colour/orientation states is too easy to win, because colour and presentation are already assigned to different phenotypic modules in the frozen Azami phenotype architecture.

The decisive test must therefore distinguish:

F — whole-organ inherited integration;
M — conventional phenotypic-module inheritance;
R — reconfigurable integration in which even the frozen phenotypic-module partition is not the obligatory unit of genomic inheritance.

Nature-scale support requires R to beat M, not merely F.

## Frozen phenotypic module partition

Use the Azami predeclared construct partition without regrouping after genomic results are seen:

- presentation: presentation_angle;
- colour: floral_lightness, floral_chroma, floral_hue;
- head_form: head_elongation, head_compactness;
- involucre_armature: involucre_form, projection_prominence, projection_pattern.

All same-individual focal plants should be photographed/scored so these constructs can be recovered where technically possible.

### Focal construct-admission rule

Before phenotype-linked genomic discovery, freeze measurement reliability and variation QC for each construct.

The strongest module-vs-reconfigurable test requires:
- at least two predeclared multi-trait modules;
- each admitted module contributes at least two constructs;
- at least five admitted constructs total.

If this condition fails, pairwise natural-reassembly analysis may proceed, but FIG3_R_STRONG is closed.

## Model F — whole-organ integration

All admitted focal constructs inherit through one shared genomic state or equivalent rank-1 block-by-trait architecture.

Multiple genomic blocks are allowed, but their trait effects must pass through the same whole-organ genomic axis.

## Model M — inherited phenotypic modules

Each frozen Azami phenotypic module may have its own genomic state, but constructs inside a module are constrained to share that module-level state.

This is the strongest conventional modularity baseline.

Model M is allowed the same number of admitted phenotypic modules as frozen before genomic analysis; module membership cannot be changed after seeing genotype results.

## Model R — reconfigurable integration

Admitted constructs may have component-specific genomic states and effect directions that do not have to respect the frozen Azami module partition.

Residual phenotypic covariance and present-day module cohesion are allowed.

Model R therefore asks whether present phenotypic modules are outcomes of development/function rather than compulsory units of inheritance.

## Gate prerequisites

1. Gate CS must be CS_GREEN for the strongest within-population route, or CS_GREEN_LOCAL_SET for a claim-limited local-set route.
2. The own-WGS reference-transferability pilot must be OWN_WGS_GREEN.
3. Same-individual taxonomic-confidence, phenotype, image, 2C/cytotype-QC and genomic records are required.
4. The admitted construct set and frozen module labels are locked before genomic discovery.

## Natural-mosaic pair

Gate CS primary colour × A1–A2 anthesis orientation remains the preregistered natural-reassembly pair.

However, colour and presentation belong to different frozen phenotypic modules.

Therefore colour × orientation mosaicism is necessary for the clean natural-reassembly demonstration but **cannot by itself establish FIG3_R_STRONG**.

The R-vs-M module test is the primary conceptual discriminator.

## Sample-size boundary

The 8-individual own-WGS pilot is technical only and cannot discover phenotype-linked genomic regions.

The legacy 40-individual C. sieboldii allocation is discovery-only unless expanded to leave an untouched validation set.

Strong route:
- TRAIN target >=40 high-confidence focal individuals;
- untouched VALIDATION >=20;
- strongest validation is an independent segregating population;
- spatially blocked within-population validation is allowed only with a lower transportability claim.

Minimum strong-route total after QC: 60.

**This is a design floor, not a power-derived target.**

After the phenotype-blind 8-individual WGS pilot passes reference transferability, final n is chosen prospectively from genotype-only quantities: callability, MAF/heterozygosity, LD/effective block count and relatedness/population structure, combined with a frozen grid of biologically meaningful F/M/R effect sizes.

Validation phenotypes and phenotype-linked genomic discoveries cannot enter that power calculation.

If a feasible field/sequencing n cannot discriminate M from R with the frozen target power under that simulation grid, the Nature predictive route is NOT_IDENTIFIABLE rather than an underpowered search for a positive result.

If validation n <20, prospective Fig.4 prediction is closed.

Final n may be increased after the phenotype-blind WGS pilot using genotype-only callability, LD and allele-frequency information. No phenotype-linked genomic result may justify the expansion.

## Frozen split and common covariates

TRAIN/VALIDATION individual IDs are frozen before phenotype-linked genomic discovery.

Both models receive the same global ancestry, population/spatial structure, 2C/cytotype QC, callability and pre-frozen batch covariates.

All genomic block selection occurs in TRAIN only.

## Model-comparison fairness

F, M and R use the **same TRAIN-derived candidate genomic-block pool**, genotype representation, ancestry/callability covariates, preprocessing and missing-data policy.

All tuning is nested inside TRAIN. VALIDATION phenotypes cannot select blocks, penalties, ranks or hyperparameters.

Most importantly, M and R must be matched for effective model complexity. R cannot receive a larger free parameter budget simply because it relaxes module boundaries.

The intended contrast is the structure of the same block-by-construct effect matrix:

- F: rank-1 whole-organ structure;
- M: effects constrained to the frozen Azami phenotypic-module partition;
- R: effects may cross or break those frozen module boundaries under the same effective parameter/regularization budget.

A module-label permutation control preserves the number and sizes of modules while randomizing which constructs are grouped. This tests whether the biological Azami module partition is more informative than an arbitrary grouping of equal complexity.

If a fair complexity match cannot be implemented, E1 is NOT_IDENTIFIABLE rather than evidence for R.

## Primary estimand E1 — R versus phenotypic-module inheritance

### E1 scoring scope

The core R-vs-M joint score uses only constructs belonging to prospectively admitted multi-trait modules.

A one-trait module cannot discriminate M from R because a module-level state and a construct-level state are the same model for that module. Therefore the single-trait presentation module is excluded from the E1 core model-advantage score.

This is why the strongest route requires at least two admitted multi-trait modules and at least five admitted constructs.

The colour × A1–A2 orientation pair remains the separately preregistered E3 natural-mosaic demonstration. Orientation can therefore provide the natural-reassembly kill shot without artificially inflating the R-vs-M model comparison.

For every untouched VALIDATION individual:

d_RM_i = log p_R(y_i | genotype_i, frozen training fit) - log p_M(y_i | genotype_i, frozen training fit)

Primary summary:

Delta_RM = sum_i d_RM_i

Also record:

d_MF_i = log p_M(y_i | genotype_i, frozen training fit) - log p_F(y_i | genotype_i, frozen training fit)

This separates ordinary modular inheritance from true reconfigurability.

Primary R support requires:
- Delta_RM > 0;
- one-sided paired sign-flip P_RM < 0.05;
- no admitted multi-trait module has a strongly negative held-out R-vs-M contribution that reveals the joint advantage is driven entirely by one module.

Training fit cannot substitute for the untouched result.

## Secondary estimand E2 — genomic-state separability

After TRAIN candidate blocks are frozen, compare local genomic-state/haplotype grouping among constructs.

Discordance among arbitrary genomic windows is expected because of recombination and is not evidence for R.

Mechanistic support requires discordant genomic states to be the same states that independently predict the admitted focal constructs, with matched-window controls for callability, allele-frequency/heterozygosity and LD/recombination context.

## Fig.4 kill-shot E3 — prospective natural mosaics

Before opening VALIDATION phenotype labels, register validation individuals whose frozen colour and orientation genomic scores predict a discordant combination.

Strong natural-reassembly support requires at least 3 high-confidence prospective mosaics whose colour and A1–A2 orientation are both correctly predicted under the frozen componentwise model, with reference stability and taxonomic QC.

A single striking recombinant is illustrative, not decisive.

## Decision classes

### FIG3_R_STRONG

Required:
- the construct set satisfies the full multi-module admission rule;
- untouched validation favors R over M under E1;
- at least two admitted multi-trait modules contribute stable genomic information;
- genomic-state separability is phenotype-predictive rather than arbitrary window discordance;
- at least 3 prospective colour × orientation validation mosaics are correctly predicted;
- inference is reference-stable.

This opens the strongest Nature Fig.4 reconfigurable-integration claim.

### FIG3_R_PARTIAL

R improves held-out prediction over M, but module coverage, validation transportability or mosaic count is insufficient for the full kill shot.

### FIG3_MODULE_INHERITANCE_SUPPORTIVE

M predicts held-out phenotypes as well as or better than R and better than F.

Interpretation: the observed phenotypic modules are approximately the units of genomic inheritance.

### FIG3_WHOLE_ORGAN_SUPPORTIVE

F predicts held-out phenotypes as well as or better than both M and R.

Interpretation: the organ behaves as a persistent integrated inherited unit.

### NOT_IDENTIFIABLE

Use when genomic signal, reference stability, module coverage, validation n or biological confounding prevents discrimination.

Do not reinterpret NOT_IDENTIFIABLE as support for F or M.

## Claim boundary

Different traits having different genes is not enough.

Non-overlapping GWAS peaks are not enough.

Local genealogies differing somewhere in the genome are not enough.

The Nature-scale result is specifically that a model allowed to break the predeclared phenotypic-module inheritance structure predicts untouched phenotypes and natural mosaics better than the module-constrained alternative.

Even FIG3_R_STRONG does not establish adaptation; Fig.5 remains a separate functional/fitness test.