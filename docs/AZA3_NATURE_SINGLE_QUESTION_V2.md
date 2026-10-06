# aza3 Nature single-question pivot v2 — reconfigurable integration

**Status:** active conceptual target; genomic answer not yet known  
**Date:** 2026-10-05

## The one question

> **Must an integrated complex organ evolve as an integrated unit?**

Working answer to test:

> **No. Phenotypic integration may be a repeatedly reconstructed state rather than a persistently inherited evolutionary unit: component couplings can differ across present-day covariance, ecological association, historical change and genomic inheritance, allowing an organ to remain integrated while its parts are reassembled.**

Working label: **reconfigurable integration**.

This label is a hypothesis shorthand, not a claim that the phrase itself is novel.

## What changed

The previous Nature framing assumed a positive alignment:

`ecological modules -> historical modules -> genomic modules`

and treated that alignment as the source of combinatorial evolvability.

The frozen Azami reanalysis does not support the first arrow as a general module-level statement.

Across nine predeclared capitulum constructs and six predeclared environmental blocks:

- phenotypic module cohesion is supported among taxa and within taxa;
- ecological fingerprint cohesion under those same module labels is not supported;
- phenotypic integration strength does not predict environmental-profile similarity.

Therefore the Nature question cannot be:

> do ecological modules match genomic modules?

The stronger question is whether **the same organ can be integrated at one level without preserving the same coupling structure at other levels**.

## Existing evidence before new aza3 population genomics

### Layer 1 — present phenotype is structured and integrated

From the frozen Azami v3 complete-18 common cohort:

- 1,734 observations;
- 42 taxa;
- within-taxon module cohesion: within-minus-between RV = 0.09228, permutation P = 0.0013;
- among-taxon module cohesion: within-minus-between RV = 0.08941, permutation P = 0.0365.

Thus the capitulum is not an arbitrary bag of image traits.

### Layer 2 — ecology does not reproduce the same module partition

The aza3 exact test uses six normalized environmental-block fingerprints for the nine constructs and all 5,040 label assignments preserving module sizes 1/3/2/3.

Among taxa:

- mean within-module cosine similarity = 0.88429;
- mean between-module cosine similarity = 0.87328;
- difference = 0.01101;
- exact one-sided P = 0.36587.

Within taxa:

- difference = 0.01634;
- exact one-sided P = 0.29524.

Separate frozen Azami analysis also finds no coupling between pairwise phenotypic integration and similarity of environmental profiles:

- among taxa: rho = 0.04299, QAP P = 0.8475;
- within taxa: rho = 0.23449, QAP P = 0.1652.

Decision:

> **FIXED_ECOLOGICAL_MODULE_ALIGNMENT_NOT_SUPPORTED**

This does not mean ecology is unimportant.

Two trait-specific ecological anchors survive the complete robustness sequence:

- floral chroma × short-wave radiation;
- presentation angle × annual precipitation.

Those environmental predictors are themselves almost unrelated in the complete-18 observation cohort (Pearson r about -0.003).

Thus the supported statement is:

> **specific capitulum components track distinct ecological dimensions, but those dimensions do not organize the whole capitulum into the same fixed modules seen in phenotype covariance.**

### Layer 3 — historical change is not synchronized

EAzami establishes, within its admitted topology uncertainty:

- orientation: 4–6 minimum changes;
- phyllary posture: 3;
- stickiness: 5;
- unequal relative-depth geometry;
- 0/3 discrete trait pairs pass the robust shared-transition-localization rule.

This rejects the simplest model in which a presently integrated capitulum has one persistent synchronized macroevolutionary history.

It does not establish genetic independence.

### Layer 4 — reference feasibility no longer blocks the genomic test

R1B-0 public-data work shows:

- high focal target-capture recovery;
- strong East-Asian coding-sequence compatibility;
- high clean-locus placement success across three public genome references;
- LOW within-individual hap1/hap2 reference noise.

The preregistered result is:

`HAPLOTYPE_CONTROL_GREEN__PROCEED_TO_OWN_WGS_PILOT`

A focal *C. sieboldii* chromosome reference is therefore not built first.

## Focal-system prerequisite — verify that reassembly substrate exists

The genomic mechanism cannot be tested merely because *C. sieboldii* has several states somewhere across its range.

The public audit currently supports species-wide variation but does not establish same-population independent segregation. Orientation is especially vulnerable to a phenology artefact because ordinary nodding heads become upright after flowering.

Therefore:

`docs/C_SIEBOLDII_COMBINATORIAL_SUBSTRATE_GATE_V1.md`

is required before selecting *C. sieboldii* as the Nature Fig.3–4 biological population.

A `CS_GREEN` result preserves the strongest within-population natural-reassembly route.

`CS_GREEN_LOCAL_SET` retains a weaker local reassembly route.

`CS_AMBER` or `CS_RED_*` means *C. sieboldii* may remain useful for historical/genomic-source questions but should not be forced into the Nature kill shot.

## The decisive aza3 test

The genomic result must distinguish **three**, not two, biological models.

A simple whole-organ-versus-independent-traits comparison is insufficient because Azami already supports a phenotypic module partition. If colour and presentation map to different genomic regions, that alone could simply recover ordinary modular inheritance.

### Model F — whole-organ fixed integration

All admitted focal constructs inherit through one shared genomic state or equivalent whole-organ genetic axis.

### Model M — inherited phenotypic modules

The frozen Azami phenotypic modules are themselves the units of genomic inheritance:

- presentation;
- colour;
- head form;
- involucre armature.

Different modules may have different genomic histories, but constructs inside each predeclared module must share the same module-level genomic state.

This is the strongest conventional modularity baseline.

### Model R — reconfigurable integration

Present phenotypic modules do not impose the genomic partition.

- construct-specific genomic states may cut across the frozen module boundaries;
- component-associated genomic states can assort or recombine independently;
- residual phenotypic covariance and present-day module cohesion can remain strong;
- new genomic mosaics can therefore reconstruct familiar or novel integrated head configurations.

The Nature-scale comparison is **R versus M**, not merely R versus F.

The preregistered test is in:

docs/FIG3_RECONFIGURABLE_INTEGRATION_TEST_V1.md

The strongest route requires at least two admitted multi-trait phenotypic modules, at least five admitted constructs, a training set of at least 40 usable focal individuals, and at least 20 untouched validation individuals. The legacy 40-individual C. sieboldii allocation is therefore discovery-only unless expanded after Gate CS and the phenotype-blind WGS pilot.

Colour × anthesis orientation remains the natural-mosaic demonstration, but because they belong to different frozen phenotypic modules, that pair alone cannot establish FIG3_R_STRONG.

## The Nature kill shot

The strongest single result would be:

> **the frozen phenotypic-module partition fails to describe genomic inheritance, and prospectively identified natural genomic mosaics predict new component combinations in untouched individuals.**

Even stronger:

> a genomic mosaic not used to define the model prospectively predicts a capitulum combination not represented in the training populations.

That turns the paper from retrospective history into a prediction about what combinations evolution can generate.

Fig.3–4 do **not** claim that those alternative combinations are equally fit or equally functional. Functional coherence after reassembly is tested only in Fig.5.

## Five-figure paper contract

The paper must not read as five parallel results. Fig.1–2 are the setup; Fig.3–4 answer the one question; Fig.5 establishes biological consequence.

### Fig. 1 — The organ is genuinely integrated

Show the strongest present-day phenotypic integration result using the frozen Azami module partition.

Claim:
**there is a real integrated phenotype whose evolutionary unit can be tested.**

Do not overload Fig.1 with all ecological associations.

### Fig. 2 — Present integration does not imply one persistent cross-level partition

Compress the two pre-genomic warnings into one figure:

- ecological fingerprints do not reproduce the frozen phenotypic-module partition under the exact null;
- historical component changes are unequal-depth and nonsynchronized;
- retain the two robust component-specific ecological anchors as examples, not as a new ecological-module theory.

Claim:
**the observed phenotype is integrated, but existing ecology/history already make a persistent inherited partition nontrivial rather than assumed.**

Fig.2 motivates the genomic model comparison; it does not itself prove reconfigurability.

### Fig. 3 — What is the inherited unit?

This is the primary new aza3 result.

On untouched validation individuals compare:

- F: one whole-organ genomic state;
- M: one genomic state per frozen Azami phenotypic module;
- R: a complexity-matched genomic architecture allowed to break the frozen module boundaries.

Primary claim opens only if R predicts untouched multi-construct phenotypes better than M under the frozen E1 contract.

This is the conceptual center of the paper.

### Fig. 4 — Prospective natural reassembly

Before validation phenotypes are opened, register genomic mosaics predicted from the TRAIN-fitted component states.

Show:

- the predicted component genomic states;
- the natural individuals carrying those mosaics;
- observed colour and A1–A2 orientation after unblinding;
- reference/taxonomic QC.

Strong route requires at least three correctly predicted prospective mosaics.

Claim:
**the genomic architecture predicts natural component combinations that were not used to define the model.**

This converts reconfigurability from a retrospective description into a falsifiable prediction.

### Fig. 5 — Reassembly still yields a coherent reproductive phenotype

Use one preregistered functional route, probably orientation:

orientation -> wetting/presentation mediator -> reproductive process -> viable achenes

Fig.5 is the first point at which functional or fitness coherence is claimed.

**Fig.5 is an independent extension, not an AND prerequisite for the central Nature evolutionary-unit claim.** The core claim is closed by FIG3_R_STRONG plus the preregistered Fig.4 prospective-mosaic criterion.

It does not need to prove every component adaptive.

Claim, if successful:
**breaking the inherited module partition need not destroy organ-level reproductive function.**

If Fig.5 fails or remains incomplete, the genomic reconfigurability result is not reclassified; only this functional extension remains closed.

## What is genuinely new and what is not

Not new on its own:

- organisms are modular;
- different definitions of modularity need not correspond;
- phenotypic integration can constrain or facilitate evolution;
- mosaic evolution occurs;
- standing variants can be recombined;
- developmental modules can change across lineages.

These are established literatures.

The specific empirical gap targeted here is narrower and more testable:

> **take a phenotypic-module partition defined before genomics, encode it as an explicit inheritance model, and ask whether a complexity-matched genomic model that is allowed to break those module boundaries predicts untouched natural phenotypes and prospectively registered mosaics better.**

That prediction test — not the observation that different modularity definitions can disagree — is the intended novelty.

The Nature-scale advance would be to close, in one natural radiation and one complex reproductive organ:

`present integration`
→ `cross-level decoupling`
→ `separable genomic inheritance`
→ `natural reassembly`
→ `functional integrated outcome`.

The key conceptual move is not “modules exist.”

It is:

> **integration itself need not be the unit of inheritance.**

## Why this could change the usual view

A common operational shortcut in studies of complex traits is to infer evolutionary constraint or evolvability from a detected covariance/module partition.

The literature already warns that developmental, functional, variational and evolutionary modularity need not coincide.

This study becomes high-impact only if it demonstrates the stronger consequence:

> **a statistically integrated organ can retain high combinatorial evolvability precisely because the couplings that integrate its present phenotype are not permanently locked into its historical and genomic inheritance.**

That changes the interpretation of phenotypic integration from a persistent architecture to a potentially transient outcome.

## Critical literature boundary

The conceptual background that must be acknowledged explicitly includes:

- Armbruster et al. 2014, *Philosophical Transactions B*: integration/modularity have multiple meanings and their relation to evolvability is not one-to-one.
- Cheverud & Marroig 2016, *Annual Review of Ecology, Evolution, and Systematics*: correspondence among variational and developmental modularity is not guaranteed.
- Dellinger et al. 2019, *Communications Biology*: floral functional modularity can change with pollination regime and affect evolutionary rate.
- Parins-Fukuchi 2020, *Evolution*: shifts in modularity can accompany bursts of mosaic evolution.
- Smith et al. 2020, *Frontiers in Ecology and Evolution*: evolutionary modularity can facilitate exchange of phenotype components without requiring strict developmental modularity.
- Houle & Rossoni 2022, *Annual Review of Ecology, Evolution, and Systematics*: functional/developmental modularity is often assumed to map to variational/evolutionary modularity without a demonstrated reason.
- Evans & Felice 2026, *Nature Reviews Biodiversity*: integration/modularity can shape diversification and module breakup/reassembly remains a central open problem.

Therefore the paper must not claim that cross-level mismatch itself is unprecedented.

It must demonstrate **the evolutionary consequence of that mismatch**.

## Falsifier

The reconfigurable-integration hypothesis fails as a Nature-scale mechanism if dense same-individual genomics shows that:

- the whole-organ Model F predicts untouched phenotypes as well as or better than more flexible alternatives; **or**
- the frozen phenotypic-module Model M predicts untouched phenotypes as well as Model R, so present phenotypic modules are approximately the units of inheritance; **or**
- prospective genomic mosaics do not independently predict component states.

A result in which M beats F but R does not beat M is scientifically informative modular inheritance, but it is not the proposed Nature-scale reconfigurable-integration result.

## Immediate execution order

1. Keep the Fig. 1 ecological-module null; do not search for another grouping that makes it significant.
2. Run Gate CS: verify whether *C. sieboldii* actually contains a same-population or tightly local multi-component substrate at standardized anthesis.
3. If Gate CS is GREEN, select those verified individuals/populations for the 8-individual WGS reference-transferability pilot.
4. If WGS_TRANSFER_GREEN, freeze the whole-organ vs phenotypic-module vs reconfigurable estimands before scaling.
5. If Gate CS is AMBER/RED, verify an alternative focal system rather than forcing *C. sieboldii* into the reassembly claim.
6. If Gate CS is GREEN and the focal system is retained, expand the focal mapping design only enough to preserve at least 40 TRAIN and 20 untouched VALIDATION individuals after QC; do not assume the legacy n=40 is Nature-ready.
7. Test whether Model R predicts untouched multi-construct phenotypes better than the frozen phenotypic-module Model M.
8. Open the natural-reassembly/prediction test only after phenotype-predictive genomic separation is established.
9. Keep Gate 0 acceleration as context, not a prerequisite.
