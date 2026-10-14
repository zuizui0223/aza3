# Cirsium canescens across 1990–2025: when does protection of seed heads actually matter?

**2026-10-14, PR #44. Comparative primary-literature ecological result.** No new raw field data, no new mathematical re-fit, and no phylogenetic/trait-fitness adaptation claim.

## One Cirsium species, three different empirical intervention designs

The direct focal species *Cirsium canescens* has an unusually valuable, pre-existing ecological evidence chain. It makes a pure "more thistle seeds therefore more thistle adults" headline empirically untenable as a universal law, but also prevents describing *all* recruitment as microsite-limited.

| Evidence | Location and ecological unit | Intervention and follow-up | What the original authors actually found | Boundary |
|---|---|---|---|---|
| **Louda, Potvin & Collinge 1990** — *American Midland Naturalist* 124:105–113, DOI [10.2307/2426083](https://doi.org/10.2307/2426083) | Nebraska Sandhills prairie, including open disturbed versus grassland exclosure context | Direct contrasts of **predispersal insect seed predation**, vertebrate postdispersal seed loss, and competition with established plants | About **3-fold fewer viable seeds** under insect predation associated with **6-fold fewer seedlings** and **6–37-fold fewer new adults**. Germination 1.7%; grass competition reduced seedling survival by 9.5-fold. | Original 1980s experiment and life-stage-specific comparisons; not inferred from the 2025 Colorado head CSV, and no spine/orientation variables |
| **Rand, West, Russell & Louda 2020** — *Oecologia* 193:143–153, DOI [10.1007/s00442-020-04656-2](https://doi.org/10.1007/s00442-020-04656-2) | **Arapaho Prairie Preserve and Niobrara Valley Preserve**, two long-term Nebraska Sandhills sites | **Four levels of seed addition**, seedlings **and reproductive adults** monitored for **9 years** | Both stages were strongly seed-limited. The strength of seed limitation **decreased with plant age**. Postdispersal loss mattered, with stronger density dependence at one site than the other. | Intervention was seed addition over a range of doses, not insecticide on a mother plant; no direct floral head geometry or pollination assay |
| **Russell, Taylor & Louda 2025** — *Ecosphere* 16:e70310, DOI [10.1002/ecs2.70310](https://doi.org/10.1002/ecs2.70310) | **Five Colorado sites in upper Arkansas River Valley** at the xeric western species range margin | Plant-randomized insect-herbivore reduction in 2017–2019, with maternal seeds plus seedlings within **2 m** monitored in the **following two growing seasons** | Experimental head-herbivore reduction increased plant-level viable-seed means strongly and variably by year: control vs insecticide **440.9/770.6 (2017), 361.1/960.5 (2018), 345.3/1103.9 (2019)**. However the extra seed output did **not produce a statistically detectable treatment increase in nearby seedling recruitment**. | Lack of detected treatment effect is not proof that *no* seed limitation exists over all possible seed-addition doses, longer times or microhabitats; 2018 Tanglefoot artifact corrected by authors |

**Important:** The 1990 three/six/37-fold relationships are published author results for that earlier experimental system; the 2025 mean seed-output numbers are published Table 4 author aggregates. None is calculated by our current source-code reader. The 2020 seed-addition result is also directly published. This file reports **source-backed comparisons**, NOT a statistical test of differences between the original experiments.

### What is *and is not* a contradiction

These are the **same plant species but NOT the same population/site, experiment, treatment dose, spatial sampling radius, life-history observation horizon, or response variable**.

- **1990** establishes that head-level seed predation *can* be strongly amplified in later life-history stages; vertebrate seed losses and grass competition were both important.
- **2020** establishes that increasing seed input across multiple levels *can* increase reproductive adult densities at two Nebraska Sandhills sites even when postdispersal mortality also limits establishment.
- **2025** establishes that reducing head herbivores and thereby increasing estimated maternal seed output **did not yield proportionately more locally recruited seedlings** in five Colorado montane range-edge populations, given the induced seed increment and two-year, 2-m-radius outcome.

Therefore, the real inference is **context-, dose-, scale- and life-stage-dependent transfer from head seed production to observed recruitment**. It is *not* verified evolutionary loss of seed selection at the range edge; not a demonstration of site-specific head morphology; not an effect of true spine geometry. The 2025 authors discuss **limited safe microsites**, but postdispersal mechanism was not directly manipulated as part of the head-insecticide comparison, so the causal contributions of safe-site availability, ant granivory and soil moisture need separate evidence.

## The missing transition that could make whole-head defense ecologically important

What does an additional viable seed saved by a head-level physical barrier buy its bearer, across contrasting contexts?

The data contain **two separable measurable quantities**, requiring different randomized observations:
1. `head_to_maternal_seeds`: effect of real head armature/orientation upon **confirmed seed-feeder access, pollen receipt, absolute filled/viable maternal achenes**. A Cirsium canescens insecticide experiment **supports the importance of herbivory for maternal viable seed production**, but does not identify *physical spine traits*.
2. `maternal_seeds_to_recruitment`: marginal extra seedlings and reproductive descendants per actual maternal seed increment under **replicated seed addition or measured microsite exposure**, keeping ambient soil/moisture/predator context. This transfer varied in the historical experiments. The 2025 experiment cannot in isolation quantify a transferable marginal response function at all seed doses.

**If** selective true armature/access filtering produces more viable achenes but the added seeds saturate microsite opportunities, that trait may yield a direct **maternal viable-seed output** advantage but little **local next-generation density increase** in that context. Selection on the trait is not automatically absent: progeny can disperse beyond the 2-m sample window, survive through seed banks, contribute to male reproduction in descendants, or establish in rare microsites. Plant individual lifetime fitness and local population recruitment are not interchangeable.

Conversely, when seed limitation is strong, avoiding seed-feeder oviposition in capitula **may** have substantial demographic consequences. This prospective conditional hypothesis is relevant to modular head evolution, but the idea that seed versus safe-site limitation can modify herbivore impact is classic prior ecology (1990–2025), **not itself a novel thesis claim**. Our genuinely missing biological step remains **measured orientation × real phyllary/spine architecture as a guild-specific physical gate at the same head**, which these three studies did not test.

## Prospective falsifiable comparisons that do not demand a premature giant experiment

- **H1: barrier-access mechanism, before seed outcome:** same species/host race of insect at same head stage and clock interval must demonstrably hit the *real* spine/phyllary barrier, undergo documented failed and successful attempts, with bypass/alternative portal tracked. If true contact doesn't happen, reject that particular physical-blocking mechanism for the observed setting.
- **H2: maternal-fitness route:** head cohorts enrolled *before insect arrival* retain complete damage/incidence denominators and their mature filled/viable seed outcomes, with species identity, head rank, flower sex, size, real orientation, head geometry, pollen limitation, prior pollination/scent context, and environmental context recorded. Compare absolute achenes, not only percentages among invaded heads.
- **H3: recruited-offspring extension (optional independent study):** compare extra viable seeds under independently varied safe-site availability or seed addition with tracked recruits. A two-year 2-m maternal neighborhood cannot substitute for a 9-year adult-density experiment. Avoid claiming generic evolutionary selection differences from published group means.

Before meeting these requirements, keep the highest original fitness claim as **NOT IDENTIFIABLE**. The existing local `EAzami` actual head, video and authentic contact tables still contain **zero biological observations**, and the external 2025 Dryad *capitulum* raw CSV remains inaccessible via anonymous endpoint (401/403), though its metadata, variable contract, author Table 4 means and checksum were verified.

## Citations and quality gates

- Louda, Potvin & Collinge (1990): [DOI 10.2307/2426083](https://doi.org/10.2307/2426083); [University of Nebraska archive](https://digitalcommons.unl.edu/bioscilouda/3/). The published 3-fold/6-fold/6–37-fold figure is a context-specific original result.
- Rand et al. (2020): [DOI 10.1007/s00442-020-04656-2](https://doi.org/10.1007/s00442-020-04656-2), PMID 32322985; [USDA ARS publication record](https://www.ars.usda.gov/research/publications/publication/?seqNo115=368722).
- Russell et al. (2025): [DOI 10.1002/ecs2.70310](https://doi.org/10.1002/ecs2.70310); source [Dryad v3 DOI 10.5061/dryad.kh18932dc](https://doi.org/10.5061/dryad.kh18932dc).
- **No pooled across-study effect sizes, standardized slopes, posterior trait-selection gradients, species-level raw reanalysis or Nature-level novel evolutionary result is claimed.** These studies apply to distinct sites and interventions and cannot be pooled as biological replicates.
