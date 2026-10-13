# 頭花モジュールの生態学：公開原著で既知の部分と、まだ検証されていない部分を凍結する

2026-10-13 | `aza3` PR #44 | Scientific decision, not a new fitness result. **Cirsium capitulum only** (orientation, actual spines, phyllary geometry, florets, true insect contact).

## First: the strongest existing ecological findings already contradict a simplistic attraction–defence story

**Fenner et al. (2002), `Oecologia`, DOI [10.1007/s004420100773](https://doi.org/10.1007/s004420100773).** Asteraceae survey: 3 populations of each of 20 species. Higher capitulum receptacle width was associated with higher infestation incidence among species; in three species investigated in more detail, size also predicted infestation within species. Capitulum lifespan and flowering synchrony did not explain their broad incidence pattern. These are real within-/among-species *size vs attack* results, NOT proof of adaptive spine/orientation evolution.

**Zhang et al. (2018), `International Journal of Plant Sciences`, DOI [10.1086/698766](https://doi.org/10.1086/698766).** In 34 alpine Asteraceae, larger heads suffered higher infestation incidence, but *conditional on infestation* a lower fraction of achenes were lost. Their species-level product `P(infested) × loss_fraction_given_infested` nevertheless increased with head size. This refutes a single monotone "large heads protect/don't protect seeds" reading. The original paper reports seed size and aboveground biomass as additional axes. We did not obtain/recalculate its 34 raw taxon rows, so DO NOT invent effect sizes or Cirsium subsample ranks.

**Xi et al. (2021), `Journal of Animal Ecology`, DOI [10.1111/1365-2656.13480](https://doi.org/10.1111/1365-2656.13480).** Across 24 alpine Asteraceae, the study linked plant–pollinator and seed-predator networks and found that positive density-dependent pollination/seed-set effects could overwhelm the disadvantage from higher seed predation, leaving rarer plant species with lower overall seed output. Their species-level seed output rate is constructed as `seed_set_rate × (1 - seed_loss_rate)`. **A reversal between mutualistic benefit and antagonistic cost in plant seed output is therefore already published prior art**, and it is not new by itself. Public Dryad DOI [10.5061/dryad.rr4xgxd7p](https://doi.org/10.5061/dryad.rr4xgxd7p) lists `Dataset_network_combination.xlsx` and a README. In this audit the public metadata/API were verified, but attempts to retrieve the original file stream returned access/automation restrictions (403/401). **No species-level numerical rerun has been performed.** The DOI listing and original abstract are evidence of prior work, not permission to guess row-level correlations.

**Wang et al., published online December 19, 2025, issue January 2026, `International Journal of Plant Sciences`, DOI [10.1086/738562](https://doi.org/10.1086/738562).** Among **17 Asteraceae species**, volatile esters/total volatile *concentration*, and nonvolatile phenolic/terpenoid compositional *evenness*, were negatively associated with capitulum seed-predator infestation incidence. Ordinary regression **and phylogenetic GLS** both supported the reported relationships. **Chemical attraction/deterrence is a real head-level competitor to physical orientation/spine effects**, even after phylogeny was considered. Yet this is a between-species observational chemical study, not a causal proof that the physical outer armature is redundant, and not a measurement of *Cirsium* spine touches or post-establishment larval seed consumption. Its supplementary Table A1 is listed publicly; the current retrieval returned 403, so no numeric parameter is attributed to it here.

## The remaining high-value Cirsium question (one, not ten)

> Does a head's true gravity-referenced orientation × genuine spine/phyllary access geometry **selectively change the transition from independently observed insect arrival to confirmed oviposition**, above and beyond head size, floral scent/daypart and reproductive stage, while keeping confirmed effective pollen service, and is that stage-specific function reflected in **absolute mature filled achenes** across *prespecified heads*?

This is **not** generic multi-functionality. The potentially new mechanism is the **measured physical interface a particular guild hits**, including obstruction, failure, rerouting, and final seed outcome, rather than a correlation between visually inferred angle/outline and opportunistically sampled visitors.

### An adversarial comparison that must be passed

We must show that a real `spine_tip` / `phyllary_gap` physical blockage is not just:
- larger heads attracting more egg-laying seed feeders (Fenner 2002; Zhang 2018);
- volatile compounds altering seed-feeder *arrival* and nonvolatile composition altering *performance* (Wang et al. 2025/2026);
- clock, prior pollination (Theis 2005/2007), weather or stage leading insect species into different observed samples;
- insects succeeding by alternate floret portals after a first thorn obstruction;
- parasitoids inducing extended larval feeding in some host systems (Xi et al. 2015, experimental ***Saussurea***, observational *Cirsium setosum*);
- infested-only / emerged-adult-only endpoint selection, which discards true zero/failure heads;
- a mathematically misleading mixture of *probability of infestation*, *loss conditional on infestation*, and *overall plant seed fitness*.

None of the previous 2-D Azami 1,734 image morphometrics is a calibrated botanical `spine_length_mm` or genuine 3-D head orientation. The existing head-only RV patterns are descriptive, and 0/12 predeclared image/climate context contrasts are BH significant. The old 1990 Cirsium parasitoid structural refuge study also demonstrates that "large heads shelter larvae from parasitoids" cannot independently be marketed as a new conceptual discovery.

## Exactly three quantitative quantities that a future Cirsium cohort must report

**1. Head-level infestation incidence** `p = N_infested / N_all_eligible_registered_heads`. "Infested" requires independently confirmed biology; zero must have been actually assessed. Head-level preregistered study enrollment may not depend on adult emergence.

**2. Conditional severity** `s = mean(D_i/N_i | infested)`, where `D_i` are confirmed damaged achenes and `N_i` total achenes in the same head. `p*s` is the **head-weighted** fraction lost across all heads if truly uninfested heads all have `D_i=0`; it is **not** an automatically causal estimate of thorn defence.

**3. Actual outcome** `mean(filled_achenes/head)` and, if measured separately, `germinating_achenes/head`. If asking for a pooled fraction across achenes, use `Σ D_i/Σ N_i` rather than `p*s` because head sizes differ. Preserve pollen limitation, head rank, plant identity, other flowering hosts, and full recording effort. "No visible feeding holes" is not itself independently measured embryo viability.

New executable tool: `analysis/audit_capitulum_incidence_severity_two_denominators_v1.py` directly reads the **pinned native EAzami Aim2 head-field CSV**, validates all registered biological head IDs and `seed_predator_presence`, `seed_predator_damage_fraction`, `total_achenes`, `filled_achenes`, never imputes missing values, and reports head- and achene-weighted decompositions per `population × taxon`. It also notes observed orientation, true measured spine and diameter availability. It makes **no p-value, correlation coefficient, selection coefficient or phylogenetic claim**. A constructed arithmetic example gives the two valid but different denominators `0.30` versus `0.15`; **neither number is a measured thistle seed loss**.

New CI `.github/workflows/audit-capitulum-incidence-severity-v1.yml` pins the untouched EAzami source at `c9870e91604ed1514344ddc8f765eae5237dbf27` and tests missingness, double-counted heads, proportions outside 0–1, impossible achene counts and incidence–conditional severity identities.

## Honest stage-gated scientific conclusion as of October 13, 2026

- **Already supported by prior ecology:** infested incidence and conditional severity can respond oppositely to capitulum size; chemical compound concentration/evenness can predict head infestation; pollination gains can override seed predation costs in multispecies Asteraceae studies.
- **Already supported only descriptively in our morphometrics:** nonuniform among-taxon coupling of head image traits, but not real physical spines, direct contacts, viability, or adaptive evolution.
- **Still unidentified for focal Cirsium:** true orientation×armature access filtering, whether species-specific mutualists and antagonists are physically segregated, whether a head module raises absolute viable seed output, and whether any repeated evolutionary assembly exceeded a phylogenetic/null expectation.

**Practical decision: keep the headline in physical-head interaction ecology and within-population trait→access→fate, not a generic trade-off, parasitoid cascade, or 34-/24-/17-species effect already reported elsewhere.** The native EAzami head and bout ledgers and aza3 contact records continue to have zero real observation rows, so a paper claiming new measured adaptive selection is not yet supported. Field permissions and all frozen Azami/EAzami endpoints remain unchanged.
