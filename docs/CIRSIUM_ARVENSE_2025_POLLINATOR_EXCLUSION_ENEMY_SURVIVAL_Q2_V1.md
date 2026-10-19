# Q2の重要な既存実験：送粉者を排除すると種子食者の羽化数も減る

2026-10-20. **Source review and independent reanalysis readiness**, not a newly established flower-head armature effect.

## 原著に基づく発見（同一Cirsium種の野外操作）

Karpate, Y., Blažek, P., Kollross, J., Rech, A. R., Sam, K. & Jorge, L. R. (2025), *Pollinators, not predators, affect florivore survival in Cirsium arvense (Asteraceae) reveals an experiment based on Structural Causal Modelling*. Research Square **preprint**, DOI [10.21203/rs.3.rs-6634172/v1](https://doi.org/10.21203/rs.3.rs-6634172/v1), publicly posted June 30, 2025. Peer-reviewed journal publication was **not independently established** in this audit.

The work randomized/manipulated access/exposure in naturally occurring female plants of **functionally dioecious** *Cirsium arvense* in **four sites around České Budějovice, Czechia**, beginning June **2022**. Its four treatments were:
1. unmanipulated control;
2. predator exclusion;
3. pollinator exclusion;
4. predator-model treatment.

A generative structural causal simulation was used to assess the proposed experimental analysis and design sensitivity; **those synthetic results are not field measurements**. The paper reports around **three heads per biological mother plant** as the field-scale sampling unit for seed counting, nested in treatment/block/site.

**Direct field observations reported by authors**:
- predator-exclusion/model arms did not exhibit detectable differences from control in pollinator visitation, florivore abundance, the fraction of fertilized seeds, or the fraction of attacked seeds, under this study's detection/design;
- pollinator exclusion **substantially reduced fertilized-seed fraction**, consistent with functional female/male separation and insect pollination dependence;
- pollinator exclusion **reduced the number of florivores emerging from heads**, consistent across four field sites in their reported study;
- paradoxically the **proportion of predated seeds was higher** when pollinators were excluded, even though fewer florivorous adults emerged.

This last result is a genuinely **nontrivial multitrophic sign divergence in the same thistle species**. It is an association within a pollinator-exclusion intervention; it does not by itself demonstrate that pollinators directly facilitate florivore survival through edible fertilized achenes rather than changed plant condition, floral access, life stages or sampling denominators. The paper's phrase "florivores depend on pollinators to complete their life cycle" is the authors' biological interpretation. **Adult emergences, damaged-seed fraction, absolute damaged achenes and absolute viable achenes are not interchangeable response variables**. Sex of mother plant, possible mechanical insect exclusion side effects and four treatment levels must be kept distinct.

## What this means for bita × sch × balance × slk

**Critical correction to Q2 novelty:** The hypothesis "changing pollinator access can change seed-predator abundance and plant fitness in Cirsium" has **already been experimentally examined** (2025 preprint, and earlier work on pollination/florivory). Therefore a future paper cannot claim to discover this generic cross-guild interdependence just by re-presenting it in aza3.

A genuinely unresolved *capitulum architecture* hypothesis remains:

> Within the SAME naturally exposed head and independently identified insect taxa, does an experimentally or quasi-experimentally distinguishable **real phyllary/spine/stickiness/orientation interface** control the pollinator's access and the seed-feeder's access **separately**; and does the resulting change in fertilization affect the *feeding quality or establishment success of enemies* as well as the plant's absolute filled achenes?

A simple monotone "barrier protects plant from seed predators" model becomes too narrow if seed-feeder recruitment or survival **also depends on pollinator-created viable seed supply**. Head phenotype may change enemy pressure by affecting **pollen transfer first**, independent of direct mechanical exclusion. In this mechanism, a trait with positive pollen benefit can increase a potential resource base for enemies. The resource-mediated feedback could oppose or amplify any direct head armature defence, but **neither the feedback magnitude nor trait-dependence has been measured for the aza3 morphological modules**.

## Specific real-data causal and denominator precautions

- The four arms are **not** a fully crossed factorial of independent actual pollinator presence × predator presence × florivore identity × head morphology. Thus they cannot estimate `W00(z), W10(z), W01(z), W11(z)` as a factorial morphological fitness response surface or prove selected trait optima.
- A reduction in total emerged seed-feeder adults and an increase in **fraction** of attacked seeds can coexist without a direct per-insect feeding-rate increase. To distinguish mechanisms, the required absolute head-level counts are `N_florets`, `N_fertilized`, `N_damaged`, `N_unfertilized`, `N_emerged by taxon`, and `N_filled viable achenes`, plus their **exact named denominators**.
- If a pollinator exclusion treatment physically covers the head, it might ALSO affect ovipositor access, weather, microclimate and predator accessibility. The paper's actual treatment details and sham/control definitions must be retained before attributing florivore emergence differences uniquely to seed provisioning.
- The study's experimental origin `2022 / four Czech sites` and female-plant sex cannot be pooled as observations from *C. purpuratum* 1993/94 or *C. discolor* 1980/82.
- Flower-stage-specific access and independent guild IDs must come from observed organisms; broad insect order is insufficient.
- Head armature true spine length, orientation to gravity, and actual contact route were **not documented** in the preprint as randomized traits. Even if this dataset becomes re-analyzable, a physical-spine adaptation coefficient remains **NOT_IDENTIFIABLE**.

## Public data archive: verify before claiming any recomputation

The manuscript explicitly names a publicly deposited data/code archive [Zenodo DOI 10.5281/zenodo.15364152](https://doi.org/10.5281/zenodo.15364152). During this session, regular web direct access to that URL was unavailable, and no raw source file has yet been downloaded or inspected. A new **read-only bounded probe** at `analysis/probe_cirsium_arvense_2025_zenodo_v1.py` uses a GitHub Actions runner to retrieve metadata, verify DOI/record ID, list checksums, and report permitted small openly accessible field-file metadata. Absence/access denial is `ARCHIVE_NOT_FETCHED` / `METADATA_UNAVAILABLE`, not a zero experimental effect.

No raw specimen head file, re-fitted effect, new Cirsium field accession, genome, pollen count or spine-contact outcome is asserted on the basis of this record review.

### Novelty decision and immediate question

Q2 has **an actual head-level multi-guild intervention with a reported sign-divergent outcome in the same Cirsium species**. The key missing causal link is not "do pollinators and seed predators affect one another?" but **which head interface/stage produces that coupling, and whether a physical module shifts it toward a net gain or loss in reproductive output**.

Before a broad new field factorial, test whether the open Zenodo release contains the absolute per-head fertilized/unfertilized/predated-seed counts and emerged insect numbers needed for a **descriptive** denominator audit and field effect replication. If the data are available, preserve original experimental randomization and source licenses. If not, retain only the authors' published outcome directions and note access failure.

## The decisive unresolved mediation experiment

**The overlooked mechanism in the actual paper's METHODS:** The 2022 pollinator-exclusion treatment was not an isolated intervention on successful seed fertilization. The authors covered **at least three heads** with **organza mesh bags** supported by wire frames installed around one week **before** flowering. Thus the intervention might simultaneously alter (a) pollinator transfer and seed provisioning; (b) **adult seed-feeder encounter/oviposition through the mesh**; (c) parasitoid/predator access, head-level abiotic humidity, flowering microclimate and inspection probability. We have not measured seed-feeder adult passage through this particular mesh. **The preprint's reported reduction in emerged adult florivores is an experimental contrast for its BAG TREATMENT, but is not uniquely an effect mediated by the absence of pollination.** A positive pollinator→florivore dependency claim remains plausible yet not isolated.

The source reports **12 blocks of four treatment plants at four sites** (48 blocks, expected 192 assigned plant slots), but separately prints **196 plants**. It also reports **609 mature heads** sampled, counting fertilized, unfertilized, predated seed categories and reared adult guilds. Until actual per-plant CSV is inspected, do not invent a correction for the apparent plant-total discrepancy: it could reflect extra individuals or differing denominators. Preserve `48 blocks`, `196 stated plants` and `609 sampled heads` as **author-reported, not independently resolved**.

### A sharply discriminating biological pilot (prospective, site authorization independent)

There are TWO mechanisms to separate, not merely a single trade-off:

**Resource facilitation:** pollen transfer increases fertilized achene resources, which increases established/enclosed florivore development or successful emergence, conditional on a *known initial number of seed feeder eggs*.

**Physical entry filtration:** mesh or true phyllary/spine/sticky geometry excludes incoming adults/ovipositors, changing the number of eggs laid BEFORE larval feeding, with or without any change in actual achene nutritional quality.

A feasible, **not yet approved** mechanistic experiment requires two temporally separable transitions and genuine insect IDs:

1. Confirm the actual species of seed feeder and **whether/when its oviposition can occur relative to pollination** on the same head. Record external arrival, genuine ovipositor contact, eggs before treatment, complete video effort and head/plant IDs.
2. Where biological timing allows, establish *matched verified initial infestation* in heads, then subject them to the **same mesh/access regime** while randomizing **pollen receipt rescue (verified hand outcross pollination)** versus equivalent sham/unsupplemented floral handling. The exact same structural insect barrier in both arms isolates the additional post-oviposition **resource-supply response**, if pollen truly reaches the stigma and other manipulations are equal. If pollen rescue is impossible in the bag or seed feeders oviposit only after fertilization, this design is **HOLD**, not a failed model with invented zero.
3. In a separate access test, maintain similar flower-head sex, size, reproductive stage and verified pollination histories while varying **permitted real armature geometry vs sham**; measure independently identified adult seed feeders' pre-entry arrival, actual oviposition/blocked attempts, and **absolute** emerged insects/seeds, not only proportions among infested heads. Do not call a bag's effect on insects a Cirsium spine effect.

These designs require independent permits, plant number/power evaluation, full recording, and suitable biological stages. They do not authorize clipping real spines, applying adhesive to natural heads, manipulating pollinator or seed-feeder populations, or collecting offspring without existing approval.

### The falsifiers and maximum claim

- If under **equal verified starting egg loads, equal insect exposure and identical mesh treatment**, additional pollen receipt changes emerged feeders, evidence would favor a resource-mediated link in this species/stage. Still distinguish embryo nutritional quality from head growth/correlated stress.
- If adult approach/egg deposition differs with true head interfaces while pollen is held comparable, that supports direct **access gating**, conditional on arrival. To infer plant fitness still measure real per-head absolute filled viable achenes, not adult eclosion alone.
- If neither link is supported, then a head module can be multifunctional without creating a measurable partner-specific fitness trade-off.
- The **Q2 Nature threshold** is demonstrating the actual *same structural axis* changes relative partner effects on **head/plant fitness**, not repeating the already observed pollinator-exclusion effect on seed-feeder emergence.

**Research programme integration:** This is the most relevant currently identified *within-one-Cirsium-species four-site field intervention* for Q2, stronger than subtracting effect sizes from separate species and experiments. But it alone cannot fulfill Q1 genetic memory or Q3 lost-pigment re-expression, and the true armature structural trait was **not manipulated**.

## Independent same-Cirsium resource manipulation (2008): do not confuse soil fertilization with pollen fertilization

Walker et al. (2008), *Journal of Animal Ecology*, [DOI 10.1111/j.1365-2656.2008.01406.x](https://doi.org/10.1111/j.1365-2656.2008.01406.x), previously used **the same botanical species *Cirsium arvense*** to manipulate *host-plant nutrient supply/quality* and *bud availability*, examining *Xyphosia miliaria* tephritid and Hymenoptera parasitoids in potted/natural systems.

The original author-reported findings distinguish independent routes:
- higher host-plant **nutrient supply and bud availability** increased tephritid abundance;
- the ovipositing tephritid **did not preferentially place eggs on nutrient-fertilized plants**, under that study's direct choice test;
- parasitoid attack was **more frequent on nutrient-fertilized thistles**, and the response to host density was variable;
- bottom-up resource availability and plant quality could matter even when a direct oviposition-preference response was not detected.

**Critical bilingual scientific definition:** in Walker et al.'s original abstract, `fertilized C. arvense` means **mineral/nutrient fertilization of HOST PLANTS**, **NOT pollen fertilization of ovules**. Therefore its lack of host-choice preference on nutrient-fertilized thistles cannot refute the 2025 preprint's suggested consequences of POLLINATOR REMOVAL for developing fertilized-seed resources. The two studies are different field designs, different cohorts, and measure different processes. We did not retrieve Walker's individual plant records or fit any new effect size.

### Guild annotation caution in the 2025 preprint

The 2025 report identifies emerged floral occupants including **Cecidomyiidae, Syrphidae, Platygastridae and Tephritidae**. A family label is **not sufficient** to conclude that every emerging insect is a seed consumer: many hoverfly larvae are predators and platygastrids are parasitoids. Their particular head occupancy and feeding role must be demonstrated from species/morphospecies and host association. Consequently even a downloaded "total emerged insects" variable must not be recoded wholesale as successful **seed-feeder emergence** without an explicit source-backed family/larval-role ledger.

**Biological reconciliation:** Cirsium may change the realized insect community via at least three mechanisms: **(a)** anatomical physical access and ovipositor reach; **(b)** pollen-dependent developing food tissues; and **(c)** nutrient/bud abundance. Their patterns cannot be reduced to a single monotonic "defense strength." To demonstrate a true interaction-partner-dependent payoff for a *head module*, compare real **head structural access** while measuring **absolute species-specific eggs/larvae and plant filled seeds**, not just the treatment differences already published.
