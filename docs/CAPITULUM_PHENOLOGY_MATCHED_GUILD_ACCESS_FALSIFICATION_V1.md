# 頭花の「選択的侵入口」は本物か、それとも花期の混合で生じた錯覚か？

2026-10-08 | New source-grounded decision in PR #44. **Cirsium capitula only.** No leaf-trait programme, no field permission change and no reclassification of frozen Azami/EAzami results.

## 結論が変わった根拠：原著のアザミ移植実験

**Vanbergen et al. (2007)**, *Ecological Entomology* 32:419–427, DOI [10.1111/j.1365-2311.2007.00885.x](https://doi.org/10.1111/j.1365-2311.2007.00885.x). A 240-plant *Cirsium palustre* transplant in 24 blocks at two sites directly compared host aggregation, isolation and phenology against *Tephritis conura* and *Pteromalus elevatus*.

Read strictly from original journal abstract:
- Herbivore abundance per plant rose with plant density, but parasitoid abundance and parasitism rate did not similarly respond.
- Patch isolation did not predict insect abundance or parasitism in this study.
- Thistle phenological stage correlated **positively** with herbivore abundance, parasitoid abundance, and percentage parasitism at both patch and plant scales.
- The authors interpreted phenological quality heterogeneity as more important than spatial habitat configuration for this assemblage.
- The study **did not experimentally manipulate real involucral spines, world-referenced head orientation, first contact portals, parasitoid attack timing relative to seed damage, or ultimate viable achenes** for our hypothesis. It cannot establish that physical head gateways changed with phenology.

Source verification: [publisher abstract](https://resjournals.onlinelibrary.wiley.com/doi/10.1111/j.1365-2311.2007.00885.x); [NERC author's repository record](https://nora.nerc.ac.uk/id/eprint/666/), which explicitly says full text not deposited there.

**Russell & Louda (2005)**, DOI [10.1007/s00442-005-0204-3](https://doi.org/10.1007/s00442-005-0204-3): *C. canescens* open head density explains some *Rhinocyllus conicus* egg redistribution on neighboring *C. undulatum*. Same-species traits without synchronized alternative-flower census are insufficient.

**Song et al. (2024)** *Biological Reviews* DOI [10.1111/brv.13060](https://doi.org/10.1111/brv.13060): bract multifunctionality across pollination, defence, abiotic stress and post-pollination photosynthesis, and changes in selection agents across developmental stages, are already recognized. Therefore **'bracts serve different functions at different stages' is not novel**.

The empirical gap specific to Cirsium is narrower: whether independent insect guilds contact **different actual head interfaces** on the *same* developmentally changing organ, and whether observed morphological gates causally affect verified oviposition/pollen transfer and viable achenes.

## The hidden Simpson mechanism in earlier route summaries

The existing `analysis/summarize_capitulum_guild_route_overlap_v1.py` pools independently identified pollinator and seed-feeder candidate approaches across all *Cirsium* phenophases. It computes `Omega = sum_route min(P(route|pollinator),P(route|seed feeder))`.

That marginal overlap is biologically unsafe as a claim about segregated head portals when the insect guilds are active in different phenophases. **A perfect matching of pathways within every phenophase can still yield a low pooled overlap.**

A constructed **synthetic** example:

| Head stage | Pollinator candidate arrivals | Seed-feeder candidate arrivals | Both guilds' route |
|---|---:|---:|---|
| bud | 90 | 10 | above |
| full anthesis | 10 | 90 | below |

Within each stage `Omega=1` (no route segregation at all). Naively pooled `Omega=.2` (looks segregated). The divergence comes entirely from guild composition **by stage**, not from barriers. This is an arithmetic counterexample, not observed Cirsium behavior. The specific ten-versus-ninety values are didactic, not fitted.

## Implementation (append-only; existing script/results unchanged)

- `analysis/audit_capitulum_stage_matched_guild_routes_v1.py`: joins the existing 32-field independently identified approach ledger to the 27-field verified continuous-video-effort ledger using the authoritative head/bout ID; delegates video denominator and no-false-zero checks to the frozen shared validator.
- Produces two distinct outputs: `naive_pooled_overlap` (clearly marked potentially confounded) and **equal-weight, within-population × same-phenophase × same-head-rank** `stage_matched_equal_stratum_overlap`, plus eligible numbers of distinct plants and heads.
- Descriptive minimum per stratum and each guild: 10 independently pre-identified candidate approaches, at least 3 independent plants, and no unclassified world routes. These **logical thresholds are not statistical power calculations**.
- If a stage is represented in only one guild, then **HOLD**, not a route-separation score of zero. A detected route with unknown or post-success guild identity cannot be turned into a negative event. Same-guild repeated encounters on one head remain nested, not independent replicates.
- Same population × stage × head rank **does not guarantee the same calendar date, clock period, weather or local insect pool**. Current ledger has no absolute observation timestamp. Even a matched output is only stage-stratified description, not causal mediation or morphological selection. Prospective field observation must add synchronized date/time plus independently validated insect taxonomy and uniform video sensitivity.
- Zero observation rows must report `NO_REAL_CONTINUOUS_BOUTS`, not no pollination or no herbivory. Synthetic examples never populate real intakes.
- No multiple-testing result and no new p-value is derived; frozen 1,734 Azami images and 0/12 supported context interaction comparisons remain unchanged.

## What would falsify the two ecological models?

**M0 / phenology-first availability:** Insects are present at different stages; once synchronized head stage and independent arrival opportunity are held comparable, actual portal use or reproductive-zone access is similar between contrasting real phyllary/spine states. A visible correlation between insect abundance and head morphology need not be mechanical defence.

**M1 / anatomical gating:** At the **same head stage, same confirmed insect species/host race, and comparable arrival exposure**, authentic phyllary gap, spine direction/rigidity and gravity-referenced head orientation predict *failed access conditional on independent arrival*. This includes genuine failures, repeated attempts and route switching, not only eventual egg-layers. Pollinators may bypass the interface while egg-layers cannot.

**M2 / plant function:** Differences in verified entry/oviposition/pollen deposition then predict mature **viable filled achenes**, distinguishable from rain-exclusion effects, resource availability and alternative host blooms. If anatomy affects adult access but not seed production, mechanism remains an interaction filter, **not measured selection**.

**M3 / evolutionary origin:** Only repeated independent populations/phylogenetic transitions with heritable traits and regime-dependent selection can test whether orientation/spine modules became synchronized or deliberately decoupled by evolution. A positive stage-stratified access observation alone is insufficient.

## Most efficient go/no-go

1. **Before a new factorial study**, review the existing raspberry-Pi continuous recordings (if available and permitted) for the *presence of real independently identified seed-feeding adults and legitimate pollinator candidates*, an unambiguous head stage, and a visible pre-contact approach space. An event-triggered clip is not a valid denominator without recall calibration.
2. **Use 24 heads as feasibility only**, not a powered adaptive interaction experiment. If both relevant guilds never occur at the same phenophase in the same population, report **TEMPORAL_GUILD_SEPARATION_WITHOUT_ANATOMICAL_GATE_TEST** (which could itself be biologically interesting but is not demonstrated mechanical selection).
3. For synchronized same-stage insect encounter evidence, compare actual approach → access failure → oviposition/pollen → mature seed outcomes, then only consider approved manipulations. Species of antagonist, head sex/rank, head phenology, and nearby other *Cirsium* flowers must be kept distinct.
4. Avoid treating source-backed GloBI relationships as in situ co-occurrence: 150,638 raw plant-centric interaction API rows, but the parasitoid second-hop audit yielded **0 independently confirmed same-head host–parasitoid pairs**. Missing API data does not imply the organism is absent.

**Scientific state today:** Direct literature supports phenological context as a strong driver in one Cirsium tritrophic assemblage; existing repository photos show some nonuniform head phenotype covariance. Stage-matched guild-route effects, morphological selection and coadaptation **remain unobserved**.
