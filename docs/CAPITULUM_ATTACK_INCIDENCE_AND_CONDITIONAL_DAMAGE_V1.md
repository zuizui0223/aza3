# 頭花が「守る」とは何か：食害成立率と成立後被害量を分離する

2026-10-13 | Source-backed ecological decision and reproducible data-readiness audit in `aza3` PR #44. **NOT a new fitted Cirsium selection effect.**

## Newly verified primary-study constraint: opposing trait signs are already known in Asteraceae

Zhang, Yang, Xi, Niklas & Sun (2018), *International Journal of Plant Sciences* 179:554–559, DOI [10.1086/698766](https://doi.org/10.1086/698766), conducted a **34-species Asteraceae** alpine-meadow survey. They found:
- increasing capitulum size was associated with **higher probability of a head being infested** by predispersal seed predators;
- among infested heads, the **fraction of achenes lost was lower** in species with larger heads;
- nevertheless, their species-level overall predispersal achene loss rate (infestation rate × loss rate conditional on infestation) increased with larger heads.

**This is real observed ecological sign discordance.** However, it is a comparative observational result across Asteraceae, **not a fitted effect of true thorn length, structural orientation, mechanical obstruction, or genotype-specific natural selection in Cirsium**. The species list was not taxonomically re-analysed in this PR; do not assume the 34 spp. include a focal Cirsium population. Confounding among head size, phenology, seed number/size and phylogenetic relatedness remains possible, and the paper's abstract does not by itself identify why conditional severity decreases.

Independent Cirsium example: Gijsman et al. (2020), *Global Ecology and Conservation*, DOI [10.1016/j.gecco.2020.e00945](https://doi.org/10.1016/j.gecco.2020.e00945), documents seed losses and head-rank dependence under *Larinus planus* infestation of *Cirsium pitcheri*. Weevils leave oviposition scars on phyllaries, and experimentally documented site and head-rank differences complicate morphology-only interpretation. **This is not experimental proof of a beneficial thorn.**

## The core distinction for the head-access research

The meaningful analysis splits **at least four separate biological events**:

1. **Head arrival and choice**: insect encounters a preregistered flower head, including head-specific verified no-insect observation intervals; scent, display, timing and accessibility can influence this stage.
2. **Egg placement / initial establishment**: independently identified seed-feeder candidate successfully places eggs, measured per eligible approach **including failed approaches**; oviposition need not entail adult body entry if ovipositors penetrate phyllaries.
3. **Predation conditional on infestation**: larval survival, seed consumption, parasitoid-induced delayed feeding, and amount lost among actually infested heads.
4. **Total head and plant outcome**: absolute damaged and mature filled achenes, head rank, sex, potential ovule/floret number and pollen limitation; do not select only heads with emerged insects.

### Descriptive decomposition, under exactly stated denominator

For **head-weighted** proportional seed-feeder loss with no focal insect damage in uninfested heads:

`mean_head_fraction_damaged = P(head_infested) * E(fraction_damaged | infested)`.

If asking for a **pooled proportion across ALL achenes**, however, heads are weighted by their number of achenes, and the calculation is instead:

`sum(damaged achenes across all registered heads) / sum(total achenes across those heads)`.

Those two estimands are **NOT interchangeable** when large heads contain more achenes; this distinction is particularly important when head size is a functional module. Furthermore, actual filled-seed fitness must separately include pollen failure, abortion, physical predation, and fruit maturation. Absence of oviposition scars alone may miss small scars or cryptic feeders.

### What this implies for thorn × direction × phyllary anatomy

Do not say "large or more exposed heads are good for enemies" based solely on increased infestation probability, nor "spiny or larger heads defend seeds" based only on lower losses among **already infested** heads. Either can be true while their weighted consequence goes the opposite direction.

Potential distinct processes:
- TRUE spine/gap/toughness **limits physical first access**, among the candidate insect approaches that actually hit the tissue;
- a large visible capitulum or scent plume **increases encounter** regardless of armature, confounding a selected-on-approach analysis;
- even after establishment, head tissue volume, parasitoid accessibility, consumer crowding and maturation could change **per-larva absolute damage**; none is assumed from body-size photographs;
- world-referenced head orientation can rotate the same morphological head architecture relative to insect paths/rain/UV, but cannot itself change floral-head internal anatomy.

The key hypothesis is **which functional stage is altered** by each head module. Distinct trait combinations could increase the **arrival of legitimate pollinators** while reducing **successful egg establishment** or conditionally reducing larval damage; some configurations might instead increase both. This must be measured, not inferred from RV covariance or theory alone.

## Hard empirically auditable GO / HOLD requirements

The new `analysis/audit_capitulum_evidence_chain_readiness_v1.py` joins native EAzami head and observation-bout tables to exact aza3 video-effort bouts, independently identified approach episodes, validated parent attempted contacts, authentic anatomical contact annotations, and final filled-achene measurements. It rejects invented head/bout/episode keys, source-stage mismatches, duplicated head enrollments, and missing original head registration. It uses existing strict continuous-video and anatomical-contact validators, and reports stage-by-stage availability **without fitting an ecological effect**. The native EAzami source schemas are fetched at immutable `c9870e91604ed1514344ddc8f765eae5237dbf27` in CI.

As of 2026-10-13 the source EAzami head registry, its actual observation-bout ledger, and all associated **current** aza3 true armature/video/approach intake tables contain **zero organismal rows**. Consequently:
- `NO_REAL_FIELD_COHORT` and `HOLD_NO_JOINT_HEAD_LEVEL_ENDPOINTS` are the correct **data availability** answers;
- `0 insects`, `0 thorn effects`, `0 seeds`, `fitness benefit`, and `adaptive coevolution` are **scientifically unsupported conclusions**;
- even if actual contact attempts are documented, **head fitness cannot be inferred until preselected heads and their filled-achene outcomes are linked**.

**Current priority:** distinguish arrival, egg establishment, post-establishment severity and final achene-weighted loss on real replicated Cirsium heads *before* claiming spine–orientation evolution. The 2018 multispecies paper proves why head-defence process rates can change in opposite directions; it does **not** demonstrate our focal causal mechanism or authorize new field manipulations.
