# A capitulum is an interaction venue across time — bounded evidence and testable evolutionary hypothesis

Date: 2026-10-08. Status: **source-backed ecological inference + prospective test, NOT confirmed spine/orientation adaptation**. A live full-genus GloBI query has been submitted but not completed; counts here come from frozen **curated literature anchors** and audited Tofts species lists, not the global GloBI query.

## The ecological insight that changes the hypothesis

The current two-function slogan `flowers attract pollinators but spines protect against herbivores` is insufficient because the *same head* is used by different organisms at different times, and because a parasitoid's attack may happen after the initial seed-feeder damage.

Verified system 1: Walker, Hartley & Jones (2008), *Cirsium arvense*, doi:10.1111/j.1365-2656.2008.01406.x (https://besjournals.onlinelibrary.wiley.com/doi/10.1111/j.1365-2656.2008.01406.x).

- *Xyphosia miliaria*, *Terellia ruficauda*, *Urophora stylata* develop within the capitula. Their larvae feed on achenes, pappus and/or receptacle; *U. stylata* forms a many-chambered gall.
- *Torymus chloromerus* and *Pteromalus elevatus* attack these tephritids and account for more than 98% of parasitism in that source system.
- **Tephritid adult activity was weeks 13–16; parasitoid oviposition was from week 17**, on third-instar fly larvae. Thus the consumer's life stage and plant stage must be measured before inferring that parasitoids benefit seed production.
- Fertilized potted thistles had approximately **2×** tephritid emergences and **4×** parasitoid emergences versus controls, but **treatment did not significantly change the proportion parasitized** in its reported model. More parasitoids do not automatically mean a larger proportion controlled or more seeds saved.
- The proportion of *X. miliaria* occupied buds was about 14% across a range of density in one potted plant experiment; this is specific to its head abundance design, not a fixed ecological infection probability.

Verified system 2: Masters, Jones & Rogers (2001), *Cirsium palustre* × *Terellia ruficauda* × *Pteromalus elevatus* / *Torymus chloromerus*, doi:10.1007/s004420000569. Root-herbivory manipulations change head size and abundances but percentage parasitism was similar across groups. Authors suggested parasitoid ovipositor accessibility could improve on smaller heads, without proving direct spine/angle access or seed rescue.

Verified system 3: Tofts (1999), *Cirsium eriophorum*, doi:10.1046/j.1365-2745.1999.00369.x (full text includes Table 2).

- **Floret layer:** flower-visiting bumblebees, Lepidoptera and pollen feeders. Visiting/feeding is not necessarily effective pollination.
- **Bract layer:** Vespidae adults remove wool and drink exudate from incisions in phyllaries. Therefore an apparently defensive involucral structure can itself be a directly exploited resource.
- **Head shelter layer:** *Forficula auricularia* shelters in seed heads. Shelter occupancy does not by itself mean either protection of seeds or seed damage.
- **Head larval-feeder layer:** multiple *Larinus*, *Terellia*, *Urophora*, *Xyphosia* and other taxa are recorded on or in flower/bud/seed heads, often from DIFFERENT source regions/years, not one local assemblage.
- **Host–parasitoid layer:** *Pteromalus vibulenus* is described as endoparasitoid of *Rhinocyllus conicus* in *C. eriophorum*. **The same Table 2 contains some parasitoids with three asterisks explicitly indicating occurrence from *C. palustre* only, NOT confirmed in *C. eriophorum*.** Do not count a compiled host association as a same-plant field observation.

New source role ledger: `data/evidence/eriophorum_microhabitat_role_audit_v1.csv`.

## Prospective novel question

> Does the phenological sequence of head presentation and defence cause a **switch in the ecological sign of the same structural module** between seed-feeder access and control by seed-feeder parasitoids, thereby favouring different combinations of orientation, spine geometry, bract posture, and stickiness across communities?

There are **three** distinguishable mechanisms:

1. **Conventional barrier:** the barrier reduces adult herbivore egg placement; fewer viable eggs and fewer damaged achenes.
2. **Enemy shielding:** the barrier suppresses parasitoid oviposition more than it suppresses established herbivore eggs, potentially raising total head damage **only if parasitoid-caused mortality prevents subsequent seed feeding**.
3. **Resource/phenology-only:** nutrient status, opening dates, head position/size or stage-specific morphology increases insect abundances on both sides without a barrier-specific interaction. Published arvense/palustre studies make this a strong competing baseline, not a straw man.

### Explicit mechanistic threshold — NOT a fitted estimate

Let `E(S,O)` = established seed-feeder eggs per comparable head, and `p(S,O)` = probability that an established seed-feeder is stopped by a parasitoid **before irreversible seed damage**. Let `d` = subsequent loss in seeds per damaging individual if it survives.

An idealized expected damage accounting relation is:

`D(S,O) = E(S,O) × (1 - p(S,O)) × d`.

If raising barrier S reduces established eggs by a factor `r_E = E_1/E_0 < 1`, the barrier still **increases** expected damage iff

`r_E × [(1-p_1)/(1-p_0)] > 1`.

Therefore a decrease in egg establishment can be outweighed by a large enough decrease in early effective parasitism. This is the actual sign-reversal condition and it is **not implied** by any literature-reported parasite abundance. The stage-specific `p` cannot be replaced by number of parasitoids, percent attacked after seed damage, or adult parasitoid visits.

Hypothetical arithmetic ONLY: if eggs fall to 70% of control, and effective pre-damage parasitism falls from 60% to 20%, expected damaging larvae increase by `0.70 × 0.80 / 0.40 = 1.40`, or **40%**, despite fewer eggs. The numerical example is not a thistle observation.

### Where orientation truly enters

Simply tilting a head **does not alter phyllary-spine geometry relative to florets**. It rotates the complete structure relative to gravity and insect approach paths. The mechanism predicts orientation × spine interactions **only when guilds have systematically different world-centred access paths or their access timing overlaps different head stages**.

- World-centred routes: parasitoids probe from the exterior, flies approach/land from above or lateral trajectories, stem-crawlers move upward and other consumers exploit exposed phyllaries.
- Head-centred null: all guilds solve the head geometry identically; orientation changes neither conditional entry nor effective pollination after stage/habitat are controlled.
- Head-position confound: terminal and lateral capitula differ in age, size and exposure; *C. pitcheri* seed predator and head-rank data show why a positional term is indispensable (Gijsman et al. 2020, doi:10.1016/j.gecco.2020.e00945).

## Ecological discriminators required BEFORE selection or adaptation language

1. The same plant/head identified by stable ID, with early bud, early/full anthesis, post-anthesis and mature-head records; head rank and position documented.
2. Actual structural traits: gravity-referenced orientation, direct botanical spine length/direction/rigidity, phyllary posture/access gaps and stickiness. The Azami 2-D projection is not a botanical defence measurement.
3. Consumer **role and time**: legitimate pollen transfer; seed-feeder adult oviposition, larval established counts, feeding onset and damage accumulation; parasitoid searching, ovipositor entry, true parasitism and timing of host death; remaining guilds (shelter, bract feeding, sap-feeding, predation).
4. Final viable seed counts; herbivore and parasitoid numbers are explanatory components, never replacements for plant fitness.
5. Baselines: ontogenetic stage, resource/phenology/head size and plant position, rain/wetness effects, and source-independent sampling bias.
6. Within-population natural variation and, only after safety/authorization, randomized factorial O × S external-access proxy with stage-specific observation. The existing JPN36 sham comparison is a **feasibility** test and may not be relabelled as an authorized multi-factor experiment.
7. Evolution claim requires independent heritable/phylogenetic comparison of trait coupling and regime-dependent selection. Absence of synchronized transitions or a low image-based RV alone cannot prove adaptive independent evolution.

## Discovery-source gate

- Full-genus GloBI API: `analysis/collect_globi_cirsium_interactions_v1.py`. It must return supported positive relations, track both Cirsium directions, preserve source citation and body-part metadata, separate refuted assertions and report complete vs truncated pages.
- Parasitoid second hop: `analysis/collect_globi_cirsium_seed_feeder_secondhop_v1.py`. The parasitoid–herbivore edge by itself is NOT evidence of the same head, plant or study; require independent original source co-localization.
- Frozen microhabitat sub-archive: `data/evidence/eriophorum_microhabitat_role_audit_v1.csv` preserves Tofts' provenance footnotes including **three-star parasite records exclusive to palustre**. Such rows must not be counted in confirmed *eriophorum* interaction degrees.
- No body-part context ⇒ `NOT_HEAD_LOCALIZED`; no stage ⇒ `UNRESOLVED_TROPHIC_TIMING`; no viable-seed endpoint ⇒ `FUNCTIONAL_EFFECT_ONLY`; no genotype/ancestry ⇒ `NOT_CORRELATED_EVOLUTION`.

## Current claim ceiling

**Evidence:** several *Cirsium* species contain real multi-trophic, microhabitat-structured assemblages with distinct adult/larval stages and parasitoid hosts; source-level counts and references are auditable, while GloBI live whole-genus retrieval remains unconfirmed.

**Not evidence:** an observed orientation × spine defence fitness interaction, a parasite's rescue of viable seeds, adaptive modularity, correlated evolution, or a GloBI-wide completeness claim.

**Novelty if tested:** the same organ architecture gates different mutualist–antagonist–enemy routes in different developmental windows; the sign of its net reproductive value can be predicted from guild access and temporal order, not simply the presence or absence of a consumer.
