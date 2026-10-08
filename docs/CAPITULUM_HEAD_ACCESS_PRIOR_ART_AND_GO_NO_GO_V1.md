# 頭花の形と生物間相互作用：今回の優先順位を凍結する

2026-10-08. Decision supplement to PR #44. No new field manipulation, no change to Azami/EAzami frozen endpoints. **Head only**. Scope: orientation, real phyllary armature/access gaps, capitulum internal dimensions, floret presentation and target guild/stage.

## Literature correction: parasitoid-free space inside a thistle head is not a new discovery

Romstöck-Völkl 1990 (*Ecological Entomology* 15:321–331), DOI [10.1111/j.1365-2311.1990.tb00814.x](https://resjournals.onlinelibrary.wiley.com/doi/abs/10.1111/j.1365-2311.1990.tb00814.x) already measured ovipositor lengths and host larval locations in *C. heterophyllum*/*Tephritis conura*/*Eurytoma* sp. near *tibialis* and *Pteromalus caudiger*. **Parts of Cirsium capitula form parasitoid-free structural refuges**; larval access differed markedly with head characters and location, and average accessibility did not vary systematically with host patch size. The same source reports considerable variation in laboratory parasitoid success driven by host position. This is **direct prior art**, not merely a theoretical possibility.

Maletti et al. 2021, DOI [10.1111/jzs.12433](https://onlinelibrary.wiley.com/doi/full/10.1111/jzs.12433) studied ovipositor-length variation within the *Pteromalus albipennis* group and discuss host association with head size. Thus the simple claim `larger/more compact thistle heads exclude short-ovipositor parasitoids` is also not intrinsically novel.

Peterson et al. 2016, DOI [10.3389/fpls.2016.01794](https://doi.org/10.3389/fpls.2016.01794), reviews plant structures as refuges from parasitoids and biological-control interference. The previous theoretical patch-migration reversal is interesting for sensitivity, but a direct, new universal theorem about indirect defense is **not** established.

## Missing mediator not yet measured in aza3

The current protocol separates adult arrival, attempted head contact, final oviposition, parasitoid search and seed output. A crucial missing within-head link is **where the damaging larva sits relative to the reachable exterior interface and the parasitoid's effective attack range**. An antagonist can enter successfully and nonetheless be protected from parasitoids by internal placement; adult access and larval refuge are not interchangeable.

For an independently identified host–parasitoid pair, the minimum *descriptive* accessibility geometry is:
- actual head stage, involucre thickness/access portals, floret/receptacle layers, gravity-calibrated orientation, head rank;
- host larval 3D/serial-section location, number, instar and prior seed damage, with destructive stage sampling appropriately timed;
- candidate parasitoid species, ovipositor length, allowed probe angle and **effective insertion depth**; nominal ovipositor length is an upper bound, not guaranteed penetration;
- outside accessible entry surface and the **shortest anatomically admissible probe path to a living larva**, `d_path`.

An idealized geometric necessary condition for physical attack is `effective_reach >= d_path`. But this is **not sufficient**: probing behaviour, mechanical stiffness, host detection, timing, sex, tissue obstruction, parasitoid survival and host defensive response also control success. `d_path` is not obtained from a flat image or unvalidated cylindrical head assumption.

This is directly testable without an inflated multiyear, genome-first project. The **same-life-stage** head dissection plus properly verified parasitoid taxon can classify accessible/inaccessible larvae (with uncertainty). Inferred fitness requires viable-achene outcome; destroying the same head for dissection precludes also measuring its final seed set, so use prespecified matched head cohorts and do not call that same-head mediation.

## Prioritization after correction

| Priority | Question | Strongest rival explanation | Required next data | Decision |
|---|---|---|---|---|
| 1 / empirical ecology | Does *real orientation × real phyllary/armature* change access of specific seed feeders and legitimate pollinators? | adult host race, head stage, plant display, alternate host blooms, image visibility, rain | pre-entry guild, continuous video with verified approaches and failed entries, real botanical head measurements, viable achenes | MAIN |
| 2 / tritrophic mechanism | Conditional on actual established larvae, is structural parasitoid accessibility explained by **probe path versus effective ovipositor reach**? | larvae choose different internal niches; host instar and timing; parasitoid behaviour | stage-specific 3D/dissection maps, true parasitoid length/reach, verified attack and larval fate | SECONDARY (prior art must be cited) |
| 3 / space-time scale | Can late parasitoid suppression reverse patch-level damage when generations return locally? | mixed dispersal, insect density dependence, outside immigration, pollen/seed limits | mark–recapture or lineage/parentage movement, next-generation pressure, independent populations | LONG-TERM HOLD |
| 4 / evolution | Do regimes of selective access explain repeated trait-module combinations? | phylogenetic history, shared image geometry, developmental allometry, neutral drift | independent ancestries and heritable morphology with measured population-specific selection | FUTURE GATE |

The previous **2-patch** migration audit is particularly decisive: under restrictive hypothetical `r=.7, K=1.2, q0=.6, q1=.1`, exchanging only ~14.3% of offspring between equal patches is enough to extinguish the theoretical *patch-level* reversal. Thus making delayed reversal the thesis headline now would depend on a **completely unmeasured dispersal assumption**.

## Small practical pilot in the existing ecosystem

The existing PolliPi/InsePi insect order classifier (Diptera/Hymenoptera/Lepidoptera) cannot on its own distinguish `Tephritis` seed feeder from other Diptera, or `Pteromalus` parasitoid from flower-visiting Hymenoptera. Positive per-episode pre-entry role requires morphology/experts, local natural-history evidence, linked video behaviour and, where permitted, host rearing.

The existing **JPN36 24-head feasibility ceiling** is not authorization for a multifactor experiment or enough to infer phylogenetic evolution. The realistic first objective is to estimate **observable approach-camera coverage** and verify at least one insect–head interaction route with unambiguous head stage. If adult seed-feeder/parasitoid species cannot be resolved, or neither enters the same head/time window, STOP the mechanistic contrast rather than coding zeros or pooling insect orders.

If verified pollinators and seed feeders take **different true portals**, the common 'pollination vs head defense' trade-off can be absent. If the same real outer barrier is contacted by a seed feeder and a parasitoid at different stages, prioritize time-specific accessibility and seed fate, not a high-dimensional correlation of 42 image taxa. If internal larvae exceed actual parasitoid insertion reach, that is a measurable *mechanism*, but the structural refuge itself was known by 1990.

## Manuscript claim ceiling today

- **SUPPORTED as literature:** *Cirsium* heads contain genuinely different biological portals and known parasitoid structural refuges; different insects use the head at different stages; prior art explicitly addresses ovipositor size and larval position.
- **SUPPORTED as original repository data:** real Azami image-phenotype correlation and nonuniform modular association, GloBI full-genus discovery records; no direct experimental function.
- **NOT supported:** thistle-spine-induced seed rescue, seed-damage reversal, phylogenetic co-adaptation, or selection on orientation–armature combinations.
- **Best prospective novelty:** when independently identified guilds see **different functional interfaces** of the *same modular reproductive head*, which head traits regulate which transition to effective pollination, seed feeding and parasitoid attack, and whether the empirically measured effect on viable seed production shifts across stages/regimes.

This decision note explicitly **demotes** generational reversal to a falsifiable conditional extension until local insect return and plant recruitment are measured.
