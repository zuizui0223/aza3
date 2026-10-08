# 頭花研究の次の識別障壁：昆虫が羽化した頭花だけを見る選択バイアス

Date: 2026-10-08. Added to PR #44, without altering any frozen Azami/EAzami input or field authorization.

## Original source — head-level observation versus experimental intent

Xi, Eisenhauer & Sun (2015), *Journal of Animal Ecology* 84:1103–1111, DOI [10.1111/1365-2656.12361](https://doi.org/10.1111/1365-2656.12361), section **Controlled Field Experiment**, tested *Saussurea nigrescens* × *Tephritis femoralis* × *Pteromalus albipennis* in enclosed transplant groups. The experiment assigned fly and parasitoid introduction at enclosure level, but the final head sample was **conditioned differently by treatment**:
- in the fly-only group, select only heads from which adult flies emerged;
- in the fly + wasp group, select only heads from which adult parasitoids emerged;
- other treatments were sampled randomly, per published methods.

The paper reported 8.04 damaged seeds/head under fly-only and 14.83 under fly + parasitoid within the analyzed sample. The independent microcosm supports prolonged host larval development, but **this does not transform the outcome-conditioned head sample into the average causal effect of parasitoid introduction over all initially eligible flower heads**. Individual-head emergence depends on attack, survival, host suitability and other potential causes of damage. In causal notation it is unsafe to use `E[Y | treatment, emerged type]` as `E[Y(1)-Y(0)]` when treatment changes emergence/inclusion. It would be equally wrong to conclude from this risk that the original effect is *entirely* an artifact. The experiment and supporting microcosm should be discussed precisely, not dismissed.

The same article includes *Cirsium setosum* in a **five-species observational survey** only. The species-specific result for the *proportion* of damaged seeds was marginal, rather than a randomized Cirsium spine/angle intervention. Its five-species field sampling also classified heads by emerged adults; this is a realized post-infestation association, not a controlled morphological treatment.

## Constructed numerical counterexample — strictly not the original experiment

Suppose two randomized arms each start with 100 identical heads: 50 potential high-damage (20 seeds lost) and 50 potential low-damage (5 seeds lost), with **no effect whatsoever of treatment on any head's true seed loss**. The average under either treatment is 12.5 seeds. Now sample differently by post-treatment insect emergence:

| Hypothetical analysis subset | High-damage heads included | Low-damage heads included | Mean damage |
|---|---:|---:|---:|
| Fly-only heads selected by emergence | 10 | 40 | 8.0 |
| Fly + wasp heads selected by emergence | 33 | 17 | 14.9 |

The selected-head ratio is **1.8625**, while the true unselected effect is **zero**. These chosen counts do not reproduce Xi's inclusion mechanism, and the resemblance to its ~85% damage increase is **not evidence** that selection caused Xi's difference. It proves only that the estimands are logically distinct.

Executable strict guard: `analysis/audit_capitulum_post_emergence_selection_v1.py`.

## Implications for Cirsium whole-capitulum architecture

**Do not preregister 'conditional on oviposition success or parasitoid adult emergence' as the only estimand.** It could make strong phyllary/spine defences appear associated with surviving enemies, because excluded/failed organisms never enter the analysis.

1. **Head-level total effect (priority):** register head/plant IDs **before** approach/oviposition; preserve all continuously observed heads including confirmed zero arrivals and failed approaches; follow planned matched intact heads to final filled/viable achenes. If feasibility/permissions allow randomized intervention, analysis follows randomized plant/head or enclosure assignment **without conditioning on emerged adults**; missing seed outcomes are audited and reported rather than assigned zero.
2. **Stage mechanism:** among independently identified pre-entry seed feeders, count physical contact attempts at the same true involucre/spine interface and confirmed entry failures. This is an access **conditional** on approaching a head and cannot itself provide an unconditional total-effect estimator when the traits influence approach.
3. **Parasitoid effect:** separated prospective timing of confirmed early host death versus continuing-host feeding, independent of adult emergence. Do not equate wasp emergence with seed protection.
4. **Robustness:** stratify exact phenophase, head rank, flower sex and clock overlap; track host/wasp species, alternative co-flowering *Cirsium*, rainfall and prior pollination state. Head/plant biological replicates rather than repeated approach episodes or images.

**The most important falsifier** is that same-stage, same-clock, same-preidentified species approaches encounter equal true outer-armature geometry but have equal entry failure or viable seed outcomes despite visible image trait associations. Under that result there is NO evidence for selection via a true spine access barrier; hypotheses about chemical attraction, post-pollination scent or abiotic head orientation remain valid alternatives.

## GO / HOLD

- GO to observational mechanism test only with actual documented anatomical spine or phyllary contact, complete verified video opportunity, independently identified insect role *before observed success*, and unsuccessful attempts retained.
- GO to head-fitness and causal comparison only with pre-outcome head sampling, filled/viable achene outcomes and a validated design.
- HOLD repeated coevolution/selection claims until independently derived lineages and heritable population variation with measured fitness.
- Existing authorized feasibility ceiling and previous genotype/trait analyses unchanged; **zero new real event/fitness rows**.
