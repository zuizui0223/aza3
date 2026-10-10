# 真の頭花の防御を判定する最短の実測単位：トゲ接触→阻害→侵入失敗

Date 2026-10-08. Add-on to PR #44. This is an **empty prospective observational contract**, not authorization for manipulating or collecting Cirsium.

## 今回の推論ギャップ

The existing `capitulum_guild_access_event_ledger_v1.csv` tracks independently preidentified candidate guild, approach direction and final access. The `capitulum_contact_attempt_sequence_v1.csv` records sequential attempts and an observed 'obstruction'. Neither guarantees that the insect **actually touched an authentic botanical spine**, as opposed to a phyllary lamina, an involucre gap, or no structural tissue at all. Correlations of image silhouette projection with flower color cannot be used as measured `spine_length_mm`.

The append-only `data/intake/capitulum_true_armature_contact_v1.csv` and source validator `analysis/validate_capitulum_true_armature_contact_v1.py` enforce:
- Same population / plant / head / reviewed video bout / independent approach episode / ordered attempted-contact ID as the original two ledgers;
- Pre-entry candidate guild/taxon identified from independent evidence, not assigned after a successful oviposition or pollen deposition event;
- An independently reviewed **true anatomical contact surface**: spine tip, spine shaft, phyllary lamina, phyllary interstice, floret disc, stem, no surface contact, or unresolved;
- Positive spine touch only for actual spine tip/shaft contact, and **NA for unknown** rather than 0;
- `physical_blockage_verified=1` only when true armature was touched, the parent continuous episode independently reported actual mechanical obstruction, access failed at that attempt, source video documents it, and a calibrated *relevant* morphological dimension is measured. Merely not entering is NOT proven mechanical blockage;
- Duplicate individual-insect approaches, multiple attempted contacts, and repeated heads are **nested**, never independent plants or proof of statistical selection;
- A complete original **video effort ledger** is needed to estimate arrival rate or genuinely absent episodes. This annotation alone is not a population denominator.

The real CSV starts with **zero biological rows**. Synthetic deliberate corruption tests reject a false spine contact, fabricated anatomical scale, claimed physical barrier despite unverified entry failure, unknown treated as a true zero, missing evidence, post-outcome insect guild labelling, absent parent attempted contacts, and non-continuous video. These tests do not claim any detected head barrier.

## Target biological alternatives

**H_structure:** In comparable same-stage, same-time verified arrivals from the same identified insect taxon, authentic phyllary/spine contact increases *mechanical* failed access compared with independently measured different armature/gap states. Real pollen deposition and confirmed egg outcomes then need independent evaluation. Without an exogenous or quasi-experimental contrast this is primarily observational association, not natural selection.

**H_olfaction_and_time:** Different prior pollination histories, chemical cues, floral displays, clock windows or host races cause the relative frequencies of insect arrivals and initial head choices. Direct true spine touching and mechanical exclusion **do not differ** after a legitimate independent contact denominator is defined. This mechanism is compatible with apparent image-angle correlations.

**H_phenophase_avoidance:** The floret portal becomes accessible at anthesis and insect success arises from bypass/rerouting rather than direct outer armature selection. The same species can use a bud phyllary route and an anthesis floret route; do not presume fixed taxon identity dictates one portal.

**H_enemy_sign:** If confirmed parasitoids attack an already established larva, their effects must be split between true early killing before irreversible seed damage, and koinobiont-induced **extra host seed feeding**. The original *Cirsium setosum* observational association from Xi et al. (2015, DOI [10.1111/1365-2656.12361](https://doi.org/10.1111/1365-2656.12361)) is NOT an intervention on a Cirsium spine. Experimental *Saussurea* head means from that paper are conditioned on emergent insects in the experimental insect-containing arms: do not turn them into unselected population treatment averages or per-larva `g` for Cirsium.

## Decision rule and field feasibility

This contract does not require a huge 4-way experiment to start. First verify the ecological encounter **with one real insect–head contact episode and independently assessed anatomical structure**; then accrue naturally varying, properly synchronized heads within the same population/stage.

1. Confirm full approach-zone/whole-head continuous-video coverage, with a reviewable unambiguous stage and independently taxon-identified insect candidate. A Raspberry Pi order-level Diptera/Hymenoptera/Lepidoptera classification is *not enough* to decide pollinator, seed feeder or parasitoid roles.
2. Confirm true spine tip/shaft or other phyllary contact in the video. For genuine morphology, measure true `spine_length_mm` and/or `minimum_phyllary_gap_mm` with calibration; avoid claiming an inaccessible micro-spine is absent if video resolution is inadequate.
3. Count actual **failed and successful contacts**, explicit unresolved episodes, the same insect's repeated approaches (when identity can be ascertained), and the total complete observation effort. Collect final seed data from permitted and appropriately identified matched heads, not only heads with insects that emerge.
4. **STOP adaptive claims** at this stage if no physical contact occurs, if most approaches have unresolvable contact, if guilds cannot be independently identified, if time/stage exposure does not overlap, or if final viable achenes are not measured.

### Interpretation ceiling

A validated mechanical-contact row establishes **an anatomical encounter and observed failure**, not causal exclusion by an evolved trait. A causal barrier contrast requires matched or independently randomized geometry states with actual insect arrival exposure. Functional adaptive significance requires adequate head/plant replication and viable maternal seed output, with pollen function and post-pollination scent controlled. Repeated adaptive module assembly additionally needs independently resolved evolutionary origins and a credible evolutionary null. None of those can be inferred from the empty contract.

## 2026-10-10 denominator gate: blocked-only annotation is invalid

An important bug in the first version was that the source validator required every **submitted annotation** to match a parent attempt, but did not require every eligible **parent attempt** to be annotated. If a scientist annotated only spectacular blocked insects and ignored successful attempts, the result could falsely suggest that every anatomically contacted insect was obstructed.

`analysis/validate_capitulum_true_armature_contact_v1.py` now builds the authoritative eligible denominator from **all continuously reviewed parent attempts with independently qualified pre-entry guild**, not from annotated positive events. It rejects any incomplete annotation through `HOLD_INCOMPLETE_ELIGIBLE_ATTEMPT_ANNOTATION`, returning `null` for verified-contact and physical-blockage outcome counts until all eligible attempts have an annotation. Each valid annotation may explicitly say `no_surface_contact`, `floret_disc`, `unresolved` or other agreed surface; a 'no blockage' row may only have `0` when actually assessed. This also prevents mistaking unknown pre-entry taxa, unreviewed video and absent attempts for negative observations.

A second new gate checks the original *attempt portal* against the anatomically resolved site: a parent recorded as `floret_disc` may **not** be retroactively labelled `spine_tip`; this is a documentation conflict, not evidence of thorn contact. Attempts where recorded parent world/head portal was unresolved are still limited to a descriptive, documented video re-review.

**Synthetic falsifier**: one independently identified insect first contacts an outer spine and fails to enter, then makes a second documented attempt at the floret disc and succeeds. If only the spectacular spine failure is annotated, report HOLD and `n_unannotated_eligible_parent_attempts=1`, **not a barrier rate of 100%**. If both attempts are annotated, the schema can report two attempted contacts, one measured mechanical blockage, and exactly one biological head. This still is not an independent morphological treatment contrast or a seed-fitness estimate.

This denominator guard concerns **eligible observed contact attempts**, not overall insect arrival. Population-level effort and possible false-negative arrivals remain controlled separately by the continuous-video effort ledger. Even a complete observed attempt ledger can be biased if insects are missed, if pre-entry taxa are preferentially identified only for easy successes, or if the plant's geometry itself changes arrival. The valid next step is only a prospective **access process** description.

## Biological source-based context

*Gijsman et al.* (2020), `Cirsium pitcheri`, DOI [10.1016/j.gecco.2020.e00945](https://doi.org/10.1016/j.gecco.2020.e00945), found that heads with weevil infestation had 60% fewer mature seeds, and weevil attack was more frequent in secondary and tertiary heads than terminal ones. That is direct seed-fate leverage but **not** a direct spine-length effect. Rank must therefore be kept in comparisons; a successful head-annotation validator cannot rescue missing head-rank matching or missing viable-achene outcomes.

**Final ceiling:** no direct field record has been submitted to these current intake CSVs. A PASS fixture means only the software knows how to reject certain invalid biological claims; it does not mean thorn selection, modular adaptation, or antagonistic-access tradeoffs have been observed.
