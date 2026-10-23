# aza3 Q2 — the first operational head-access intake after real morphology and seed studies

**2026-10-23.** Clean mainline port from [unmerged PR #44](https://github.com/zuizui0223/aza3/pull/44), *not* a merge of that large stacked branch. Two independent, already formulated public-research observational contracts plus their source validators are now available directly in aza3/main.

## Biological reason for the intake

The observed evidence is **not one coherent plant-specific morphology→interaction→fitness cohort**:

* The real 114 source plants in [PR #85](https://github.com/zuizui0223/aza3/pull/85) contain internal-size-matched pairs with distinct external phyllary-tip spread (including 13 nonoverlapping pairs when body width is matched). **No filmed seed-feeder encounters, pollen assay or achene outcomes** are present in those individuals.
* The 2022 original Cirsium arvense 615-head study in [PR #82](https://github.com/zuizui0223/aza3/pull/82) shows bagged heads often lose viable seed and emerged insects together (25/36 complete site-blocks), but never measured true spine and head gravity orientation. The composite physical bag is not a selective pollinator-only intervention.
* The 2007–2008 Cirsium undulatum USDA 105 coherent mothers in [PR #84](https://github.com/zuizui0223/aza3/pull/84) show apical and lateral seed-success routes can jointly cover different observed maternal failures, yet five of twelve apically-damaged plants under ambient herbivory produced no scored good seed in either location. That is **not a manipulative test of route insurance, phyllary protection or genetic re-expression**.

The **missing biological link** is same-prespecified-head true outer morphology → independently identified pollinator/seed-feeder approach including failures → physical entry and realized service/damage → final viable achenes.

## Files actually made executable on main

* `data/intake/capitulum_video_effort_denominator_v1.csv` (header only): unique plant × head × bout, complete time range, actual continuous camera coverage, explicit verified zero-approach bouts, video evidence, other Cirsium head resources within a census radius.
* `data/intake/capitulum_guild_access_event_ledger_v1.csv` (header only): every observed approach episode with a guild assigned **before entry outcome** or recorded unresolved; world-frame vs head-frame routes; zone contact, pollen assay and oviposition separately; parasitoid host killing before seed damage only with the required timing evidence.
* `data/contracts/capitulum_video_effort_denominator_v1.json`, `data/contracts/capitulum_guild_access_event_contract_v1.json`: original validated semantics, keys and missingness. The related EAzami head/phenology tables remain separate existing records, not synthetic local head observations.
* `analysis/validate_capitulum_video_effort_denominator_v1.py`, `analysis/validate_capitulum_guild_access_event_ledger_v1.py`: source-port validators, exercised with synthetic positive and falsification checks. **Synthetic episodes are not entered as observations.**

The scripts reject unverified zero approaches, event-trigger-only footage reported as continuous, camera blind spots, false organism/head IDs, guessing a pre-entry guild based on successful oviposition, and a land/flower-contact event mislabelled as pollen deposition. A 5-second event-triggered camera clip alone cannot estimate true insect *approaches per unit filming time* without independently audited missed detections. Never use machine-labelled Lepidoptera/Hymenoptera/Diptera as a proven pollinator or seed feeder.

## Mainline verification gate

The CI runs both validators on the **actual header-only intake** plus their internal synthetic tests. It MUST report:

* real video bouts **0**, real approach episodes **0**;
* unconditional arrival rate **null / not identifiable**, **not** zero arrivals;
* no seed fitness, anatomical physical barrier or genetic LEGO result;
* synthetic falsification tests pass while no synthetic rows enter source CSVs.

A future observed head can be admitted **only after field access, sampling/manipulation and conservation permissions are individually checked**. Each head needs measured gravity orientation, true spine length/rigidity, involucre gap and stage/sex, identifiable effort and encounters, independently assayed pollen/oviposition and its own eventual filled/germinable achene fate. Validated approach rates are observational even then; the trait's selection benefit additionally requires legitimate experimental or quasi-experimental identification.

**Scientific goal:** distinguish a true architecture-by-partner *access filter* from signal-only attraction, phenological/olfactory cues, background Cirsium host availability and simple camera visibility. The current datasets alone cannot answer that causal question.
