# 共有撮影時間で比較するための修正と、受粉履歴による未解決交絡

2026-10-08. Patch review for PR #44. **No real biological observations and no approved new field intervention.** This corrects *which recordings are comparable*, not the status of causal adaptive evidence.

## 1. Exact time-boundary equality incorrectly excludes usable recordings

Previous `analysis/audit_capitulum_concurrent_clock_guilds_v1.py` grouped videos by **identical UTC start and end timestamps**. Two continuous 15-minute recordings at 09:00–09:15 and 09:05–09:20 have 10 minutes of genuinely simultaneous recorded opportunity, but were placed in different, incomparable groups.

`analysis/audit_capitulum_common_coverage_bins_v1.py` repairs this using **preregisterable UTC-aligned nonoverlapping 5-minute bins**, not a sliding window selected after viewing visitors. A head contributes to a bin **only if all five minutes are entirely covered by valid continuous, camera-recall-checked video**. In that same head, event time is reconstructed from its authoritative `time_from_bout_start_s`, `video_window_start_s`, and the timezone-aware contiguous clip start. Events lying on partial, missing, event-triggered, or unverified intervals cannot become zero events or inflate overlap.

The 09:00–09:15 and 09:05–09:20 recordings share *two* valid 5-minute bins: 09:05–09:10 and 09:10–09:15, rather than being falsely labelled incompatible.

**Synthetic-only adversarial verification:** both independently identified insect guilds use the *same world approach route within each of these two common bins*. In the first bin, candidate pollinators/seed-feeders are 90/10 and all approach from above; in the second bin 10/90, all approach from below. **Marginal pooled route overlap = 0.2; within-shared-bin overlap = 1.0.** The ecological inference reverses without any biological change in how a guild approaches a head. These are mathematically constructed counterexamples, not Cirsium head observations or estimates.

## 2. Valid denominator and analytical selection are distinct

Five-minute bin inclusion depends on **continuous photographic effort**, not on insects showing up. The output must retain:
- fully covered bins where neither independently identified focal guild approaches;
- fully covered bins where only one of the two focal guilds approaches;
- bins where both occur, but individual head coverage, specimen identification, or replicate number is insufficient;
- bins eligible for **conditional descriptive** route overlap, requiring >=10 validated candidate approaches/guild, >=3 independent plants/guild, and >=3 **same heads actually visited by both guilds**.

The 10/3/3 limits are usability/overlap gates, **not statistical power or generalizability thresholds**. Restricting output to bins with both guilds arriving conditions on a downstream event and selects unusually favorable contexts. Consequently, even a route contrast from a matched bin **cannot estimate the magnitude of defence filtering for the entire Cirsium population**. Unobserved guild ≠ ecological absence if identity/effort is uncertain. Also report the full-coverage 5-minute exposure distribution and selection status.

The authoritative EAzami `aim2_capitulum_observation_bout_ledger_v1.csv` already has `observation_date`, `start_time_local` and `end_time_local` columns. A pinned native original-bout version (`c9870e91604ed1514344ddc8f765eae5237dbf27`) is compared to the new sidecar (which adds timezone and reviewed-video footprint); **both currently have zero organismal observation records**. No frozen data or field authorization is altered.

## 3. Even synchronized and phenophase-matched heads need NOT have the same scent

A newly emphasized direct original study matters here:

**Theis & Raguso 2005, *Journal of Chemical Ecology*, DOI [10.1007/s10886-005-7615-9](https://doi.org/10.1007/s10886-005-7615-9)** experimentally hand-pollinated heads of *Cirsium arvense* and *C. repandum*. Both species reduced total floral-scent emission **approximately 89% within 48 hours** after pollination. In the *C. arvense* study population, *Apis mellifera* was **nearly three times more likely to visit an unpollinated than a pollinated head**. *C. repandum* showed less pronounced discrimination by its main swallowtail visitor. This is real manipulation/observation of **post-pollination scent/visitation**, not a real physical-spine access result. The large scent change occurs at a finer physiological timescale than the broad `full_anthesis` label.

Thus identical clock/stage does **not** equate the floral advertisement history. If a head had been pollinated earlier, it could attract fewer insects despite identical exterior armature. Additional prospective record: for each head×bout, verified `prior_effective_pollen_deposition`, `time_since_evidence`, the assay/video source and `unassessed` distinct from negative. If not measured, report physiological history as **unknown**, not unchanged. Given no real data, do not expand the frozen EAzami source schema or fabricate prior pollen treatment; this is a required future study design check.

These prior works directly undermine any claim that stage/time-adjusted guild arrival differences alone identify structural trait effects:
- Theis 2006 *Journal of Chemical Ecology* DOI [10.1007/s10886-006-9051-x](https://doi.org/10.1007/s10886-006-9051-x): fragrance components can draw floral herbivores and pollinators;
- Theis, Lerdau & Raguso 2007 DOI [10.1086/513481](https://doi.org/10.1086/513481): floral scent and insect guild abundance vary by time-of-day and development;
- Vanbergen et al. 2007 DOI [10.1111/j.1365-2311.2007.00885.x](https://doi.org/10.1111/j.1365-2311.2007.00885.x): C. palustre host-parasitoid variation strongly associates with plant phenology.

## 4. Primary scientific target remains unchanged

**Conditional on the same real insect taxon, head stage, clock-time observation, effective prior pollination state, validated arrival opportunity and site, does independently measured true involucre/spine geometry × gravity-referenced orientation change physical access failure and confirmed oviposition more than pollen delivery?**

First-contact failures, rerouting, actual oviposition/pollen deposition and eventual filled/viable achenes must be observed; the last is the direct maternal-fitness proxy. If only candidate activity differs, the correct claim is **temporal/olfactory sorting**, not an anatomical access gate. If access differs but final seeds do not, call it an interaction filter, not demonstrated adaptation. Evolutionary reassembly requires multiple independent lineages, heritability/ancestry and regime-specific selection, none of which follows from image RV.

**Decision as of 2026-10-08:** The common-coverage code fixes a real logic bug. The previously frozen Azami 1,734 image cohort and the EAzami/GloBI original assays supply zero real matched, independently identified, continuous time–portal–viable-seed results. Keep the anatomical-adaptation claim **NOT IDENTIFIABLE** until observation; do not infer either a zero effect or a sign reversal.
