# Chapter 3 conservation and P6 acquisition audit v10

This document sits between the species-level hole-fill ledger and exact locality freezing.

The nationwide tree requires **representation slots**, not a fixed number of destructive wild collections. Acquisition mode is chosen from taxon status and evidence before a field site is frozen.

## National Red List screen

Source: Environment Ministry, 5th Red List, vascular plants (2025). The source table is frozen in `data/planning/chapter3_national_redlist_cirsium_source_v10.csv` and retains full-species versus infraspecific scope.

The 128-species operational core currently resolves to:

- **101** `STANDARD_WITH_LOCAL_REVIEW`;
- **13** `CONSERVATION_GATE_REQUIRED` (full-species EN/VU or accepted authority crosswalk);
- **6** `HERBARIUM_OR_MINIMAL_AUTHORIZED_WILD` (full-species CR);
- **5** `CONSERVATION_REVIEW` (full-species NT or species slot with NT plus an infraspecific higher-risk taxon);
- **2** `LOCALITY_AND_SUBTAXON_GATE` (`C. japonicum`, `C. tashiroi`);
- **1** `HISTORICAL_MATERIAL_ONLY` (`C. toyoshimae`, EX).

`NOT_LISTED_AT_CORE_SPECIES_LEVEL` means only that no national Red List hit is assigned to the operational species slot after the current authority crosswalk. It **does not** mean unrestricted collection. Prefectural red lists, protected-area rules, designated plants, land-manager restrictions and collection/tissue permits are still checked for every slot.

### Important crosswalk rule

Do not turn an infraspecific listing into a blanket species-level threat category. For example, the current `C. japonicum` species slot contains Red List signals for `var. australe` (CR) and `var. ibukiense` (NT); the target population/subtaxon must be resolved before acquisition mode is finalized.

Conversely, spelling/rank differences that clearly refer to the same current NMNS species slot are crosswalked explicitly rather than silently lost. Current examples include the NMNS/RL spellings around `C. ashinokurense`, `C. pseudosuffltum`, `C. austrokiusianum`, and `C. yakushimense`.

## Queue consequence

The exact-site queue still contains **229 execution rows** (228 required tree-representation slots + one optional P7 trait-link slot), but the latest full-screen build gives:

- **180** rows with no national nonstandard gate yet (`NOT_FROZEN`; local review still required);
- **47** rows at `CONSERVATION_GATE_REQUIRED`;
- **2** `C. toyoshimae` representation rows at `WILD_COLLECTION_BLOCKED` and historical/herbarium-only mode.

These counts are slot counts, not species counts. A two-representative species contributes two queue rows.

## Acquisition hierarchy

1. `HISTORICAL_MATERIAL_ONLY`: authenticated herbarium/historical DNA or compatible existing sequence only; no wild collection slot.
2. `HERBARIUM_OR_MINIMAL_AUTHORIZED_WILD`: seek authenticated herbarium/ex situ material first; live wild tissue only if conservation authorities explicitly approve a minimal sample.
3. `CONSERVATION_GATE_REQUIRED`: current occurrence plus explicit conservation/land/collection/tissue review before any live sampling.
4. `CONSERVATION_REVIEW`: complete local and protected-area review before normal field-freeze logic.
5. `LOCALITY_AND_SUBTAXON_GATE`: resolve which infraspecific entity/population is being sampled before applying the relevant conservation mode.
6. `STANDARD_WITH_LOCAL_REVIEW`: still requires current occurrence, local/prefectural status, land manager and permits.

## P6 — public singleton plus own complement

P6 has 26 species. The own sample is meant to add nonredundant geography and individual-linked phenotype, so public provenance must be known before an opposite/complementary range sector is chosen.

The retry-safe BioSample audit is integrated with the frozen Moreyra/voucher context. Current best evidence is:

- **8** `PRECISE_SECTOR` public localities;
- **16** `REGION_ONLY`;
- **1** `COUNTRY_REGION_INCOMPLETE`;
- **1** `COUNTRY_ONLY`.

After conservation integration:

- **6** species can proceed to complementary-sector selection after current-occurrence and permission audit: `C. ishidzuchiense`, `C. yezoense`, `C. maritimum`, `C. spicatum`, `C. gyojanum`, `C. dipsacolepis`;
- **4** species are conservation overrides where a simple opposite-sector wild sample is not the next action: `C. magofukui`, `C. pseudosuffltum`, `C. nippoense`, `C. aidzuense`;
- **1** (`C. japonicum`) needs subtaxon/locality review despite a precise public locality;
- **15** remain geographically unresolved and must not be assigned a complementary sector by guesswork.

For an unresolved P6 public singleton, the fallback is voucher/supplement/herbarium reconciliation. If that fails, use a two-own-sector sensitivity design rather than pretending the public accession supplies one geographic side.

## Public locality safety

The audit records whether a public BioSample contains coordinates but does not copy those coordinates into this public repository. Precise coordinates and rare-species microhabitat directions remain outside the public repository in a private field ledger linked by deidentified locality IDs.

## Reproducibility

- `analysis/build_core_conservation_screen_v10.py`
- `analysis/audit_p6_biosample_localities_v10.py`
- `analysis/build_p6_best_evidence_v10.py`
- `.github/workflows/build-conservation-and-p6-audit-v10.yml`
- `analysis/build_exact_site_freeze_queue_v8.py`
- `.github/workflows/build-exact-site-freeze-queue-v8.yml`

The current workflows are fail-closed: no acquisition is authorized merely because a species or sector appears in a planning ledger.
