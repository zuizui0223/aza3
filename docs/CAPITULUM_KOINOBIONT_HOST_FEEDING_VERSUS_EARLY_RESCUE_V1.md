# アザミ頭花の防御は、寄生蜂を遮断すると有利にも不利にもなり得る

2026-10-08 | PR #44 primary-literature correction + exact conditional identity.
**これは Cirsium の新しい野外効果量ではなく、検証可能性と先行研究上の上限を更新する文書。**

## 1. 発見ではないものを明示する：寄生蜂は「植物を守る」とは限らない

Xi, Eisenhauer & Sun (2015), *Journal of Animal Ecology* DOI [10.1111/1365-2656.12361](https://doi.org/10.1111/1365-2656.12361), **“Parasitoid wasps indirectly suppress seed production by stimulating consumption rates of their seed-feeding hosts.”**

原著の三種類の証拠を混ぜない。

| Source component | Focal botanical taxon | Actually shown | What it does NOT establish |
|---|---|---|---|
| 5-species observational capitulum comparison | **Includes *Cirsium setosum*** plus four other Asteraceae | Capitula with seed feeder + parasitoid emergence often had greater seed damage; in *C. setosum*, the damaged-seed **proportion** response was only marginal in species-specific tests | Did NOT manipulate true *Cirsium* spine or head orientation; parasitoid presence and seed-feeder condition were not randomized for *C. setosum* |
| Controlled field-enclosure experiment | ***Saussurea nigrescens*, NOT Cirsium** | *Tephritis femoralis* alone: 8.04 damaged seeds/capitulum; with *Pteromalus albipennis*: 14.83, ~85% higher; parasitoid-only did not directly damage seeds | Not a *Cirsium* experimental effect; note that capitula within fly/wasp treatments were filtered by insect emergence, which may affect generalization from experimental enclosures |
| Host developmental microcosm | *T. femoralis* larvae fed *Saussurea* seeds | Parasitoid-infected larvae prolonged development ~5 days, resulting in larger pupae and greater seed consumption | No direct measurement that involucre geometry gates the parasitoid or that these effects are shared by all guilds and all *Cirsium* |

This is **direct pre-existing evidence that a third trophic level can reduce, not increase, short-term plant seed fitness.** The earlier Cirsium-specific 1990 intracapitulum ovipositor-refuge literature (Romstöck-Völkl, DOI [10.1111/j.1365-2311.1990.tb00814.x](https://doi.org/10.1111/j.1365-2311.1990.tb00814.x)) remains additional prior art. Neither a refuge alone nor parasitoid-feeding stimulation alone is a new discovery.

## Source-data availability audit — not an individual Cirsium dataset

The public [Dryad record](https://doi.org/10.5061/dryad.m195j) for Xi et al. (2015) identifies three deposited spreadsheets: `Data_Figure1.xlsx`, `Data-Figure2.xlsx`, `Data_Figure3.xlsx`. The **published usage notes explicitly describe figure-level minimum, quartile, median, maximum and mean per treatment**, not observation-level insect identities, true spine measurements, exact parasitism/seed-damage timing or within-head approach events. The five-species `Cirsium setosum` figure's available summaries are therefore **NOT enough** to identify an individual-level morphological effect or to fit `q,p,g` to Cirsium. The file-stream retrieval returned HTTP 403 in this review; **no new numeric re-analysis of the spreadsheets is claimed**. An inaccessible file is not a negative biological result.

The original 2015 paper's field trial also selected individual capitula for some post-treatment groups based on fly/wasp adult emergence. This selection and the multilevel enclosure design must be retained when discussing generalizability; the reported +85% aggregate mean damage is not a per-larva causal coefficient.

## 2. Two functionally opposite parasitoid pathways

The former simplified same-head model `D = E(1-q)d` treats parasitoid action as removal of seed feeders **before irreversible seed damage**. It misses **koinobiont parasitoids**, whose living host continues feeding and may eat even more.

For a given head group, let:
- `E0>0` = baseline established damaging larvae per comparable head, `r=E1/E0 in [0,1]` = changed establishment after exterior barrier;
- `q0,q1` = independently measured fraction killed before causing irreversible achene damage (early beneficial control);
- `p0,p1` = of remaining larvae, the fractions successfully parasitized **while they continue seed feeding**, NOT late adult parasitoid abundance;
- `d>0` = damage by an unparasitized feeding larva, `g>0` = damage multiplier of a parasitized larva relative to one unparasitized larva, here held fixed across barrier states for identifiability.

Then a limited expectation identity:
```
D0 = E0 * d * (1-q0) * [1+p0(g-1)]
D1 = E0 * r * d * (1-q1) * [1+p1(g-1)]

D1/D0 = r*(1-q1)/(1-q0) * [1+p1(g-1)]/[1+p0(g-1)]
```
when `q0<1`, E0,d positive. `q0=1` makes baseline damage zero and ratio undefined. This is a **bookkeeping identity under restricted assumptions**, not a fitted causal model; density dependence, compensatory plant reproduction, ovule potential, parasitoid guild mixtures, time of attack, altered rates of pest oviposition and differential host/parasitoid movement are omitted.

**Crucial biological contrast:** If `g>1` and a barrier excludes late parasitoids `p1<p0`, it *reduces* the damage added by parasitoid-extended host feeding. This makes blocking a parasitoid potentially **beneficial** for plants in the same seed season. In contrast, reducing early effective parasitoid killing `q1<q0` can increase damage. Their net effect is indeterminate without both stage-specific channels.

A purely hypothetical example (`r=.7, q0=.6, q1=.2, p0=.6, p1=.1`):
- when `g=1`: damaged-seed ratio `D1/D0=1.4` (**barrier apparently harmful**, because it stops early helpful parasitoid killing);
- when `g=2`: `D1/D0=0.9625` (**barrier apparently beneficial**, even though its early parasitoid effect is unchanged).
- exact hypothetical sign boundary: `g=1.869565...`, assuming *all other parameters fixed*.

**These parameter values do not come from the Xi paper.** In particular its ~85% difference in mean damage **per capitulum** is **not** an empirical estimate of `g`, damage per parasitized larva, in any *Cirsium* species.

This calculation and an extensive no-real-data gate are in `analysis/check_capitulum_parasitoid_feeding_sign_v1.py`.

## 3. How the ecological result could explain head modules, if observed

The interesting Cirsium-specific question is **not** simply "does more parasitism help plants?". It is:

> Under the same head stage, ambient guild pressure, head orientation and actual phyllary/spine architecture, does barrier geometry preferentially admit (a) egg-laying seed feeders, (b) early-killing parasitoids, or (c) late koinobiont parasitoids? Does the resulting change in real host feeding and seed production make one structural combination better than another?

This can generate three falsifiable outcomes:
1. **Selective barrier**: same-species candidate seed feeders contact real involucre/spines but fail to oviposit more often, while effective pollen deposition is unaffected; actual filled/viable achenes increase.
2. **Helpful-enemy exclusion**: barrier stops parasitoids which would have killed damaging larvae before seed loss, driving higher damage even with fewer established seed feeders.
3. **Harmful-enemy exclusion**: barrier stops koinobiont parasitoids that would otherwise prolong larval feeding, reducing seed damage and potentially benefitting plants despite blocking an enemy of a pest.

Failure to detect any of the above when route-identified insects do not contact the real structural surfaces supports a **no realised physical conflict** explanation. Body length, head-size allometry, plant display, microclimate, phyllary toughness, larval internal position and floral scent remain alternative pathways.

## 4. Most conservative empirical milestone

Before any further multigeneration spatial model, record the causal sequence **by identified taxon or host–parasitoid pair**, never order-level automatic insect classification alone:

`continuous video effort → independently classified arrival → actual involucre/true-spine touch → physical entry attempt and failure/avoidance → confirmed oviposition / parasitoid probing → timing of larval feeding/damage → viable achenes`.

The current event ledger can record observed parasitoid searching and host attacks, but **cannot equate search or adult wasp emergence with early host death, nor classify a koinobiont's per-host subsequent feeding**. Retain `NA` whenever timing is not independently linked, and separately document whether recorded 'obstruction' had physically verified spine contact or was behavioural avoidance. Future stage-specific larval tracking and mature seeds may require different matched destructive/non-destructive head cohorts, not falsely claimed same-head mediation.

**What can be claimed today:** The original Xi et al. paper directly establishes that the sign of parasitoid effects on seed production can be negative, including field experimental support in **Saussurea** and limited observational *Cirsium setosum* evidence. This makes a one-direction "parasitoid access benefits thistles" prior biologically unsafe.

**What is NOT claimed:** new field effects of thistle armature or orientation; a fitted `g`; evolution of a head module; or any measured fitness-sign switch in the existing **zero-row** current field intake. Existing Azami/EAzami frozen results and permission ceilings remain unchanged.
