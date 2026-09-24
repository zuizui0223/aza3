# Chapter 3 transition-first sampling plan v3

Status date: 2026-09-02

This is the authoritative new-sampling plan for Chapter 3. It supersedes the zero-baseline v2 acquisition logic where v2 still treated Japan38 as the practical discovery frame. Earlier 38×3, eight-anchor RAD-pilot and inventory-first designs remain as design history only.

## 1. Starting point: Japan38 is not the sampling universe

Chapter 2 was defined on the 38 Japanese paper concepts sampled by Moreyra et al. (2025). That set is the **hypothesis-origin subset**: it is where the current P01–P05 uncertainties were discovered. It is not a complete inventory of Japanese *Cirsium* and is not the default Chapter 3 field frame.

The current NMNS `日本のアザミ` index returns 161 authority records. Those are database records, not 161 guaranteed independent species: varieties, synonyms and other infraspecific/duplicate concepts must remain explicit until taxonomic reconciliation. Moreyra also notes that the Japanese archipelago harbours more than 100 *Cirsium* species and sampled only 38 of them.

Therefore the Chapter 3 acquisition order is:

> authority-wide trait screen → phylogenetic placement where needed → transition-neighbourhood selection → population RAD

not

> Japan38 → three individuals each → one all-Japan RAD matrix.

## 2. Two kinds of overlap

Sampling overlap has two different meanings and must not be confused.

### Efficient overlap — retained

The same voucher-linked individual should score orientation, phyllary posture, stickiness, colour and cytotype. The same JPN36 individual can therefore inform P01, P03 and P04. This is desired reuse, not pseudo-replication, because P04 explicitly asks whether module histories co-localize on the same individual/population ancestry.

### Evolutionary-replicate overlap — restricted

For P03, repeatedly sampling one local clade does not establish breadth across recurrent orientation histories. P03 therefore requires at least **four non-overlapping candidate transition neighbourhoods** before a broad recurrent-history claim is attempted. JPN36 or any other single neighbourhood may occupy at most one of those four slots.

`non-overlapping` means distinct focal transition edges/neighbourhoods on the admitted backbone, not merely different collection sites.

## 3. A0 — authority-wide digital screen, no field campaign

Do not physically survey 38 or 161 entries first.

Build a reproducible candidate universe from the NMNS index and retain, per authority record:

- authority taxon concept and Japanese name;
- infraspecific status;
- distribution summary;
- categorical orientation from the authority catchphrase where explicit;
- categorical phyllary posture where explicit;
- categorical stickiness where explicit;
- a source-text hash rather than redistributing the catchphrase;
- whether the concept is already represented in the Moreyra/Comp1061 scaffold;
- whether direct individual/population observations are still needed.

The authority state is a **screening prior**, never a same-individual measurement and never an ancestral-state reconstruction.

A0 generates two pools:

1. **backbone candidates** already represented in Moreyra/Comp1061;
2. **placement candidates** absent from the current backbone but potentially informative because their authority state would bracket or add an otherwise missing trait contrast.

Geographic proximity alone cannot promote a taxon to `close relative`. Non-backbone candidates must first be phylogenetically placed.

## 4. A1 — fixed focal population layer from the locked P01/P02 + M01 questions

These systems are sampled regardless of the later authority-wide screen because they already carry frozen Chapter 2 discrimination value or the embedded M01 depth role.

Default population bank:

- **12 primary RAD individuals**;
- **3 preregistered reserves**;
- 15 physical plants per population.

All 15 receive the same individual ID and phenotype/voucher package; reserves are not substitutes for biologically inconvenient genotypes.

### Minimum defensible launch

| system | populations | physical | primary RAD | reason |
|---|---:|---:|---:|---|
| JPN36 *C. sieboldii* | 3 | 45 | 36 | P01 phyllary + P03/P04 overlap; within-species polymorphism must remain visible |
| JPN06 *C. dipsacolepis* | 2 | 30 | 24 | P02 nonsticky side |
| JPN15 *C. lineare* | 2 | 30 | 24 | P02 sticky side |
| *C. brevicaule* | 2 | 30 | 24 | M01 discovery: Okinawa Honto + Amami Oshima |
| *C. irumtiense* | 2 | 30 | 24 | M01 discovery: Miyako + Ishigaki |
| **total** | **11** | **165** | **132** |  |

### Target focal replication

Before treating mainland population structure as well represented, target:

- JPN36: **4 populations**;
- JPN06: **3 populations**;
- JPN15: **3 populations**;
- *C. brevicaule*: 2 discovery populations at this stage;
- *C. irumtiense*: 2 discovery populations at this stage.

This yields **210 physical plants / 168 primary RAD individuals**. The number is a consequence of population replication, not a Chapter 3 quota.

M01 floral RNA/pigment sampling remains nested within the same four discovery populations: collect six developmentally matched floral samples/population, sequence five primary RNA libraries and retain one reserve.

## 5. A1-RAD — focal assay pilot is nested inside biological sampling

Do not collect unrelated high-ploidy taxa solely to test enzymes before their biological inclusion is justified.

Initial assay systems:

- JPN36;
- JPN06;
- JPN15;
- *C. brevicaule*;
- *C. irumtiense*.

Stage A-focal:

- 5 systems × 1 qualified DNA × 3 candidate enzyme protocols = **15 shallow libraries**.

Stage B-focal:

- 5 systems × 2 independent biological templates × top-2 protocols = 20 primary libraries;
- +5 independent-batch technical repeats = **25 libraries**.

Retain the previously frozen technical thresholds:

- callable overlap genotype concordance ≥0.95;
- core-locus recovery ≥0.90;
- depth target chosen from a locus-recovery saturation curve;
- trait outcomes cannot tune enzyme, Stacks or missingness settings.

This pilot qualifies only the focal cytotype/genome-size stratum. If later transition neighbourhoods introduce tetraploid/high-ploidy material, run a second stratum-specific assay check before production RAD.

## 6. A2 — backbone augmentation for informative non-Moreyra taxa

A Japanese species absent from Moreyra may improve transition discrimination, but its phylogenetic role is unknown until placed.

A0 may nominate an initial batch of **up to 12 non-backbone authority concepts** for placement. Nomination requires all of the following:

1. an explicit authority trait state relevant to P01/P03/P04 or a directly documented polymorphism;
2. a state contrast not already redundantly bracketed by existing backbone tips;
3. field/voucher identity that can be established;
4. no claim of close relationship based on geography alone.

For each nominated concept, collect **two biological placement representatives** where feasible. Use a Comp1061-compatible target-capture assay (or a separately justified compatible nuclear placement assay), not RAD, to attach the concept to the cross-species scaffold.

A narrow endemic may use two individuals from one verified population; otherwise prefer two populations.

A non-backbone concept cannot enter population RAD as a transition comparator until its placement is stable enough to define which candidate transition neighbourhood it belongs to. Ambiguous/hybrid/polyploid placements remain an explicit network/uncertainty result rather than being forced into a sister pair.

## 7. A3 — transition-neighbourhood selection

After A0 and any required A2 placement, select sampling units by **history discrimination**, not by species count.

### P01 phyllary posture

Keep JPN36 population replication as the fixed focal test. Add neighbouring comparator taxa only where the admitted nuclear scaffold shows that their inclusion distinguishes terminal-JPN36 from internal/alternative-edge histories. Do not add taxa merely because they share a phyllary state.

### P02 stickiness

JPN06–JPN15 remains the fixed ancestry-matched contrast. Population replication is the primary new information. Additional comparator species require a specific network/topology ambiguity they resolve.

### P03 orientation

Select at least **four non-overlapping candidate transition neighbourhoods** from the admitted expanded backbone.

Each neighbourhood must:

- contain an explicit orientation-state contrast;
- bracket a different candidate transition edge/neighbourhood from the other slots;
- include the closest placed comparator lineage(s) needed to distinguish terminal versus deeper placement;
- avoid counting two populations of one species as two evolutionary-transition slots;
- avoid counting JPN36 more than once.

Because Chapter 2 did not assign posterior probabilities to equally parsimonious histories, v3 does **not** describe this as “80% transition probability mass.” Candidate neighbourhoods are selected by discrete history discrimination across the admitted topology/minimum-history ensemble.

### P04 cross-module linkage

P04 reuses A1 and A3 individuals and deliberately measures all three modules on the same plants. No separate P04 taxon quota is opened unless the reconstruction-aware null shows an identified coverage gap.

### P05 continuous phyllary dimensions

P05 is measurement-first. Standardized calibrated images and same-voucher dimensions are added to A1/A3 plants. RAD coverage is not expanded merely to rescue a continuous-trait sample size.

## 8. Population sampling inside an added transition neighbourhood

For a non-focal comparator taxon that reaches A3:

- default: 2 populations × (10 primary RAD + 2 reserve) = **20 primary / 24 physical per taxon**;
- narrow endemic: one verified population may be used but the single-population limitation remains explicit;
- within-species polymorphism: prefer ≥3 populations and ensure ≥10 primary individuals per observed state across ≥2 populations where feasible, so morph state is not identical to locality.

Individuals are spread across the occupied patch; use ≥10 m spacing where biologically sensible and avoid obvious connected/clonal ramets.

## 9. M01 island expansion is still conditional

Do not collect Iriomote, Yonaguni and extra Amami populations merely to maximize island coverage.

Open the M01 expansion only if:

- ancestral-state/pigment history is interpretable enough to distinguish H1/H2 from H3 at the stated ceiling; and
- the corolla expression association replicates across both discovery populations within each lineage.

Then add one intermediate Amami-group population, one southern Amami-group population, Iriomote and Yonaguni. Default bank remains 12 primary +3 reserve/population.

RAD remains the population/background layer; candidate-region selection evidence requires an independent non-RAD confirmation route.

## 10. Cross-species target capture and RAD have different jobs

- Moreyra/expanded Comp1061-compatible target capture: cross-species placement and scaffold;
- RAD: shallow population ancestry, admixture, reticulation and morph/population history inside qualified neighbourhoods;
- same-individual phenotype: transition localization;
- M01 RNA/pigment: molecular decomposition;
- candidate resequencing/WGS/capture: selection-consistent locus evidence.

Do not ask RAD to solve all five jobs.

## 11. Claim ceilings

- 161 NMNS records are an authority candidate universe, not 161 independent species or transitions.
- 38 Moreyra concepts are a hypothesis-origin subset, not the Chapter 3 sampling frame.
- repeated use of JPN36/JPN06/JPN15 across P01–P04 is allowed as data reuse, but cannot be counted as multiple independent transition neighbourhoods.
- new non-Moreyra taxa are not called close relatives before nuclear placement.
- four P03 neighbourhoods are candidate transition replicates, not automatically four independent evolutionary origins.
- no RAD-only adaptation or pollinator-agent claim.

## 12. Current state

- physical biological samples: **0**;
- authority-wide A0 snapshot frozen: **false**;
- A1 field collection authorized: **false**;
- focal RAD pilot authorized: **false**;
- A2 target-capture augmentation authorized: **false**;
- A3 transition-neighbourhood RAD authorized: **false**;
- M01 expansion authorized: **false**;
- selection / pollinator-agent claims authorized: **false**.
