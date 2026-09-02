# aza3 — own-data discrimination after EAzami Chapter 2

This repository starts Chapter 3 from the uncertainty left by the completed public-data Chapter 2. It does not treat future RAD-seq, target-capture, transcriptomic or field data as retroactive confirmation.

## Locked starting point

EAzami Chapter 2 remains frozen at merge `62fa8c5c913c2b236e710f6bad366e80676aa78f`.

Within its admitted public topology ensemble:

- orientation requires at least four state changes;
- phyllary posture requires at least three;
- stickiness requires at least five;
- minimum counts are better resolved than individual event placements;
- species-tip compression hides within-species state multiplicity.

## One question, breadth and depth

> Which trait histories still survive when species-tip states are returned to linked individuals, populations and nuclear ancestry, and how far can a retained transition be decomposed into mechanism without confusing history, selection and causation?

### P01-P05 breadth

The frozen core asks where apparent transitions in orientation, phyllary posture, stickiness and related morphology remain identifiable after adding same-individual phenotype, population replication, cytotype and nuclear ancestry.

### M01 depth — embedded worked case, not P06

The white versus bluish-purple floral-pigmentation contrast in `Cirsium brevicaule` and `C. irumtiense` is one **embedded worked case** of the same identifiability problem. M01 follows one retained contrast through ancestral/transition history, pigment chemistry, same-individual corolla expression and population-genomic evidence consistent with selection after background controls.

Loss, regain/re-evolution and selective-agent attribution are output-dependent claims, not assumptions.

## Transition-first sampling v3 — authoritative new-sampling plan

Physical biological samples currently equal **0**. New sampling is governed by `docs/CHAPTER3_TRANSITION_FIRST_SAMPLING_PLAN_V3.md` and `data/planning/chapter3_transition_first_sampling_v3.json`.

### Japan38 = hypothesis-origin subset, not sampling universe

Moreyra et al. sampled 38 Japanese species, and Chapter 2 P01-P05 were discovered on that subset. The current NMNS `日本のアザミ` index returns **161 authority records**; those records include infraspecific/duplicate taxonomic units and are not treated as 161 independent species or transitions.

Therefore Chapter 3 no longer starts from `38 concepts × n individuals`.

The acquisition order is:

> authority-wide digital trait screen → nuclear placement where needed → transition-neighbourhood selection → population RAD

A reproducible NMNS candidate-universe builder is provided in `analysis/build_nmns_transition_candidate_universe_v3.py`. Authority catchphrases are converted only to short screening states plus source hashes; authority states are not same-individual observations or ancestral states.

## Two kinds of overlap

**Same-individual multitrait overlap is desirable.** Every RAD individual scores orientation, phyllary posture, stickiness, colour and cytotype so P04 can test shared history without unlinking phenotype from ancestry.

**Evolutionary-neighbourhood overlap is not replication.** P03 requires at least **four non-overlapping candidate orientation-transition neighbourhoods**. Two populations of one species do not count as two transition replicates, and JPN36 can occupy at most one P03 neighbourhood slot.

## Fixed focal population layer

Default focal population bank = **12 primary RAD +3 preregistered reserves =15 physical plants/population**.

Minimum defensible launch:

- JPN36 `C. sieboldii`: 3 populations =45 physical /36 RAD;
- JPN06 `C. dipsacolepis`: 2 populations =30 /24;
- JPN15 `C. lineare`: 2 populations =30 /24;
- `C. brevicaule`: Okinawa Honto + Amami Oshima =30 /24;
- `C. irumtiense`: Miyako + Ishigaki =30 /24.

**Minimum =165 physical /132 primary RAD.**

Target mainland replication is JPN36=4 populations and JPN06/JPN15=3 populations each while M01 remains at two discovery populations/lineage. This gives **210 physical /168 primary RAD**.

For M01, floral RNA/pigment sampling is nested within the same discovery plants: six floral collections/population, five primary RNA-seq libraries and one reserve.

## RAD assay pilot is nested inside the five mandatory systems

Do not collect unrelated taxa solely to test enzymes.

- Stage A-focal: 5 systems ×1 DNA ×3 enzyme candidates = **15 shallow libraries**;
- Stage B-focal: 5 systems ×2 independent biological templates ×top-2 protocols =20 primary +5 independent-batch repeats = **25 libraries**;
- technical gates remain genotype concordance ≥0.95 and core-locus recovery ≥0.90.

This qualifies only the focal stratum. If later transition neighbourhoods introduce tetraploid/high-ploidy material, a second stratum-specific assay check is required before production RAD.

## Non-Moreyra species enter through placement, not assumption

A taxon outside the current Comp1061 backbone can be scientifically valuable, but geography or similar morphology does not make it a close relative.

The authority screen may nominate an initial batch of **up to 12 non-backbone concepts**. For each, collect two placement representatives where feasible and use a Comp1061-compatible target-capture nuclear placement assay. Only after stable placement may that concept enter population RAD as a transition comparator.

Ambiguous/hybrid/polyploid placement remains explicit; it is not forced into a sister pair.

## Transition-neighbourhood sampling

P01 keeps JPN36 as the fixed phyllary focal system. P02 keeps JPN06-JPN15 as the fixed stickiness contrast.

P03 does not sequence 38 or 161 units uniformly. It selects at least four **distinct orientation-transition neighbourhoods** from the admitted expanded nuclear scaffold. Each neighbourhood must contain an explicit state contrast and the placed comparator lineages needed to distinguish terminal from deeper placement.

Chapter 2 did not assign posterior weights to equally parsimonious histories, so v3 does **not** use the old `80% transition probability mass` language.

For an added non-focal comparator taxon, the default population design is:

- 2 populations;
- 10 primary RAD +2 reserve/population;
- 20 primary RAD /24 physical per taxon.

For within-species polymorphism, prefer at least three populations and at least 10 primary individuals per observed state across at least two populations where feasible, so morph is not identical to locality.

P04 reuses the same A1/A3 plants. P05 is measurement-first and does not open extra RAD merely to increase continuous-trait n.

## Method roles

- Moreyra/expanded Comp1061-compatible target capture: cross-species placement and scaffold;
- RAD: population ancestry, admixture, reticulation and morph/population history inside qualified neighbourhoods;
- same-individual phenotype: transition localization;
- M01 RNA/pigment: molecular decomposition;
- candidate resequencing/capture/WGS: selection-consistent locus evidence.

RAD is therefore an instrument for shallow ancestry/history discrimination, not the definitive Japanese species tree or an adaptation test.

## Individual field package

Every admitted RAD plant links one immutable ID to taxon/population, voucher/diagnostic images, orientation, phyllary posture and calibrated image, stickiness/gland state, flower colour, developmental stage, silica leaf, cytotype/genome-size evidence status and deidentified authorization/conservation IDs.

Spread individuals across the occupied patch; use ≥10 m spacing where biologically sensible and avoid obvious connected/clonal ramets. Reserves replace technical failures, not inconvenient biological outcomes.

## M01 expansion remains conditional

Do not add Iriomote, Yonaguni and extra Amami populations merely for island coverage. Open that expansion only after ancestral-state/pigment history is interpretable at the stated claim ceiling and the corolla-expression association replicates across both discovery populations within each lineage.

RAD supplies the population/background layer; final candidate-region selection evidence requires an independent non-RAD confirmation route. Genomics alone cannot identify the selective agent.

## Start here

1. `docs/CHAPTER3_TRANSITION_FIRST_SAMPLING_PLAN_V3.md` — authoritative transition-first sampling design.
2. `data/planning/chapter3_transition_first_sampling_v3.json` — machine-readable v3 contract.
3. `analysis/build_nmns_transition_candidate_universe_v3.py` — authority-wide candidate-universe builder.
4. `docs/CHAPTER3_SCOPE_AND_HANDOFF_V1.md` — frozen Chapter 3 question and claim boundaries.
5. `data/planning/chapter3_sampling_priorities_v1.csv` — frozen P01-P05 discovered from Chapter 2.
6. `data/planning/chapter3_mechanistic_demonstration_v1.json` — M01 evidence ladder.
7. `docs/CHAPTER3_RADSEQ_END_TO_END_V1.md` — general RAD laboratory/analysis safeguards where compatible with v3.
8. `data/intake/chapter3_individual_intake_v1.csv` — empty biological intake ledger.
9. `data/intake/chapter3_radseq_library_intake_v1.csv` — empty RAD library/QC ledger.

v1/v2 inverse-sampling, 38×3, eight-anchor and physical-inventory documents remain visible as design history but do not authorize or constrain new field acquisition.

## Current state

- own biological data admitted: **0**;
- physical samples: **0**;
- A0 authority snapshot frozen: **false**;
- A1 field sampling authorized: **false**;
- focal RAD pilot authorized: **false**;
- A2 target-capture augmentation authorized: **false**;
- A3 transition-neighbourhood RAD authorized: **false**;
- M01 expansion authorized: **false**;
- regain claim authorized: **false**;
- selection claim authorized: **false**;
- pollinator-agent claim authorized: **false**.

Run the Chapter 3 validators before accepting changes to these boundaries.
