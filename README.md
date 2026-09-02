# aza3 — own-data discrimination after EAzami Chapter 2

This repository starts Chapter 3 from the uncertainty left by the completed public-data Chapter 2. It does not treat future RAD-seq, genomic, transcriptomic or field data as retroactive confirmation.

## Locked starting point

EAzami Chapter 2 remains frozen at merge `62fa8c5c913c2b236e710f6bad366e80676aa78f`.

Within its admitted public topology ensemble:

- orientation requires at least four state changes;
- phyllary posture requires at least three;
- stickiness requires at least five;
- minimum counts are better resolved than individual event placements;
- species-tip compression hides within-species state multiplicity.

## One Chapter 3 question, breadth to depth

> Which trait histories still survive when species-tip states are returned to linked individuals, populations and nuclear ancestry, and how far can a retained transition be decomposed into mechanism without confusing history, selection and causation?

### P01-P05 breadth

The frozen core asks where apparent transitions in orientation, phyllary posture, stickiness and related morphology remain identifiable after adding same-individual phenotype, population replication, cytotype and nuclear ancestry.

### M01 depth, embedded case not P06

The white versus bluish-purple floral-pigmentation contrast in `Cirsium brevicaule` and `C. irumtiense` is used as one worked example of the same identifiability problem. M01 follows one retained contrast through:

1. ancestral-state / transition history;
2. pigment chemistry;
3. same-individual corolla expression;
4. population-genomic evidence consistent with selection after background controls.

Loss, regain/re-evolution and selective-agent attribution are output-dependent claims, not assumptions.

## Zero physical samples: authoritative sampling plan v2

The current sampling design starts from **zero physical biological material**. The authoritative field/genomic plan is `docs/CHAPTER3_ZERO_BASELINE_SAMPLING_PLAN_V2.md` with machine-readable contract `data/planning/chapter3_zero_baseline_sampling_v2.json`.

The earlier `38 x 3` design is **not the first sampling step**. It is retained only as historical design context / an optional late all-Japan RAD sensitivity.

### S0 — direct trait and occurrence reconnaissance

Before broad RAD expansion, repair individual/population phenotype coverage across Japan38. Target two populations x five flowering plants per concept where feasible; use one population x >=5 for narrow endemics. Destructive collection is not required. Bank up to three silica leaves per concept only when easy and authorized.

S0 freezes the direct trait/polymorphism table used to nominate P03/P04 transition neighbourhoods.

### S1 — minimum focal launch

Default population bank = **15 physical plants = 12 primary RAD +3 preregistered reserves**.

Minimum S1:

- JPN36 `C. sieboldii`: 3 populations = 45 physical /36 initial RAD;
- JPN06 `C. dipsacolepis`: 2 populations = 30 /24;
- JPN15 `C. lineare`: 2 populations = 30 /24;
- `C. brevicaule`: Okinawa Honto + Amami Oshima = 30 /24;
- `C. irumtiense`: Miyako + Ishigaki = 30 /24.

**S1 minimum = 165 physical plants /132 initial RAD individuals.**

This number is a consequence of population replication, not a study-level quota.

For M01 discovery, collect developmentally matched floral material from six plants/population, sequence five primary RNA replicates and retain one molecular reserve. Record sexual morph in `C. irumtiense`.

### RAD assay pilot is nested inside S1

Do not collect biologically irrelevant taxa merely to test enzymes.

The first focal-stratum pilot uses the five mandatory S1 systems: JPN36, JPN06, JPN15, `C. brevicaule`, `C. irumtiense`.

- Stage A-focal: 5 systems x1 qualified DNA x3 enzyme candidates = **15 shallow libraries**;
- Stage B-focal: 5 systems x2 independent biological templates xTOP2 protocols =20 primary +5 independent-batch technical repeats = **25 libraries**;
- technical gates remain genotype concordance >=0.95 and core-locus recovery >=0.90.

This qualifies only the focal production stratum. If later S2/S4 introduces polyploid/high-ploidy concepts, run a second stratum-specific protocol check before those samples are sequenced.

### S2 — P03/P04 transition-neighbourhood RAD

After S0 is frozen, use the fixed target-capture topology ensemble and direct trait states to select the smallest concept neighbourhood set covering **>=80% cumulative transition-placement uncertainty** for each trait, plus every directly observed polymorphic concept. Take the union across orientation, phyllary and stickiness.

Default S2 per concept = two populations, six primary RAD +two reserves/population = 12 primary +4 reserves/concept. The number of S2 concepts is an output of the frozen uncertainty-cover rule, not a quota.

### S3 — conditional M01 expansion

Only after M01 E1/E2 succeeds, add:

- `C. brevicaule`: one intermediate Amami-group + one southern Amami-group population;
- `C. irumtiense`: Iriomote + Yonaguni.

S3 adds 60 physical /48 initial RAD. Final M01 bank = eight populations x15 =120 physical, with 96 initial RAD and 24 reserves.

RAD supplies population/background structure. Candidate-region selection evidence must use a separate non-RAD confirmation assay.

### S4 — all-Japan RAD sensitivity, optional and last

Only if still needed after S1/S2 and only if shared-locus/ploidy gates pass, use the S0 silica bank to build a Japan-wide RAD sensitivity. This product is not the primary species tree.

The Moreyra 2025 target-capture framework remains the cross-species nuclear scaffold.

## Individual field package

Every admitted RAD plant links one immutable ID to:

- taxon and population;
- voucher or diagnostic images;
- orientation;
- phyllary posture / calibrated image;
- stickiness / gland-exudate state;
- flower colour;
- developmental stage;
- silica leaf DNA;
- cytotype/genome-size evidence status;
- deidentified authorization/conservation IDs.

Spread sampling across the occupied patch; use >=10 m spacing where feasible and avoid obvious connected/clonal ramets. Reserves may replace technical failures, not inconvenient biological outcomes.

## RAD analysis boundary

RAD is an instrument for population ancestry, structure and network sensitivity, not an unconditional species-tree or adaptation instrument.

- primary assembly: Stacks 2 de novo within qualified strata;
- m/M/n and missingness/locus filters are tuned on technical metrics before trait outcomes;
- mixed ploidies are not forced through one diploid caller;
- RAD-only two-species FST or one outlier scan is not final M01 selection evidence.

## Start here

1. `docs/CHAPTER3_ZERO_BASELINE_SAMPLING_PLAN_V2.md` — authoritative zero-material field/RAD plan.
2. `data/planning/chapter3_zero_baseline_sampling_v2.json` — machine-readable v2 contract.
3. `data/planning/chapter3_zero_baseline_population_targets_v2.csv` — S1/S3 population targets.
4. `docs/CHAPTER3_SCOPE_AND_HANDOFF_V1.md` — Chapter 3 question and claim boundaries.
5. `data/planning/chapter3_sampling_priorities_v1.csv` — frozen P01-P05.
6. `data/planning/chapter3_mechanistic_demonstration_v1.json` — M01 evidence ladder.
7. `docs/CHAPTER3_RADSEQ_END_TO_END_V1.md` — general RAD laboratory/analysis safeguards retained where compatible with v2.
8. `data/intake/chapter3_individual_intake_v1.csv` — empty biological intake ledger.
9. `data/intake/chapter3_radseq_library_intake_v1.csv` — empty RAD library/QC ledger.

Earlier inverse-sampling, eight-anchor pilot and physical-inventory documents remain visible as design history but are superseded for new field acquisition by zero-baseline v2.

## Current state

- own biological data admitted: **0**;
- confirmed physical samples: **0**;
- S0 started: **false**;
- S1 sampling authorized: **false**;
- focal RAD pilot authorized: **false**;
- S2 open: **false**;
- M01 S3 open: **false**;
- all-Japan S4 open: **false**;
- regain claim authorized: **false**;
- selection claim authorized: **false**;
- pollinator-agent claim authorized: **false**.

Run the Chapter 3 validators before accepting changes to these boundaries.
