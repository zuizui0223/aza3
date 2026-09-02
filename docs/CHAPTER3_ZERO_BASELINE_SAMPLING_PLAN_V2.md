# Chapter 3 zero-baseline sampling plan v2

Status date: 2026-09-02

## Why v2 exists

This plan assumes **zero physical biological material at the start**. It supersedes the earlier habit of beginning from a fixed `38 concepts x 3 individuals` RAD quota or from a physical-inventory audit. The sampling unit is chosen from the inference first.

The cross-species Japanese Cirsium scaffold is not assigned to RAD-seq. The public Moreyra et al. 2025 target-capture framework supplies the primary species-level nuclear scaffold. RAD-seq is reserved for the shallower information gaps that require multiple individuals: population ancestry, admixture/reticulation, within-concept trait polymorphism and localization of selected trait-history alternatives.

The plan therefore has four biological stages plus one optional breadth product.

## Universal individual package

Whenever a flowering plant is admitted as a RAD biological individual, assign one immutable `individual_id` before any molecular result and record:

- taxon determination and population ID;
- private exact locality / public deidentified locality key;
- voucher ID or voucher-linked diagnostic images;
- standardized whole-plant, capitulum lateral and involucre close-up images;
- orientation state;
- phyllary posture and calibrated phyllary images;
- stickiness / gland-exudate state;
- flower-colour record;
- phenological/developmental stage;
- silica-dried leaf tissue for DNA;
- cytotype/genome-size evidence status;
- authorization / conservation-review identifiers.

For `C. irumtiense`, record sexual morph (female versus hermaphroditic when diagnosable) because floral molecular measurements must not silently confound colour-lineage comparison with sex expression.

A population is a spatially distinct reproductive stand or habitat unit, not an arbitrary cluster of nearby plants. Exact separation criteria are site-specific and frozen before sampling. Within a stand, spread sampling across the occupied patch; use >=10 m spacing where feasible and avoid obviously connected/clonal ramets. Later kinship filtering may exclude close relatives but may not trigger outcome-driven replacement.

## S0 — trait and occurrence reconnaissance before genomic expansion

### Purpose

Repair the major pre-genomic uncertainty: direct individual trait coverage and within-concept polymorphism. Chapter 2 species-tip coding cannot represent polymorphic populations and many Japan38 concepts still lack complete direct capitulum states.

### Design

For each Japan38 concept that can be visited or documented under rights constraints:

- target: 2 distinct populations x 5 flowering plants photographed/scored;
- narrow endemic / single verified population: 1 population x >=5 flowering plants;
- no destructive collection is required for S0;
- where leaf collection is easy and authorized, bank silica tissue from up to 3 outcome-blind individuals per concept for possible later S4 use, but do not sequence it merely because it was banked.

S0 is distributed work: collaborators and rights-cleared voucher/herbarium evidence may repair trait states. The investigator is not expected to personally visit all 38 concepts.

### Output and decision rule

Freeze an updated individual/population trait table before selecting the P03/P04 RAD expansion. Using the fixed target-capture topology ensemble, calculate trait-history uncertainty and nominate transition-neighbourhood concepts by a predeclared coverage rule, not by convenience.

## S1 — launch focal RAD panel

S1 contains only systems that directly answer the locked high-priority history questions or the embedded M01 depth demonstration. All numbers below refer to **physical plants collected** and **initial RAD libraries** separately.

### Common per-population rule

Default field bank = **15 physical individuals/population**:

- 12 primary RAD individuals;
- 3 predeclared reserves.

The 3 reserves replace only extraction/library failures or preregistered close-kin exclusions. They do not replace biologically inconvenient genotypes.

For strongly conservation-limited populations, a reduced 10-plant field bank (8 primary +2 reserve) may be used only as a documented reduced panel; it cannot silently inherit the claim ceiling of a 12-primary population.

### S1A — JPN36 `Cirsium sieboldii`, P01 phyllary-history system

Minimum launch design:

- **3 geographically distinct populations x 15 physical = 45 plants**;
- initial RAD = **3 x 12 = 36 individuals**.

Target design if access is straightforward:

- 4 populations x 15 = 60 physical;
- 48 primary RAD.

Select populations to span the species' broad Honshu/Shikoku wetland distribution and, where known before genotyping, include populations expressing documented orientation/colour polymorphism. Do not construct a post-genotype morph-balanced sample.

### S1B — JPN06 `Cirsium dipsacolepis`, P02 stickiness-history system

Minimum:

- **2 populations x 15 = 30 physical**;
- 24 primary RAD.

Target:

- 3 populations x 15 = 45 physical;
- 36 primary RAD.

Use geographically separated verified extant grassland populations. Historical locality records are not sufficient without current confirmation.

### S1C — JPN15 `Cirsium lineare`, P02 comparator

Minimum:

- **2 populations x 15 = 30 physical**;
- 24 primary RAD.

Target:

- 3 populations x 15 = 45 physical;
- 36 primary RAD.

The three-population target should span separated Japanese regions when possible (for example western Honshu/Yamaguchi, Shikoku and Kyushu), but conservation/access constraints take precedence over geographic symmetry.

### S1D — M01 discovery, `C. brevicaule`

Discovery populations:

1. Okinawa Honto;
2. Amami Oshima.

Per population: 15 physical; 12 primary RAD. Discovery total = **30 physical / 24 initial RAD**.

Within each population collect developmentally matched floral material from 6 of the 15 plants for RNA/pigment work; sequence 5 biological RNA replicates and retain one molecular reserve.

### S1E — M01 discovery, `C. irumtiense`

Discovery populations:

1. Miyako;
2. Ishigaki.

Per population: 15 physical; 12 primary RAD. Discovery total = **30 physical / 24 initial RAD**. Record sexual morph on all focal flowering individuals.

### S1 totals

Minimum launch panel:

- JPN36 45 physical /36 RAD;
- JPN06 30 /24;
- JPN15 30 /24;
- `C. brevicaule` discovery 30 /24;
- `C. irumtiense` discovery 30 /24;
- **total = 165 physical plants, 132 initial RAD individuals**.

Target mainland replication (JPN36=4 populations, JPN06/JPN15=3 populations) while M01 remains discovery-only:

- **210 physical plants, 168 initial RAD individuals**.

These totals are consequences of the population design, not study-level quotas.

## RAD assay pilot is nested inside S1

Do not collect JPN03, JPN05 or JPN35 merely to test enzymes before the focal programme exists.

The first RAD assay pilot uses material from the five mandatory S1 systems:

- JPN36;
- JPN06;
- JPN15;
- `C. brevicaule`;
- `C. irumtiense`.

Stage A-focal: one qualified DNA template/system x three candidate enzyme protocols = **15 shallow libraries**.

Stage B-focal: two independent biological templates/system x TOP2 protocols = **20 primary libraries**, plus **5 independent-batch technical repeats**, one/system = **25 libraries**.

The existing 0.95 genotype-concordance and 0.90 core-locus-recovery technical gates remain. This pilot qualifies the focal diploid/Ryukyu production stratum only. It does **not** qualify high-ploidy/polyploid Japan-wide production.

If S2 later admits polyploid/high-ploidy concepts, run a second ploidy-stratum protocol check before sequencing those concepts. A technical requirement may not force biologically irrelevant field sampling into S1.

## S2 — P03/P04 transition-neighbourhood expansion

### Why not 38 x 3 first

Three individuals from every concept spend sequencing effort on lineages that may contribute little to transition localization while providing too little depth to characterize the lineages that matter most. The target-capture scaffold already supplies cross-species placement. S2 therefore allocates RAD to the concepts that cover the unresolved transition-placement set.

### Predeclared nomination rule

After S0 is frozen, for each of orientation, phyllary posture and stickiness:

1. reconstruct the admissible transition-edge set over the fixed public topology ensemble;
2. assign each concept/edge a contribution to transition-placement uncertainty;
3. include every directly observed polymorphic concept;
4. then choose the smallest set of additional concepts whose sampled neighbourhoods cover >=80% of cumulative transition-placement uncertainty for that trait;
5. take the union across the three modules before opening S2.

Do not substitute a convenient taxon after seeing its RAD topology.

### Per-concept S2 design

Default localization panel:

- 2 populations/concept where distribution allows;
- **6 primary RAD +2 reserve per population**;
- 12 primary /4 reserve per concept.

For single-population/narrow endemic concepts, use the maximum authorized non-destructive/tissue sample and lower the claim ceiling rather than inventing a second population.

Expected S2 size if 8 concepts are nominated: **96 primary RAD +32 reserve physical samples**. If 12 concepts are nominated: **144 primary +48 reserve**. The number of nominated concepts is an output of the frozen uncertainty-cover algorithm, not a target chosen in advance.

P04 reuses the same S1/S2 individuals; orientation, phyllary and stickiness must be scored before genotypes are inspected.

## S3 — M01 population expansion only after E1/E2 promotion

Do not collect or sequence the full southern/central Ryukyu panel solely because it is geographically attractive.

Open S3 only if ancestral-state/pigment analysis is interpretable and the corolla expression mechanism replicates across both discovery populations within each lineage.

Add:

`C. brevicaule`:

- one intermediate Amami-group population;
- one southern Amami-group population.

`C. irumtiense`:

- Iriomote;
- Yonaguni.

Each added population: 15 physical, 12 primary RAD +3 reserve.

S3 adds **60 physical plants /48 initial RAD individuals**. Final M01 bank becomes 8 populations x15 = **120 physical**, with 8 x12 = **96 initial RAD** and 24 predeclared reserves.

RAD supplies population structure and demographic/background information. If E2 nominates a pigmentation candidate, candidate-region selection evidence is validated with a non-RAD assay (targeted capture/amplicon/qualified low-pass or WGS). The RAD sample bank can be reused for that assay, and reserves can be sequenced prospectively if power calculations justify it.

## S4 — optional all-Japan RAD sensitivity, last not first

S4 is opened only if the dissertation still requires a single-protocol Japan-wide RAD sensitivity after S1/S2 and the technical shared-locus gate passes.

Use the S0 silica bank where available. For concepts not already represented by S1/S2, admit up to three wild individuals split across >=2 populations when feasible. This is a **secondary topology/network sensitivity** only, not a definitive species tree or population-genomic panel.

If restriction-site dropout or mixed ploidy prevents a common matrix, cancel S4 and retain the target-capture scaffold plus focal RAD population layers.

## Cytotype/genome-size sampling

Do not require fresh-flow material from every RAD plant if logistics make that impossible.

For a presumed uniform diploid focal population, target direct flow-cytometry/genome-size confirmation on **3 individuals/population** when fresh transport is feasible. Increase this if within-population cytotype variation is suspected.

For polyploid/high-ploidy S2/S4 concepts, cytotype qualification is mandatory before dosage-dependent pooled analysis; unresolved individuals remain usable for phenotype/voucher work but not for a falsely diploid genotype matrix.

Genome size and chromosome/ploidy are separate variables.

## Personal versus distributed field work

The design does not require one investigator to collect every plant.

Priority for direct investigator effort:

1. the four M01 discovery populations because RNA/pigment developmental matching is sensitive to protocol;
2. JPN36 because it links P01 and known orientation/colour polymorphism;
3. JPN06/JPN15 populations that require careful identity/conservation handling.

S0 breadth reconnaissance and many S2/S4 leaf/photo packages should be distributed to trained collaborators under one field sheet and ID convention.

## Earliest practical seasonal sequence

- JPN36 flowers mainly August-October; JPN06 mainly September-November; Japanese `C. lineare` is an autumn-flowering system (roughly September-October in published Japanese cultivation/field references). These make an autumn mainland collection block.
- `C. brevicaule` is reported mainly December-April in the Ryukyus; M01 discovery should use a local pre-survey rather than a fixed date because flowering can extend beyond the core period.
- `C. irumtiense` exact site-level phenology should be confirmed before travel; northern/southern island populations must be collected at matched capitulum developmental stages for RNA.

## Claim ceilings by stage

- S0 only: direct trait distributions/polymorphism; no own genomic history.
- S1: focal population ancestry and P01/P02/M01 discovery history; no Japan-wide RAD species tree.
- S2: transition-neighbourhood history localization for P03/P04.
- S3: M01 multi-island population background; RAD alone still not a selection result.
- S4: optional all-Japan RAD sensitivity only.

The recommended first commitment from a zero-material start is therefore **S0 reconnaissance + S1 minimum launch panel**, not 38 x3 and not a full eight-island M01 population-genomic project.