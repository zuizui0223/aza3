# Chapter 3 field priority and calendar v5

Status date: 2026-09-02

This document operationalizes the authoritative nationwide species-tree design v4. It does not replace the Level 1–4 architecture. Its purpose is to decide **what to collect first, in which season, and how to avoid sacrificing nationwide taxon coverage to premature population depth**.

## Governing rule

The primary field objective is the Japan-wide nuclear species-tree/network.

For every field block:

1. **tree breadth first** — recover the Level-1 target-capture representatives for as many admitted species as possible;
2. **focal top-up second** — only after the tree representatives are secured, extend focal species to population-RAD sample sizes;
3. **enrichment last** — third/fourth individuals for widespread, variable, polyploid, hybrid-suspected or taxonomically difficult species are collected when they do not displace missing primary species.

The operational census must be frozen before any claim that the final Japanese sampling universe equals 120 species. The current planning baseline remains approximately 120 species; 161 NMNS database rows are not treated as 161 independent species.

## Sample package by tier

### T1 — nationwide tree core

Default: **2 independent biological individuals/species**.

Preferred design:

- two different populations or separated geographic occurrences;
- if only one verified population exists, two clearly separate plants may be used and the species is flagged geographically limited;
- rare/conservation-limited taxa may remain `SINGLE_SAMPLE_LIMITED` rather than being substituted by another taxon.

Each individual links:

- immutable individual ID;
- field taxon determination and later accepted operational concept;
- population/locality key;
- voucher or diagnostic image set;
- head orientation;
- phyllary posture and calibrated image;
- stickiness/gland-exudate state;
- flower colour;
- developmental stage/phenology;
- silica-dried leaf for DNA;
- fresh tissue for cytotype/genome-size work when logistically feasible;
- authorization/conservation identifiers kept deidentified in the public repository.

### T2 — widespread / visibly variable species

Default: **3 target-capture individuals/species**, preferably from three geographic occurrences or from two populations spanning the observed state variation.

### T3 — complex/polyploid/hybrid-suspected species

Default: **4 target-capture individuals/species**, preferably 2 populations ×2 individuals, with cytotype/genome-size qualification where feasible. Individuals are retained through gene-tree/network analyses rather than collapsed a priori to one species consensus.

### F — focal population-RAD systems

Population bank: **12 primary RAD +3 reserves =15 physical plants/population**.

Full target:

- `Cirsium sieboldii`: 4 populations =60 physical /48 primary RAD;
- `C. dipsacolepis`: 3 populations =45 /36;
- `C. lineare`: 3 populations =45 /36;
- `C. brevicaule`: 2 discovery populations =30 /24;
- `C. irumtiense`: 2 discovery populations =30 /24.

Full focal bank = **210 physical plants /168 primary RAD**.

The first two nationwide-tree representatives for each focal species are chosen from geographically separated focal populations. They are not additional plants.

## Priority classes

### Priority A — seasonal rescue / hard-to-repeat species

Highest field priority because missing the flowering window delays the nationwide tree by a full season or year.

Includes:

- high-alpine species with July–September flowering;
- single-mountain or very narrow-range taxa;
- island endemics with narrow accessible populations;
- conservation-limited taxa for which collection requires a coordinated authorized visit.

Do not insist on extra enrichment individuals until at least one or preferably two valid tree representatives are secured.

### Priority B — focal hypothesis systems

Second priority because one collection supplies Level 1 tree, Level 2 traits and Level 3 population genomics.

- `C. sieboldii` — flower August–October; Honshu/Shikoku wetlands;
- `C. dipsacolepis` — flower September–November; dry grasslands;
- `C. lineare` — flower September–October; Yamaguchi/Shikoku/Kyushu grasslands;
- `C. brevicaule` — flower December–April, sometimes outside that window; Amami Oshima–Okinawa Honto;
- `C. irumtiense` — Miyako/Yaeyama; flowering individuals documented in late February–March; exact local peak must be verified before travel.

### Priority C — route-efficient nationwide backbone species

Common or moderately widespread admitted species encountered in the same geographic block. Collect their T1 two-individual package **before** spending extra time topping a focal species from 12 to 15 plants.

### Priority D — replication/enrichment completion

Second/third/fourth individuals, extra focal populations and difficult-network enrichment. These are scientifically important but should not displace a missing Level-1 species from the nationwide tree.

## Field calendar from 2026-09-02

### Campaign 1 — September 2026: northern and alpine seasonal rescue

**Priority A first.**

Target all operational-census taxa whose documented flowering closes in August–September and that cannot be recovered in the late-autumn mainland campaign.

Illustrative high-priority examples from the current NMNS resource include:

- `C. ugoense` — July–August, Tohoku alpine grassland;
- `C. chokaiense` — August–September, Mt Chokai endemic;
- `C. zawoense` — August–September, Zao range;
- `C. yezoalpinum` — August–September, Daisetsu/Shiretoko alpine zone;
- `C. okamotoi` — August–September, Mikuni mountains;
- `C. otayae` — August–October, Myoko/Northern Alps/Hakusan.

For every reachable taxon: secure the two tree individuals first. Focal RAD depth is not the objective of this campaign.

Remote Hokkaido/high-alpine taxa should preferentially be assigned to authorized collaborators when that avoids a dedicated single-species trip.

### Campaign 2 — late September to October 2026: central/mainland autumn + `C. sieboldii`

Objectives in order:

1. T1 representatives for all admitted flowering taxa encountered in the Honshu/Shikoku route;
2. `C. sieboldii` two-population RAD starter: 2 populations ×15 =30 physical /24 primary RAD;
3. collect additional T1 species along the same wetland/low-mountain blocks;
4. defer the third/fourth `C. sieboldii` populations until nationwide breadth is substantially filled or Level-2 remapping confirms the added value.

`C. sieboldii` flowers August–October, so it can follow the September alpine-rescue window.

### Campaign 3 — October to November 2026: western Honshu–Shikoku–Kyushu autumn block

Objectives in order:

1. T1 representatives for route-efficient western Japanese species;
2. `C. dipsacolepis`: 2-population RAD starter =30 physical /24 primary;
3. `C. lineare`: 2-population RAD starter =30 /24;
4. capture other autumn species in the same grassland/forest-edge blocks before collecting third focal populations.

`C. dipsacolepis` flowers September–November. `C. lineare` flowers September–October, so lineare sites must be scheduled earlier within this block.

### Campaign 4 — late February to late March 2027: Ryukyu discovery block

This is the preferred joint M01 window because `C. brevicaule` is within its documented December–April flowering season and flowering `C. irumtiense` is documented in late February and March.

Required discovery populations:

- `C. brevicaule`: Okinawa Honto + Amami Oshima;
- `C. irumtiense`: Miyako + Ishigaki.

Per population:

- 15 physical plants;
- 12 primary RAD +3 reserves;
- six developmentally matched floral collections for RNA/pigment, with five primary RNA libraries + one reserve;
- record sexual morph in `C. irumtiense`.

Also collect T1 two-individual representatives for every other admitted southern-island Cirsium encountered under the same permits and routes.

Do not add Iriomote/Yonaguni or extra Amami islands merely for geographic completeness before M01 E1/E2 promotion.

### Campaign 5 — July to September 2027: northern/high-alpine completion

Primary objective: fill every Priority-A species missed in September 2026 and obtain second geographic representatives where only one was recovered.

This block also handles taxa whose normal flowering is genuinely summer rather than autumn. The target is **zero missing high-alpine operational species at the end of the campaign**, subject to permits and conservation constraints.

### Campaign 6 — September to November 2027: nationwide mainland gap fill

Primary objectives:

1. complete two-per-species T1 representation;
2. obtain missing second-population representatives;
3. only then complete focal RAD to the full 4/3/3 mainland population design;
4. collect T2/T3 enrichment individuals;
5. close taxonomic or phenotype gaps identified by the first own target-capture tree/network.

## RAD depth is staged, not front-loaded

### RAD starter

Before the nationwide tree is complete, limit mainland focal population depth to:

- `C. sieboldii`: 2 populations ×12 primary =24;
- `C. dipsacolepis`: 2 ×12 =24;
- `C. lineare`: 2 ×12 =24;
- `C. brevicaule`: 2 ×12 =24;
- `C. irumtiense`: 2 ×12 =24.

Starter total = **120 primary RAD** in 10 populations, banked as 150 physical plants including reserves.

### RAD completion

After Level-1 tree admission and Level-2 trait-history remapping, add only if still justified:

- `C. sieboldii`: +2 populations =+24 primary;
- `C. dipsacolepis`: +1 =+12;
- `C. lineare`: +1 =+12.

Completion adds **48 primary RAD**, reaching the v4 target of 168.

This prevents 48 extra population-genomic samples from displacing missing species-tree taxa.

## Target-capture sequencing waves

### Technical qualification batch

As soon as a representative spread of material is available, use **12 species ×2 individuals =24 samples** spanning geography, DNA quality and cytotype/genome-size classes to qualify extraction/capture/QC with the chosen Comp1061-compatible assay. This is a technical check, not a separate biological sample quota; all 24 are intended to remain in the nationwide tree.

### Wave 1 — breadth

Sequence the two-per-species panel, nominally **240 individuals at 120 species**. If field sampling is incomplete, prioritize taxon breadth over third/fourth individuals.

### Wave 2 — complexity

After the initial tree/network is inspected under frozen QC, sequence the preregistered third/fourth individuals for widespread/variable and complex/polyploid/hybrid-suspected taxa, nominally up to **50 additional individuals**, giving the v4 full panel of 290.

## Trip-level stopping rule

On any trip, if time or collection limits force a choice between:

- a missing admitted species for the nationwide tree; and
- an additional focal RAD individual beyond the first 12 primary +3 reserves in an already represented population,

**collect the missing species-tree taxon.**

The dissertation's first product is the nationwide species tree/network.