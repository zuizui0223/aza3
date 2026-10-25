# Cirsium sieboldii Gate CS0 field card v1

**Purpose:** verify whether a natural population contains a directly observed combinatorial substrate before biological WGS sample selection.

**This is observation-only.** No tissue collection, tagging, off-path entry or manipulation is authorized by this card.

## Taxonomic rule

Use `docs/C_SIEBOLDII_CS0_TAXONOMIC_DIAGNOSTIC_V1.md`.

Primary Gate CS occupancy uses **high-confidence focal C. sieboldii individuals only**. Medium-confidence records are sensitivity-only; unresolved records are descriptive only.

Never use upright or nodding orientation itself to decide that the plant is focal *C. sieboldii*.

## Before counting

1. Confirm *C. sieboldii* identity.
2. Confirm flowering heads are available at A1/A2.
3. Walk the permitted observation route once before scoring.
4. Define the census direction/route before choosing individuals.
5. Do not start from conspicuous white or upright plants.

## Individual census

Target **30 flowering individuals** if visible and distinguishable.

Primary classification requires **>=20 complete A1/A2 individuals** with both colour and orientation recorded.

For every individual:

- immutable observation ID;
- population/site ID;
- flower stage: A0 / A1 / A2 / A3 / P;
- standardized whole-head photo;
- diagnostic involucre photo;
- colour calibration image where feasible;
- continuous colour/chroma value when derivable;
- descriptive white/non-white state;
- gravity-referenced orientation angle:
  - 0° vertically upward;
  - 90° horizontal;
  - 180° vertically downward;
- phyllary rows;
- outer-phyllary posture/angle if measurable;
- stickiness only when a separately authorized non-destructive observation protocol permits it;
- habitat microzone;
- taxonomic confidence;
- visible reason for any missing value.

## Primary C x O gate

Only **A1-A2** individuals enter the primary gate.

Orientation substrate:
- observed angle range >=30°;
- >=3 individuals <=60°;
- >=3 individuals >=120°.

Colour substrate:
- white/non-white may be used for the transparent 2x2 display only when each state has >=3 individuals;
- continuous colour remains the stored primary variable.

Combinatorial substrate:
- >=3 of four white/non-white x upright/nodding cells occupied;
- every occupied counted cell >=3 individuals;
- |Spearman(colour, orientation)| <0.8;
- no pattern created by flower stage.

## Immediate field interpretation

**Potential CS_GREEN:** one population visibly satisfies all primary thresholds.  
Do not declare final GREEN until the CSV is run through the frozen classifier.

**Potential CS_GREEN_LOCAL_SET:** no population passes alone, but a predeclared nearby local set appears to contain the combinations.  
Do not pool sites post hoc. Each colour state and each orientation extreme must occur in >=2 component populations.

**NOT_IDENTIFIABLE:** too few A1/A2 individuals, insufficient visibility, or access prevents an unbiased census.

## Site-specific priorities

### Kurumayama Moor / Kirigamine

Primary goal: verify white-form frequency and whether ordinary coloured plants are present in the same observable population.

Secondary goal: record A1/A2 orientation on **all** scored white and non-white individuals.

Boundary: white flowers at this locality do not themselves establish combinatorial reassembly.

### Kiyooka-Mukaiyama Wetland

Primary goal: standardized A1/A2 orientation census across the boardwalk-observable population.

Secondary goal: prospectively score colour for every censused individual, including absence of white individuals.

Boundary: do not preferentially select unusually upright heads.

### Kashibaru Wetland

Role: orientation-validation site only.

Primary goal:
- test whether the published predominance of upward-facing **flowering** Maa-azami can be reproduced;
- verify focal *C. sieboldii* identity on every scored individual using the diagnostic characters available from permitted observation positions.

Mandatory downgrade:
- if *C. sieboldii* cannot be separated confidently from Satsuma-maa-azami on the scored individuals, classify the site as taxonomically NOT_IDENTIFIABLE for Gate CS;
- do not pool uncertain individuals with Kurumayama or Tsukude observations.

This site cannot become the primary colour × orientation substrate merely because upright flowering is common.

## Hard stops

- P/post-anthesis upright heads never count as upright-flowering morphs.
- Do not leave public access infrastructure where entry is restricted.
- Do not combine distant populations to manufacture 3/4 cells.
- Do not infer absence when fewer than 20 complete A1/A2 individuals are observable.
- Do not choose Nature WGS individuals until Gate CS is formally classified.
