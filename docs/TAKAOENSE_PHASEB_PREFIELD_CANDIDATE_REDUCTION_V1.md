# takaoense Phase-B pre-field candidate-reduction plan v1

**Status:** active pre-field molecular screen  
**Date:** 2026-10-05

## Why this exists

*Cirsium japonicum* var. *takaoense* is the strongest current molecular anchor for the adaptive-variant programme because six public transcriptome individuals have recovered colour labels:

| code | run | morph |
|---|---|---|
| FC | SRR35152718 | BP |
| WY | SRR35152717 | W |
| FB | SRR35152738 | W |
| TJ | SRR35152736 | BP |
| NH | SRR35152735 | BP |
| LT | SRR35152734 | W |

These are **one plant per locality** and the RNA was not collected as a matched floral-expression experiment.

Therefore this panel can reduce hypotheses before new field sampling, but it cannot establish colour causation or adaptation.

## Primary pre-field question

> Do the six public morph-labelled transcriptomes contain shared coding/haplotype differences that nominate a tractable molecular route, or does the public panel leave regulatory/ancestry explanations unresolved?

## What the public panel may test

### T1 — sample genealogy / ancestry context

Reanalyse the six W/BP-labelled samples using the published/public orthologous transcriptome information where available.

Ask:

- do W samples form one supported cluster?
- do BP samples form one supported cluster?
- is morph grouping stable across gene/topology uncertainty?
- do candidate colour genes have local genealogies discordant with the genome-wide/sample background?

This is ancestry context, not a colour-gene association test.

### T2 — expressed coding-variant screen

Use the public reads only for genes that are actually covered sufficiently in the young-leaf transcriptomes.

Priority gene families:

- CHS;
- CHI;
- F3H;
- F3'H / F3'5'H where distinguishable;
- DFR;
- ANS/LDOX;
- UFGT/glycosyltransferase candidates where defensible;
- anthocyanin regulatory MYB/bHLH/WD40 candidates only when orthology/family assignment is defensible.

For each callable gene:

- recover expressed coding haplotypes;
- flag fixed W-versus-BP differences in the six-sample panel;
- flag nonsynonymous, truncating or splice-disrupting candidates;
- retain within-morph polymorphism;
- compare candidate-gene genealogy with the background sample genealogy.

A six-sample fixed difference is a **candidate discriminator only**, because morph and geography are confounded.

### T3 — pathway-retention screen

Existing DFR/ANS targeted searches already show homologous reads in both colour classes.

This supports pathway-retention plausibility.

The pre-field screen may extend retention checks to additional pathway components, but it must not turn gene presence into a floral-function claim.

## What is explicitly prohibited

- no differential-expression claim from these young-leaf libraries;
- no floral regulatory-restoration claim;
- no adaptive-gene claim;
- no causal-gene claim from a 3-vs-3 fixed difference;
- no treating morph clustering as proof that colour caused the genomic structure;
- no treating absence of a transcript as gene deletion;
- no post-hoc expansion to large unrelated gene families just because the first screen is negative.

## Decision classes

The pre-field result must end in one of:

### P1 — STRUCTURAL_CODING_CANDIDATE

A pathway/regulatory gene shows a strong morph-linked coding/haplotype candidate that survives coverage and orthology checks.

Next step:
population replication in the planned 20–30 W + 20–30 BP panel, followed by floral RNA/pigment and fitness.

### P2 — REGULATORY_ROUTE_FAVOURED

No convincing shared coding difference is found, while the pathway appears retained.

This does **not** prove regulatory causation.

Next step:
population-level dense genotype plus floral-stage expression and regulatory-region analysis.

### P3 — ANCESTRY_CONFOUNDED

Morph grouping is inseparable from genome-wide geography/ancestry in the six public samples.

Next step:
mixed or paired nearby W/BP populations become mandatory before any gene-level association.

### P4 — MULTIPLE_CANDIDATES

Several coding/regulatory candidates remain and cannot be ranked.

Next step:
freeze a population-replicated association/segregation design rather than selecting one preferred gene.

### P5 — NOT_IDENTIFIABLE

Coverage/orthology/sample structure does not permit a defensible candidate reduction.

Next step:
skip additional public-data fishing and proceed directly to the new population panel.

## How this connects to the AV evidence ladder

This public six-sample analysis cannot pass AV1.

At best it produces a **pre-AV1 candidate set**.

AV1 requires an ancestry-controlled trait association in new population-replicated data.

AV2 requires molecular mediation in floral tissue/pigment.

AV5 requires reproductive-fitness evidence.

AV6 requires recurrent reuse across independent histories.

## Stop rule

The public takaoense panel is allowed to reduce the candidate search space once.

After this screen, no additional untargeted SRA/BLAST/transcriptome fishing is opened unless it resolves a prespecified candidate decision.
