# aza3 reactivation-history gate v1

Status: required before any focal event is called evolutionary reactivation or regain.

## Purpose

A molecular result can look like pathway reactivation even when the historical polarity is wrong.

Therefore latent-pathway biology and trait-history inference are separated.

## Historical classes

### H_REGAIN_SUPPORTED

Use only when a loss -> regain history is supported across the admitted topology / branch-length / transition-model uncertainty.

Minimum requirements:
- ancestral state before the loss is identified under the frozen state model;
- the loss interval is identifiable enough to distinguish loss from ancestral polymorphism/sorting;
- the descendant recurrent state is not better explained by introgression from a lineage that retained the trait;
- regain support remains under the predeclared topology/rate sensitivity set.

### H_REGAIN_COMPATIBLE

Some admissible histories contain loss -> regain, but equally admissible histories permit no-regain explanations.

Allowed language:
- regain-compatible;
- candidate re-expression;
- molecular evidence consistent with reactivation.

Forbidden language:
- demonstrated reversal;
- restored ancestral pathway;
- re-evolved state.

### H_SUPPRESSION_ONLY

A derived phenotype shows regulatory/pathway suppression and retained machinery, but there is no identified later regain event.

This can still test whether phenotypic loss preserves future evolvability.

### H_RETICULATE_REACQUISITION

The apparent regain is better explained by introgression / hemiplasy from a lineage retaining the state.

This is not within-lineage reactivation, even if the imported allele turns the same pathway back on.

### H_NO_REGAIN

The best-supported history contains no loss -> regain sequence.

Do not use reversal framing.

### H_NOT_IDENTIFIABLE

Ancestral state, topology, transition polarity or taxonomic coding prevents classification.

## Evidence order

Historical classification is frozen before opening focal molecular reactivation results.

Order:
1. freeze tips/population states;
2. freeze topology ensemble and branch-length treatment;
3. freeze transition model / root-state prior or sensitivity grid;
4. classify historical event;
5. only then open pathway integrity / expression / regulatory evidence.

## Population-aware rule

Species-tip compression can hide recent colour reversibility.

When a taxon contains verified white and coloured populations or individuals, population/sample-aware states are required for the primary recent-history analysis.

Species-level ambiguity may be retained as a sensitivity analysis but cannot replace observed morph-linked states.

## Reticulation rule

Before calling reactivation, test whether the relevant haplotype/local genealogy is inherited from the focal lineage or reacquired from another lineage.

If introgression provides the intact program, classify the mechanism as reticulate reacquisition rather than spontaneous reactivation of a suppressed program inside the focal lineage.

## Relationship to molecular evidence

Historical and molecular evidence cross-classify.

Examples:

H_REGAIN_SUPPORTED + retained structural genes + loss-state suppression + restored homologous expression
-> REACTIVATION_SUPPORTIVE / potentially REACTIVATION_STRONG.

H_REGAIN_COMPATIBLE + same molecular pattern
-> pathway reactivation is molecularly plausible, but historical regain remains unresolved.

H_NO_REGAIN + suppression
-> regulatory suppression / latent capacity may still be real, but not evolutionary reversal.

H_RETICULATE_REACQUISITION + restored expression
-> imported pathway/haplotype reactivation, not within-lineage reversal.

## Current Cirsium consequence

The existing East Asian colour systems are candidates, not yet universal reversal examples.

EAzami already shows that:
- the six morph-linked takaoense samples are compatible with regain on some topologies, but regain is not required across the full topology set;
- Arenicola loss-versus-regain polarity remains transition-rate / ancestral-state sensitive.

Therefore aza3 should first test pathway persistence/suppression without presupposing reversal, while the historical gate is refined.

## Nature consequence

The Nature paper does not require reversal.

If reversal is supported, it provides a powerful mechanistic example that latent ancestral programs can be reactivated.

If reversal is absent, the broader reconfigurable-integration result can still stand, and suppression/persistence can still explain why components remain evolutionarily available.