# P04 — source-verified Cirsium evidence readiness for evolutionary afterlife

**Audit: 19 October 2026.** Completed source-only provenance audit. No new plant collection, floral experiment, recovered trait, or fitness result.

## Present source-derived conclusion

Can current published Cirsium materials identify a *white floral lineage* with an intact and re-inducible anthocyanin toolkit, and show that its regain is heritably maintained?

**Not identifiable from the presently registered sources.** This is a source adequacy result. It is neither biological evidence that white-flower pigment competence was lost, nor a negative estimate of recovery frequency.

| Evidence source | Existing material | Permitted interpretation | Crucial missing evidence |
|---|---|---|---|
| PRJNA1311153 | 25 RNA-Seq runs / 25 BioSamples / 8 titled taxa, all **young leaf** | accession- and specimen-code-backed leaf RNA opportunity | 0 flower-stage RNA-Seq runs in this source; no joint floret pigment/enzyme/fitness assay |
| Ryukyu subset | *C. brevicaule* 3, *C. irumtiense* 3 leaf runs | six original run identities | no ancestral pigmentation polarity, functional reactivation or phenotype proven per individual |
| *C. japonicum* var. *takaoense* | W-labelled leaf libraries 3 and BP-labelled leaf libraries 3 | morph-coded source contrast | geographic confounding, no flower RNA/metabolites or independent population replication |
| Park (2020) *C. japonicum* | 7 NJ/FastTree reference-family comparisons of five distinct transcript records | all 7 between reference-family anchors, uncertain DFR/ANR and ANS/FLS placement | zero verified catalytic assignments, no measured Ryukyu flower alleles |
| Future Phase A | 6 targeted taxa, **298 intended individuals** (75/75 focal Ryukyu) | prospective sampling plan | none newly collected by P04; no heritable regain events |

The union of the eight actual source taxa and six planned Phase-A taxa contains **12 distinct named taxa**. They are NOT 12 independent origins, field experiments or clades. The exact [machine receipt](../data/anto_source_handoff/p04_afterlife_public_source_readiness_v1.json) records the statuses for all 12.

The source inputs are the existing [25-run metadata manifest](../data/anto_source_handoff/public_sra/prjna1311153_run_manifest_20261017.csv), [QC](../data/anto_source_handoff/public_sra/prjna1311153_qc_20261017.json), [takaoense provenance](../data/anto_source_handoff/public_sra/takaoense_aza3_morph_crosswalk_20261017.csv), [NJ/FastTree molecular result](../data/anto_source_handoff/molecular/cirsium_2020_nj_fasttree_method_comparison_v5.json), and [Phase-A panel](../data/planning/aza3_nature_genomic_sampling_panel_v1.csv). All five Git blob SHA values are pinned in the [source-only audit script](../analysis/anto_empirical/p04_afterlife_source_readiness_v1.py) and frozen output.

## Why additional phylogeny algorithms would not yet answer this question

The same genuine Park (2020) C. japonicum transcript candidates were investigated with MAFFT/NJ and WAG FastTree. The 7/7 matching unresolved family placements do **not** imply a true intermediate enzyme, nonfunctional pigment protein, reversible floral gene, or any independently inferred evolutionary transition. The recorded residue/site-bootstrap comparisons similarly do not identify catalytic activity. Repeated analysis of the same source cannot substitute for new discriminating observations.

See the existing [real gene-family audit](CIRSIUM_FLAVONOID_GENE_FAMILY_UNCERTAINTY_V5.md).

## Required biological bridge from ANTO v10/v10b

Both [ANTO v10](https://github.com/zuizui0223/anto/blob/main/docs/EVOLUTIONARY_AFTERLIFE_V10_EXACT_FRONTIER_AND_UNIVERSALITY.md) and [v10b](https://github.com/zuizui0223/anto/blob/main/docs/AFTERLIFE_V10B_INDUCTION_IDENTIFIABILITY_AND_EMPIRICAL_GATES.md) are **conditional model findings**: a visually lost trait may hide a functional toolkit, but even one induction experiment does not identify heritable continuation.

The following **five distinct measurement layers** are required, on properly matched biological units:

1. Flower colour itself, by calibrated floret VIS/UV, organ/stage and anthocyanin chemistry, not RGB-white alone.
2. Same-individual functionally verified modules: structural coding integrity, DFR/ANR and ANS/FLS biochemical/product tests with suitable controls, flower-specific regulatory capacity; not transcript family matches alone.
3. Inducibility: a controlled stage-matched floral perturbation, with chemistry and RNA confirming actual floral pathway activation and environmental plasticity controls.
4. Heritable continuity: segregation or matched descendants, independent lineage/population replication and an explicitly defined multiple-generation persistence outcome; one-off flower pigment is not heritable regain.
5. Historical and ecological origin: rooted ancestral histories, local ancestry, cytotype, immigration/introgression/de novo alternatives, and final viable filled achenes under a registered ecological mechanism.

Do not equate 40 planned C. sieboldii specimens with 40 independent regained trait origins. The larger 298 target and all current sample/legal gates remain unchanged.

## Inference states that must not be silently promoted

- LEAF_RNA_METADATA_ONLY: genuine sequencing opportunity but no floral biochemical capacity.
- PUBLIC_HOMOLOGY_ONLY: protein identity/alignments without functional orthology or reaction.
- NO_FLORAL_ASSAY: stage-specific pigment/transcription not measured.
- HISTORY_NOT_POLARIZED: trait ancestry, origin and gene flow unresolved.
- INDUCTION_ONLY: pigmentation observed after a perturbation but its inheritance unmeasured.
- FUNCTIONALLY_COMPETENT_BUT_NO_REEXPRESSION: a **possible genuine adverse result only if** the functional toolkit and induction were measured.
- LOSS_OF_FUNCTION_VERIFIED: requires direct evidence that the genetic-developmental resource no longer functions.
- EVOLUTIONARY_REGAIN_SUPPORTED: separately established loss ancestry, justified molecular origin and causal inherited floral trait restoration.
- FITNESS_CONFIRMED: independently validated ecological interaction and viable-achene fitness endpoint.

Current admitted public evidence does not reach either functional or restored-trait levels; missing assay data are UNKNOWN rather than an adverse biological response.

## Reproduction and editorial boundary

Python stdlib with exact source identity checks:

    python -m unittest discover -s analysis -p test_p04_afterlife_source_readiness_v1.py -v
    python analysis/anto_empirical/p04_afterlife_source_readiness_v1.py --expected data/anto_source_handoff/p04_afterlife_public_source_readiness_v1.json --output /tmp/p04-readiness.json

The next paper-relevant biological result requires an independently established flower-stage mechanism and origin, rather than another source-provenance rearrangement. No new field sampling, local permissions, or genetic intervention is authorized by this audit.
