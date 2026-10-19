# Cirsium DFR/ANR and ANS/FLS gene-family ambiguity — real-data audit v5

**Study owner: aza3. Published sequences: Park et al. (2020), *Cirsium japonicum*; **not** the 2026 Ryukyu *C. brevicaule* or *C. irumtiense* plants.**
**Status: actual protein alignments analysed; no enzyme function, floral pigment re-expression, orthology or evolutionary trait reversals established.**

## Scope and successful reproducible analysis

The prior aza3 public-source screens found 26 author-annotated DFR/ANS-like transcript records in 23 assembly groups, reduced by similarity/coverage to five transcript candidates in four assembly groups. Two DFR-like candidates retained their coarse reference-family nearest-match class after reference-taxon omission, but the functional family remained unverified.

v5 adds nine GenBank reference proteins: three DFR, two ANR, two ANS and two FLS, split into independent **DFR/ANR** (short-chain reductase-related) and **ANS/FLS** (2-oxoglutarate/Fe(II) oxygenase-related) protein alignments, never one cross-superfamily tree. The original 2020 transcripts and annotation-sheet provenance remain linked by accession and checksum.

The [actual seven source evaluation records](../data/anto_source_handoff/molecular/cirsium_japonicum_2020_gene_family_topology_evidence_v5.json) were extracted from an independently successful GitHub Actions run after correcting an invalid bootstrap minimum threshold. Alignments were built with MAFFT. Unrooted neighbor joining used pairwise p-distances from nongap shared amino acids, with **200 accepted site-bootstrap samples** per tree. Additional query/reference family identity contrasts were resampled across alignment sites **1,000 times**.

## What the real data support and do not support

| Tree comparison | Cirsium 2020 source transcripts | NJ family-side classification | Relative AA identity, 95% site-bootstrap interval |
|---|---|---|---|
| DFR vs ANR, primary | TR35526|c0_g1_i1 | Between reference family anchors | −0.0238 to +0.0488 |
| DFR vs ANR, primary | TR38875|c0_g4_i5 | Between anchors | −0.0455 to +0.0726 |
| DFR vs ANR, primary | TR39766|c0_g1_i2 | Between anchors | −0.0174 to +0.0770 |
| DFR vs ANR, alternate isoform | TR35526|c0_g1_i1 | Between anchors | −0.0199 to +0.0617 |
| DFR vs ANR, alternate isoform | TR38875|c0_g4_i7 | Between anchors | −0.0538 to +0.0716 |
| DFR vs ANR, alternate isoform | TR39766|c0_g1_i2 | Between anchors | −0.0098 to +0.0796 |
| ANS vs FLS | TR58302|c0_g1_i1 | Between anchors | −0.0490 to +0.0398 |

**Seven out of seven** contrast intervals include zero. The NJ topology classifier returned the same BETWEEN_FAMILY_ANCHORS category in all **seven** query/configuration tests and throughout each associated 200-resample trial.

**Crucial interpretation:** This DOES NOT constitute 100% evidence that a genuinely intermediate enzyme exists! It means that under the sparse set of two biochemical anchor families and this particular unrooted tree classification, the tested query did not nest unambiguously on only one of the two sides. Branch support here conditions on the method and selected source proteins; it cannot assign enzyme function or orthology. Because the identity-difference intervals overlap zero, the protein-family evidence does not reliably distinguish DFR from ANR or ANS from FLS on its own.

The alternative TR38875 source transcripts (i5/i7) belong to one assembled transcript group, so these are **not** two independent biological trait origins or independent functional loci.

## Prespecified reference residue audit — direct evidence, not knockout claims

The source's positions were aligned to the annotated **Saussurea** reference proteins and inspected as individual residues, *not* read as causal variants in the focal 2026 Ryukyu specimens.

- TR39766|c0_g1_i2 (DFR/ANR): five of five inspected reference positions are covered, with **4/5 matching** (G13, G16, Y164, K168) and a G19→A alignment substitution. This is a position-level alignment observation, NOT evidence that an ancestral DFR was knocked out, revived, or catalytically validated.
- TR35526|c0_g1_i1: only the Y164 and K168 positions are aligned, with H and F, respectively. Missing coverage of the glycine-containing N-terminal positions makes a full DFR-motif claim unjustified.
- TR38875 alternative transcripts: Y164 and K168 align and match, but the three N-terminal reference residues lack coverage; enzyme-family identity remains uncertain.
- TR58302|c0_g1_i1 (ANS/FLS): H235 and D237 match the chosen Saussurea ANS reference. These metal-binding-related positions may be broadly conserved across related oxygenases; matching them is not an exclusive ANS functional test.

The protein alignments cannot establish gene origin dates, causal promoter/regulatory changes, cofactor flux, anthocyanin accumulation, ploidy, recombination or effects on seed fitness.

## Second inference method completed: NJ and FastTree agree on unresolved placement

A source-backed second comparison has **now completed successfully** in [GitHub Actions run 39214299802](https://github.com/zuizui0223/aza3/actions/runs/39214299802). Using the SAME published 2020 *C. japonicum* amino-acid alignments, an independent WAG FastTree approximate-maximum-likelihood analysis gives **7/7 identical unrooted reference-family placement labels** to the original NJ: BETWEEN_FAMILY_ANCHORS. The seven evaluations cover **five distinct published transcripts**, not seven genetic loci or trait transitions.

This is an important *negative/identifiability* result: the result is **method-independent ambiguity**, NOT high bootstrap evidence for a distinct enzyme function or proof that any source transcript catalyzes an intermediate reaction. A more complex optimizer cannot rescue missing phylogenetic/functional information by itself.

The [complete, independently reported source comparisons](../data/anto_source_handoff/molecular/cirsium_2020_nj_fasttree_method_comparison_v5.json) are frozen and verified with [strict source and inference gates](../analysis/validate_cirsium_2020_method_concordance_v5.py). The validator also rechecks new FastTree outputs against the seven observed source identities, refusing an upgraded orthology/function statement.

**Decision after two inference methods:** do not yet commit another tree re-estimation without independent additional reference taxa and/or a calibrated catalytic assay. The next discriminating evidence is direct substrate/product enzymology and same-individual floral versus phyllary pathway flux, not increased algorithmic complexity.

## What is scientifically worth doing next (aza3 only)

1. Independently test the **tree inference method** (approximately-ML WAG FastTree vs NJ on *identical MAFFT alignments*). Concordance of ambiguous family sides does not make those sides diagnostic; discordance reveals further source-method sensitivity. Store full source alignment and Newick version/uncertainty metadata.
2. For the relatively complete TR39766 source homolog, expand functionally characterized DFR/ANR references, include site-wise protein domain and structure diagnostics, and avoid function inference from one position.
3. For TR58302, add experimentally validated FLS vs ANS enzyme paralogs and relevant catalyst/substrate discriminators; 2 shared oxygenase sites are not sufficient.
4. **Then** evaluate whether same-population floral pigmentation observations and genomics can link any locus to an actual trait state in a Cirsium individual. Only independently matched cytotype, genetic ancestry, pigmentation chemistry and floral-stage expression/functional tests can address heritable pigment restoration and molecular LEGO reuse.
5. *anto* stays the general history-dependent resource-accessibility programme. Its general mathematical hypotheses are not confirmed or rejected by the inability to assign a 2020 thistle transcript to a specific enzyme reaction.

## Reproducibility

Code: [source-anchored protein trees](../analysis/anto_empirical/cirsium_flavonoid_gene_tree_v5.py), [frozen seven-source record evidence validation](../analysis/validate_cirsium_2020_gene_family_evidence_v5.py), [independent FastTree comparison](../analysis/compare_cirsium_fasttree_v5.py).

Run from aza3 after installing Biopython, MAFFT and FastTree:

    python analysis/test_anto_cirsium_gene_tree_v5.py -v
    python analysis/anto_empirical/cirsium_flavonoid_gene_tree_v5.py --bootstraps 200 --output /tmp/cirsium_gene_families_v5.json --alignment-directory /tmp/cirsium_gene_family_alignments_v5
    python analysis/validate_cirsium_2020_gene_family_evidence_v5.py --analysis /tmp/cirsium_gene_families_v5.json
    python analysis/compare_cirsium_fasttree_v5.py --tree-json /tmp/cirsium_gene_families_v5.json --alignment-directory /tmp/cirsium_gene_family_alignments_v5

The source alignment tests and full results are evidence of **limited functional identifiability under the present public reference panel**, not an experimentally observed limit on trait re-evolution. No new plant field sampling or molecular intervention has been authorized or performed.
