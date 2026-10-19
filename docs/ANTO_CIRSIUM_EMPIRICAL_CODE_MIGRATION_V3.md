# Canonical Cirsium empirical execution moved from anto to aza3 (v3)

**2026-10-19.** This follows the merged [empirical design/data transfer](https://github.com/zuizui0223/aza3/pull/45) and [higher-order head LEGO assay plan](https://github.com/zuizui0223/aza3/pull/46).

anto retains **general evolutionary accessibility theory and cross-angiosperm comparisons**, not an active Cirsium gene-analysis pipeline. All new Cirsium genotype/phenotype/ecology results must be added here with source provenance. Historical anto code and snapshots remain archived to preserve the original research history.

## Existing public source data and executable analysis now have one owner

- [Cross-checked 2026 Cirsium accession sources](../data/anto_source_handoff/public_sra/) — 25 registered young-leaf RNA runs including six Ryukyu source samples (3 C. brevicaule, 3 C. irumtiense) and 6 *takaoense* white/blue-purple external morph labels. Neither project measured same-plant floral pathway induction in these 25 libraries.
- [Real 2020 Cirsium japonicum 26 transcript sequences](../data/anto_source_handoff/molecular/cirsium_japonicum_2020_26_source_sequences_v6.fasta) — original publisher BLAST-annotated DFR/ANS-like **candidates**, source-verified against an XLSX annotation table and sequence checksums. These are *not* alleles in the six Ryukyu samples or automatically 26 orthologous functional loci.
- [Molecular reference candidates](../data/anto_source_handoff/molecular/asteraceae_reference_cds_v5.fasta) — *Saussurea medusa* DFR and ANS and *Gerbera hybrida* GMYC1 independent GenBank references (NOT focal Ryukyu genes).
- [Organ-specific 2020 Cirsium gene/metabolite observations](../data/anto_source_handoff/molecular/cirsium_japonicum_organ_evidence_v6.json) — CHS flower/leaf ratios and DFR root/leaf ratios remain distinct, with no invented floral anthocyanin functional claim.

### The active empirical analysis modules

The canonical processing modules are under [analysis/anto_empirical](../analysis/anto_empirical/):

- cirsium_run_intake_v3.py: reject swappable run/BioSample/taxon or false per-plant flower pigment values;
- cirsium_molecular_reference_v3.py: reference CDS reading frame and candidate validation, never call true coding gene loss from a leaf RNA absence;
- cirsium_organ_evidence_v3.py: source qRT-PCR/metabolite contrasts, no pooling of different tissue denominators;
- cirsium_source_homology_v3.py: six-frame protein similarity vs two independent *Saussurea* DFR/ANS references, candidate ranking without orthology/function claim;
- inspect_public_cirsium_unigene_v3.py and extract_2020_26_genes_v3.py: optional bounded public supplementary archive retrieval, with the original 26 exact author annotation/source entries and uncompressed file hashes. No unapproved 80.66 GB Ryukyu read download.

New **single entrypoint:** [analysis/run_anto_cirsium_empirical_sources_v3.py](../analysis/run_anto_cirsium_empirical_sources_v3.py). It executes full public metadata intake and 26-source transcript similarity triage, with explicit role separation and biological claim ceilings. Unit tests and CI reject source mixing or upgrade to an unsupported Ryukyu anthocyanin re-expression conclusion.

## Reproducibility

    python -m pip install "biopython>=1.80,<2"
    python analysis/run_anto_cirsium_empirical_sources_v3.py --out results/aza3_cirsium_sources_v3.json
    python analysis/test_anto_cirsium_empirical_sources_v3.py -v

**Decision boundaries:** Even if many *C. japonicum* source transcripts show strong DFR-like alignment against *Saussurea*, neither functional orthology, the number of independent inherited pigment genes, nor white-floret reactivation is established without locus genealogy, biochemical activity and genotype-by-flower expression on same individuals. Apparent stop codons from RNA assemblies are validation candidates, not a verified causal gene disruption. Genome history and morphological LEGO will need the aza3 Phase-A and Phase-B sampling/experiment programme (which remains prospective and unauthorized until permissions).

### Upstream-to-downstream bridge

- anto's mathematical result: even full first- and second-order module availability data can leave a three-way required-resource probability unidentified. Plant movement and region-dependent loss affect relative composite ability.
- aza3's direct test: identify joint functional resource retention in a single capitulum and its source genome; demonstrate independent or coupled module origin and developmental compatibility; evaluate head-specific ecological/fitness outcomes under controlled replication.
- Neither colour recurrence nor capitulum morphology covariance alone establishes the ancestral-genetic LEGO reuse mechanism. Further research might find **no** functional reusable toolkit — an informative falsification rather than an automatic failure of the theory.

## Literature/source provenance

- [Chang et al. 2026](https://doi.org/10.1186/s12870-026-08097-6), accession PRJNA1311153
- [Park et al. 2020](https://doi.org/10.1021/acsomega.0c04001), public supplementary FASTA and BLAST XLSX via Europe PMC PMC7675961
- [anto v7 general theory](https://github.com/zuizui0223/anto/blob/main/docs/JOINT_MODULE_ACCESSIBILITY_AND_DISPERSAL_V7.md)
- [aza3 functional LEGO assay plan](ANTO_HIGHER_ORDER_LEGO_MEASUREMENT_PLAN_V2.md)
