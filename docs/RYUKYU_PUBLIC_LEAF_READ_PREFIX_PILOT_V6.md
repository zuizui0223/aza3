# aza3: bounded direct Ryukyu public RNA read test (v6; October 19, 2026)

**Owner:** aza3, NOT anto. anto remains the general theory of history-dependent latent evolutionary potential. This is an **actual public-data access protocol**, not a flower-colour, genotype, enzyme-catalysis or evolutionary-reversal result.

## Current physical source and exact sample grain

Chang et al. (2026, DOI 10.1186/s12870-026-08097-6) registered 25 young-leaf paired RNA-seq experiments under [NCBI BioProject PRJNA1311153](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1311153). ENA has **six** source-identified Ryukyu runs: three from the white-flowered taxon *C. brevicaule* and three from coloured-taxonomic *C. irumtiense*. Those published taxon-level flower descriptions must **not** be assigned as six experimentally measured individual corolla pigment outcomes.

| Published species | Sample code | ENA Run | R1 compressed GB (provider) |
|---|---|---|---:|
| *C. brevicaule* | OKN-4295 | SRR35152725 | 1.605 |
| *C. brevicaule* | YRN-4166 | SRR35152729 | 1.413 |
| *C. brevicaule* | AMM-4163 | SRR35152730 | 1.531 |
| *C. irumtiense* | YGN-4296 | SRR35152724 | 1.589 |
| *C. irumtiense* | ISG-4111 | SRR35152731 | 1.640 |
| *C. irumtiense* | MMY-4078 | SRR35152732 | 1.876 |

The six full R1 gzip files total ~9.66 GB; both paired mates together are about double that. The earlier full 25-run ENA pair listing totals about 80.66 GB. Nothing in this pilot obtains those full data.

The official URL, BioSample, specimen code, full file size, provider-declared MD5, read-count metadata, and explicit tissue limits are frozen in [public R1 source manifest](../data/anto_source_handoff/public_sra/ryukyu_six_R1_prefix_source_v6.json). The full gzip MD5 cannot validate an 8 MiB prefix; the code records a SHA256 of the actual downloaded prefix instead.

## Predeclared bounded run

- Exactly **first 8 MiB compressed bytes from R1** per approved SRR: maximum 48 MiB combined, with HTTP Range requested and strict stream caps even if Range is ignored.
- Incrementally decode the truncated gzip stream, discarding any final incomplete FASTQ record; record counts, observed read lengths and base-quality summaries. A partial gzip stream not reaching checksum EOF is **expected**, not a corrupt whole sample.
- Align *only* those complete R1 reads to the **26 source-confirmed 2020 *C. japonicum* transcripts** previously annotated as DFR/ANS-like, **plus 26 nucleotide-composition-matched shuffled sequence decoys**. The decoys assess a narrow artificial sequence-mapping null, not complete false discovery or orthology.
- Map using minimap2 short-read preset (-ax sr -k 15 -w 5), record all primary alignments, MAPQ and edit distances, and store only transcript-level aggregates.
- Define a **stringent sequence similarity candidate** as a primary minimap2 alignment with at least 85 aligned read bases, at least 85% approximate edit-distance identity and MAPQ >=20. Preserve weak/partial alignments separately and do not suppress the decoy results.
- No unplanned download of R2, no complete public FASTQ, no genotype or plant origin inference, no raw reads committed, no private GPS/site data. All temporary reads and SAM files are deleted after analysis. Run-specific failures produce explicit status, never negative gene-calls.

The method uses **first-in-file** read prefixes, *not* true random sampling and *not* coverage-balanced read subsampling. Read order, library treatment, expression bias, mappability to only 26 known transcripts, sequence divergence and paralogs can all change detection. A zero cannot establish structural gene loss or a lack of flower-specific anthocyanin pathway competence. Likewise, observed mapping can be caused by related protein-family genes, contaminants, RNA assemblies or introgressed alleles; it cannot yet demonstrate intact plant enzyme activity.

## Hypotheses and STOP labels

- **Target detectable:** a subset of young-leaf reads from a given specimen aligns to a source-verified 2020 DFR/ANS-like transcript under frozen nucleotide similarity thresholds. Highest claim: LEAF_PREFIX_HOMOLOGY_CANDIDATE.
- **No strict candidate:** there were no high-quality alignments in this first 8 MiB R1 prefix. Highest claim: NOT_DETECTED_IN_BOUNDED_PREFIX, **not** no gene.
- **Decoy mappings:** if shuffles attract strict mappings, the naive prefix reference search is not specific; must revise controls before interpreting any source signal.
- **No data retrieved:** server may reject HTTP Range or connection. Highest claim: DATA_ACCESS_BLOCKED_OR_FAILED, not no transcript.
- **All samples:** no differential gene expression estimates between species may be made from these partial single-mate prefix counts or without a proper full-library model.
- **Functional distinction:** the previous real gene-family comparisons (NJ and FastTree) were unresolved between DFR/ANR and ANS/FLS, with seven AA-identity 95% intervals spanning zero. Even a read homolog of the 2020 *C. japonicum* candidates does **not** change this conclusion.

## If the capped read test yields signal

Advance **only in aza3** to a source-aware, feasible small cohort with proper whole-read mapping on paired reads if resource cost/authorization allows; verify the focal sequence family by independent gene phylogeny against ANR/FLS controls and genomic allele ancestry. Follow with same-individual corolla (flower-stage) and involucre/leaf pigment chemistry, structural enzyme assays and regulatory expression, separate from ancestral colour polarity and gene-flow reconstruction. Validate viable-achene/ecological hypotheses only after head modules become functional/genomic measurements.

This remains one **mechanism candidate test** within aza3's larger head LEGO experiment. anto can eventually use a frozen, uncertainty-bearing lineage-level latent resource estimate, but cannot take a leaf-prefix read match as an indicator that a silent flower can re-evolve a phenotype.

## Reproduce

    # Ubuntu / GitHub Actions
    sudo apt-get install minimap2
    python analysis/test_ryukyu_bounded_read_prefix_v6.py -v
    # Actual public source, capped to 8 MiB per R1 file; six runs max ~48 MiB:
    python analysis/anto_empirical/ryukyu_bounded_read_prefix_v6.py --fetch-public --cap-mib 8 --output results/ryukyu_2026_leaf_prefix_v6.json

The GitHub CI runs the public read retrieval only from this frozen experiment with explicit resource caps and stores only source-compliant numeric evidence as an artifact. Full biological inference remains on hold.
