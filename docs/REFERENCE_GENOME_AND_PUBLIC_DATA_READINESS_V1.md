# aza3 Gate R1 — reference genome and paper-linked public-data readiness v1

**Status:** R1A RESOURCE AUDIT COMPLETE — R1B TRANSFERABILITY PILOT PENDING  
**Date:** 2026-10-05

## Why this gate exists

The Nature-scale genomic claim is not a generic population-structure result. It requires defensible inference that different capitulum modules can follow separable local genomic histories in the same populations or individuals.

That inference is vulnerable to cross-species reference bias. A congener reference can be adequate for chromosome coordinates and genome-wide ancestry while still creating false window-level differences in callability, allele balance, local PCA, local genealogy or structural-variant discovery.

Therefore the first genomic gate is not:

> is there a Cirsium reference genome?

It is:

> can the existing public Cirsium resource set support the exact local-history estimands required for the focal Japanese system without reference-induced artefact?

## R1A — public resource audit

R1A was expanded beyond genome databases to include article-linked supplementary and archival data.

Canonical machine-readable inventory:

`data/planning/reference_public_resource_inventory_v1.csv`

### 1. Chromosome-scale congener resources

**Cirsium heterophyllum**
- hap1 chromosome assembly: `GCA_965225835.1`
- hap2 scaffold assembly: `GCA_965225975.1`
- raw genomic/transcriptomic project: `PRJEB84092`
- umbrella project reports about 174 Gb of SRA data
- Ensembl Genebuild is available for the hap1 assembly

**Cirsium dissectum**
- hap1 chromosome assembly: `GCA_965276805.1`
- hap2 scaffold assembly: `GCA_965276745.1`
- raw genomic/transcriptomic project: `PRJEB85033`
- umbrella project reports about 167 Gb of SRA data
- Ensembl Genebuild is available for the hap1 assembly

NCBI Datasets resolves both hap1 assemblies to **17 chromosome-scale scaffolds**, consistent with a haploid complement of 17:
- `C. heterophyllum` hap1: 944,019,485 bp; contig N50 2,674,182 bp; scaffold N50 50,939,392 bp;
- `C. dissectum` hap1: 1,069,227,339 bp; contig N50 3,154,544 bp; scaffold N50 59,465,939 bp.

The earlier 19-chromosome description was incorrect and is superseded by the accession-level NCBI assembly report.

The hap2 assemblies are retained as **same-individual reference controls**. This is particularly important for `C. dissectum`, whose NCBI assembly comments explicitly report haplotypic inversions between hap1 and hap2 on chromosomes 1, 6 and 7. A locus-placement result that changes merely by switching haplotypes of the same reference individual cannot be treated as robust evidence of lineage-specific genomic structure.

These are high-value chromosome-coordinate and synteny backbones. They are not assumed to be unbiased focal haplotype references.

### 2. East-Asian whole-genome resource

The 2024 `C. nipponicum` genome paper provides more than a database assembly record:

- raw Illumina + Oxford Nanopore sequencing deposited under `PRJNA1127082`;
- assembled genome size about 929.4 Mb;
- contig N50 about 0.7 Mb;
- BUSCO complete about 95.1%;
- repeats about 70.94%;
- assembled genome and annotations deposited at Figshare DOI `10.6084/m9.figshare.26927092`.

This is currently the most useful public East-Asian whole-genome sequence anchor even though it is not chromosome-scale and is not the focal Japanese radiation.

### 3. Moreyra et al. 2025 genus-wide phylogenomics

Moreyra et al. (Molecular Phylogenetics and Evolution 204:108285) sampled 266 `Cirsium` accessions representing 248 species and inferred the phylogeny from 350 nuclear target-capture loci.

Public raw reads are deposited under `PRJNA957074`.

aza3 use:
- broad taxon identity and nuclear orthology;
- Japan-radiation phylogenetic scaffold;
- topology sensitivity;
- calibration/recovery support where directly justified by the published analysis.

Claim ceiling:
- not a fine local-ancestry panel;
- not a haplotype-age panel;
- not a structural-variant panel;
- not sufficient on its own for module-specific recombination blocks.

### 4. East-Asian phylotranscriptomic public data

The 2025 Nipponocirsium study deposited 13 young-leaf RNA-seq experiments under `PRJNA1158676` (about 89 Gb reported by BioProject).

The 2026 `C. japonicum` complex / `C. brevicaule` group study deposited 25 new young-leaf RNA-seq experiments under `PRJNA1311153` (about 175 Gb reported by BioProject). The paper combines these with previously generated Nipponocirsium/`C. lineare` transcriptomes to analyze 33 `Cirsium` accessions plus outgroups.

The 2026 focal set includes:
- Japanese `C. japonicum` var. `japonicum`;
- Taiwanese `C. japonicum` complex accessions including white and bluish-purple var. `takaoense`;
- `C. brevicaule`;
- `C. irumtiense`;
- allied East-Asian taxa.

aza3 use:
- coding orthology;
- expressed coding haplotype comparison;
- focal-lineage sequence-divergence audit;
- ancestry context for the existing six morph-labelled `takaoense` public samples;
- annotation support if a focal genome is built.

Claim ceiling:
- young-leaf RNA is not a floral regulatory assay;
- transcriptome data do not provide intergenic regulatory sequence;
- they do not provide genome-wide local ancestry;
- geography-confounded colour samples do not pass AV1.

### 5. Focal C. sieboldii status

No public chromosome-scale or whole-genome `C. sieboldii` reference was identified in the current NCBI / Ensembl / Darwin Tree of Life + literature-linked audit.

This is an audit-bounded statement, not a claim that no private or unindexed sequence resource exists.

A focal assembly is therefore a possible downstream requirement, not an automatic first action.

## What changed after inspecting article-linked public data

The earlier binary framing

`public reference exists / does not exist`

is inadequate.

The resource stack is now:

`chromosome backbones (heterophyllum, dissectum)`
+
`East-Asian genome sequence/annotation (nipponicum)`
+
`genus-wide target capture (Moreyra)`
+
`East-Asian coding sequence diversity (PRJNA1158676, PRJNA1311153)`.

This is enough to run a serious transferability pilot before paying for a focal assembly.

## R1B — C. sieboldii reference-transfer pilot

### Pilot material

Target 6–10 `C. sieboldii` individuals spanning:
- more than one population;
- directly measured capitulum states;
- the strongest available colour/orientation combinations;
- same-individual 2C/cytotype when feasible.

The first mapping-QC pass should remain phenotype-blind where operationally possible.

### Run the same reads against

1. `C. heterophyllum GCA_965225835.1`;
2. `C. dissectum GCA_965276805.1`;
3. the public `C. nipponicum` assembly.

Reference-internal sensitivity controls:
- `C. heterophyllum GCA_965225975.1` (hap2);
- `C. dissectum GCA_965276745.1` (hap2).

The hap2 controls do not replace the three primary references. They quantify how much locus placement can change within one diploid reference individual before any among-species reference effect is interpreted.

### Required outputs

Global:
- unique mapping fraction;
- mapped-pair fraction;
- depth and coverage uniformity;
- callable fraction;
- allele balance;
- genotype-likelihood PCA;
- pairwise kinship/relatedness.

Window level:
- missingness;
- local PCA or local relatedness;
- reference-specific signal gain/loss;
- allele-balance distortion.

After QC is frozen, open phenotype labels and test:
- phenotype-linked missingness;
- phenotype x reference callability interaction;
- population-linked mapping failure;
- whether candidate module-associated windows survive reference substitution.

## Decision classes

### GREEN_EXISTING_REFERENCE_SUFFICIENT

Use existing public resources if genome-wide ancestry **and the local-history signal used for the Nature inference** are materially reference-invariant and there is no phenotype-linked mapping/callability artefact.

This authorizes discovery/triage. It does not automatically authorize:
- absolute haplotype age;
- fine causal variant identity;
- structural-variant reuse.

### AMBER_FOCAL_REFERENCE_REQUIRED

If global ancestry is stable but local genealogy, haplotype blocks or genotype-phenotype localization changes materially with the reference, build a focal `C. sieboldii` reference before Tier-2 Nature inference.

Minimum focal-reference requirements:
- reference individual from the focal biological system;
- direct module phenotypes recorded;
- flow-cytometry checked;
- PacBio HiFi or equivalent long-read sequence;
- Hi-C or equivalent chromosome scaffolding;
- repeat annotation;
- gene annotation supported by the public Cirsium transcriptomes.

### RED_PANGENOME_OR_GRAPH_REQUIRED

If reference bias tracks phenotype/ancestry classes, or a focal linear reference still changes the local-history conclusion, use at least two divergent focal assemblies and a graph/mini-pangenome representation before claiming combinatorial genomic reuse.

### NOT_IDENTIFIABLE

If pilot coverage or biological sampling cannot distinguish mapping artefact from local-history signal, do not open module-specific genomic-history claims.

## Nature claim boundary

The central Figure-3 claim is only open when:

`same-individual phenotype`
→ `reference-stable genomic signal`
→ `module-specific local history`
→ `novel module combination`

is supported.

A high overall mapping rate is not enough.

## Immediate execution order

1. Freeze this public-data inventory.
2. Recover exact public sequence assets required for the three-reference pilot.
3. Design 6–10-individual `C. sieboldii` low-cost pilot.
4. Freeze transferability metrics before phenotype-linked genomic results are opened.
5. Classify GREEN / AMBER / RED / NOT_IDENTIFIABLE.
6. Only then authorize either existing-reference Tier 2 or a new focal reference build.

## Primary public sources audited

- Moreyra LD et al. 2025. *A thorny tale: The origin and diversification of Cirsium (Compositae).* Molecular Phylogenetics and Evolution 204:108285. DOI: 10.1016/j.ympev.2025.108285. Raw data: PRJNA957074.
- Cheon et al. 2024. *De Novo Genome Assembly and Phylogenetic Analysis of Cirsium nipponicum.* Genes 15:1269. Raw data: PRJNA1127082. Genome/annotation: DOI 10.6084/m9.figshare.26927092.
- East-Asian Nipponocirsium phylotranscriptomic study, 2025. Raw data: PRJNA1158676.
- Chang et al. 2026. *Phylotranscriptomics and genome-size evidence clarify the Taiwanese Cirsium japonicum complex and delimit C. brevicaule and allied East Asian thistles.* BMC Plant Biology. Raw data: PRJNA1311153.
- Darwin Tree of Life / Wellcome Sanger resources for `C. heterophyllum`: PRJEB84092, GCA_965225835.1, GCA_965225975.1.
- Darwin Tree of Life / Wellcome Sanger resources for `C. dissectum`: PRJEB85033, GCA_965276805.1, GCA_965276745.1.
