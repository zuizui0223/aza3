# GloBI Cirsium genus-wide interaction screen: verified subset and audit gate

**Status:** 2026-10-08. Partial source-backed verification; global bidirectional API query submitted to GitHub Actions, but raw all-Cirsium result **not yet confirmed**. Any count below is for one specified archived GloBI source, never the global *Cirsium* interaction network.

## Why this screen comes first

Prior EAzami and aza3 hypotheses ask whether capitulum orientation, phyllary posture, true spines and stickiness filter pollinator, ovipositing seed predator, florivore, parasitoid, predator and aphid/ant pathways. GloBI offers a discovery layer across organisms, but its edges are **claims under a source-specific relation ontology**; they do not automatically demonstrate head-local use, pollen transfer, seed damage or the adaptive effect of a floral trait.

The full-bidirectional extractor and artifacts are specified in `analysis/collect_globi_cirsium_interactions_v1.py` and `.github/workflows/census-cirsium-globi-v1.yml`. It queries `sourceTaxon=Cirsium` and `targetTaxon=Cirsium` with observation provenance, writes full and provenance-deduplicated CSV tables and a page-ledger of URLs and SHA256 digests, and prints `COMPLETE_DISCOVERY_QUERY` only if both directions are exhausted.

## Confirmed GloBI-indexed source-level examples (NOT a global census)

| Dataset/review identity | Cirsium focal | Indexed partner/type | Scope |
| --- | --- | --- | --- |
| USGS Pollinator Library, indexed GloBI review `zenodo.org/records/16416566`, source DOI 10.5066/P9DSS3VL | *C. arvense* | *Apis mellifera* `visitsFlowersOf`, **423 source-indexed claims** of this pair | USGS source's 5,687 indexed visit records; *C. arvense* listed as 454 target records in this source only |
| Same USGS source | *C. vulgare* | *Bombus* sp. `visitsFlowersOf`, **33 source-indexed claims**; genus-level bee identification | *C. vulgare* has 78 target records in that dataset, not a global count |
| UT–UTBFL, GloBI review `zenodo.org/records/20600916` | *C. texanum* | *Melissodes coreopsis* `interactsWith`, **9 source-indexed claims** of this pair | 41 `C. texanum` associated taxon rows of a 566-row mixed data source; "interactsWith" is not formal demonstrated effective pollination |
| Pocock 2021 indexed GloBI review `zenodo.org/records/18429892` | *C. arvense* | `eatenBy` *Harpalus rufipes*; **2 source-indexed claims** of this pair | Plant parts/stage unspecified by this page; cannot classify as head florivory or pre-dispersal seed predation |
| Dorey 2023 indexed GloBI review `zenodo.org/records/16416973` | *Cirsium* sp. | Bee-associated `hasHost`/related archive, *Cirsium* sp. has **470 associate records** | Plant identified only to genus; host relation does NOT itself identify flower/foliar use or effective pollination |

Do not sum dataset totals: sources may overlap, treatment of record count varies, and the GloBI live cross-source query can contain mirrors/adapters.

Additional GloBI source archives verified in the independent cross-check:
- Bees of Ireland (https://zenodo.org/records/20547080), `Cirsium vulgare` has 12 source-indexed `visits` claims, with `Bombus lapidarius visits C. vulgare` appearing **4** times. `visits` is a general interaction claim, not confirmed effective pollination.
- US Mountain West floral visitor programme (https://zenodo.org/records/20548388), `C. arvense` is present with **76** source-level `visitsFlowersOf` records; partner-specific breakdown is not established by the top-20 pairs on the review page.
- `globalbioticinteractions/gbif-us-bees` (https://zenodo.org/records/20547850) has **1,283** `Cirsium sp.` associated source rows, explicitly lacking species resolution. Neither the 1,283 nor the visitor abundance may be treated as *Cirsium* species richness or a population-level interaction rate.

**Crucial provenance firewall:** `globalbioticinteractions/refuted-biotic-interactions-by-eol` (https://zenodo.org/records/16417164) is an archive of *refuted*, not affirmative, claims. It includes hundreds of `Cirsium arvense`-labelled statements, but **must never** be counted as confirmed pollination, herbivory, parasitism or any positive ecological relation. Exclude that source from all affirmative partner summaries unless examining contradictions separately. `globalbioticinteractions/mangal` may mirror independent Mangal records and is not an independent biological source without a matching DOI/study-level deduplication.

**Source pages (review- and source-identifying):**
- https://zenodo.org/records/16416566
- https://zenodo.org/records/20600916
- https://zenodo.org/records/18429892
- https://zenodo.org/records/16416973
- GloBI data & recommendation for stable versioning: https://www.globalbioticinteractions.org/data
- GloBI official field schema: https://api.globalbioticinteractions.org/interactionFields
- GloBI official interaction-type ontology: https://api.globalbioticinteractions.org/interactionTypes
- GloBI concept archive DOI: https://doi.org/10.5281/zenodo.3950589

## Independent ecological verification (NOT claims of occurrence in GloBI)

The following are **publication-confirmed focal interaction names** to cross-check against GloBI; do not assert their inclusion in a given archived GloBI source until a retrieved row verifies it.

- *Cirsium arvense* capitula: seed-feeding tephritids *Xyphosia miliaria*, *Terellia ruficauda*, *Urophora stylata*; parasitoids *Torymus chloromerus*, *Pteromalus elevatus*. Walker, Hartley & Jones 2008, DOI:10.1111/j.1365-2656.2008.01406.x. The parasitoids attack the **tephritids**, not Cirsium directly. A GloBI claim connecting parasitoid to Cirsium alone, if present, would need original-paper adjudication.
- *Cirsium palustre*: tephritid seed-predator *Terellia ruficauda*; parasitoids *Pteromalus elevatus* and *Torymus chloromerus*; experimental root herbivory shifts seed predator and parasitoid pressure through plant-mediated mechanisms. DOI:10.1007/s004420000569.
- *Cirsium pitcheri*: non-native *Larinus planus* oviposition/seed damage; infested heads about 60% fewer mature seeds. Gijsman et al. 2020, article `S2351989419307243`.
- *Cirsium purpuratum*: visits by *Bombus diversus* are well documented; adult visitation alone is not an estimate of pollen transfer or reproductive fitness. Makino et al. 2007, DOI:10.1111/j.1365-2435.2006.01211.x.

## Source-verified head network anchor matrix (current completed product)

The reproducible **bounded source-check** is now committed at `data/evidence/cirsium_capitulum_interaction_source_anchors_v1.csv`, validated separately by `analysis/validate_cirsium_interaction_source_anchors_v1.py`.

This is a hand-adjudicated source subset, **not** a full GloBI census:
- **28** taxon–source-role rows across **7** named *Cirsium* taxa;
- **4** rows from individually identified GloBI source review pages and **24** independently verified primary-literature rows;
- **21** rows explicitly within the capitulum, in **5** *Cirsium* species; counts are bibliographic anchors, **not** network degrees or replicated observations;
- two indispensable organ controls: `C. arvense × Urophora cardui` is a **stem** gall, `C. arvense × Harpalus rufipes` is a **seed-eating carabid** claim, not certified pre-dispersal head feeding.

The strongest within-head parasitoid systems currently supported by *primary literature*, not asserted present in GloBI:
- `C. arvense`: *Xyphosia miliaria*, *Terellia ruficauda* and *Urophora stylata* in heads, with *Torymus chloromerus* and *Pteromalus elevatus* as parasitoids of the tephritid assemblage (Walker et al. 2008, DOI 10.1111/j.1365-2656.2008.01406.x). Parasitoid taxa are **not** assigned to every specific fly species without a source-verified per-host statement.
- `C. palustre`: *T. ruficauda* plus parasitoids *P. elevatus* and *T. chloromerus*; Masters et al. 2001, DOI 10.1007/s004420000569. Root-herbivore treatments altered parasitoid and fly numbers but did **not** demonstrate a significant difference in parasitism percentage.
- `C. eriophorum`: a diverse source-compiled assemblage of head tephritids and weevils; *Pteromalus vibulenus* reported parasitizing *Rhinocyllus conicus* in this thistle, plus *Bombus lapidarius* as a visitor and wasps physically attacking phyllaries (Tofts 1999, DOI 10.1046/j.1365-2745.1999.00369.x).

### Why a second GloBI query is essential

Querying *Cirsium* on either side identifies a **plant–partner** relation, but a record such as `Torymus chloromerus parasitoidOf Terellia ruficauda` contains **no plant** and therefore will not appear in either plant-centric query.

The new `analysis/collect_globi_cirsium_seed_feeder_secondhop_v1.py` uses the anchored **10** binomial seed-feeder names and queries GloBI in both directions for their parasitoid host relations. Crucially:
- a host–parasitoid GloBI record cannot independently certify that the pair co-occurred in *Cirsium*;
- a primary study must verify the same head, location, season and interaction before the three-trophic plant–host–parasitoid chain is promoted;
- absent second-hop records cannot establish absent parasitoids;
- the workflow `.github/workflows/validate-cirsium-globi-bridge-secondhop-v1.yml` contains separate offline source/semantic validation and a bounded live second-hop acquisition with archived outputs.

### Strong positional control: a capitulum is not exchangeable with another head on the same plant

Gijsman, Havens & Vitt (2020), DOI 10.1016/j.gecco.2020.e00945, found that *C. pitcheri* **terminal** capitula were largest and most productive; weevil *Larinus planus* evidence was concentrated disproportionately on **secondary and tertiary** heads, and infested capitula produced **60% fewer mature seeds** than heads without evidence of weevil infestation. These observations do not prove an angle or spine effect.

Therefore an O × S natural-history analysis must also record **head rank on the flowering stem, branching position, emergence order, within-plant capitulum size and flowering stage**, before treating a spatially exposed or drooping head as an adaptive orientation morph. Head position is a biologically meaningful allocation/predation axis, not a disposable nuisance.

### Strong ecological discriminator: stage-specific defence

A head may exclude adult seed-feeder oviposition while simultaneously restricting the later ovipositors of parasitoids of established larvae. The **fitness sign** of spines, phyllary access and orientation then depends on the *relative timing* of (i) seed-feeder egg placement, (ii) successful early parasitoid attack, and (iii) irreversible seed damage. In particular:

`net damaging larvae ≈ established eggs × [1 − parasitoid mortality before damage] × damaging potential`,

only under the explicit approximation that early parasitism prevents future feeding. A barrier that suppresses initial egg laying might nonetheless increase damage if it reduces **early effective parasitism** sufficiently more. This is not demonstrated in *Cirsium* and must be directly falsified.

An **adult flower visitor** and a **larval seed consumer** are distinct life stages, even if they belong to one insect species. The same named partner must not be assigned a static beneficial/harmful sign.

## Analysis tiers after the full GloBI export

**Level 0: evidence existence.** Retain all exact source and target taxon strings, relation, provenance, date, geography, source namespace, body-part and life-stage fields. Distinct partner names are merely *species-rank name candidates*, not validated accepted species richness.

**Level 1: source semantics.** Distinguish flowersVisitedBy/visitsFlowersOf and pollinatedBy/pollinates (claim only); `interactsWith` and `coOccursWith` are nonspecific; `hasHost` does not resolve larval feeding, and `parasitoidOf` must be interpreted as a consumer–host relationship distinct from the plant.

**Level 2: organ localization.** Direct head/bract/floret/receptacle/achene involvement is required for capitulum-level inference. Leaf miners, root feeders and post-dispersal seed eaters remain relevant to whole-plant ecology but cannot fill focal head-guild cells. Many GloBI records lack body-part and life-stage annotation.

**Level 3: three-trophic chains.** Join verified original sources:
`Cirsium capitulum -> insect seed predator -> insect parasitoid/predator`, preserving access timing and evidence that parasitoid suppression can precede seed damage. Never infer a plant–parasitoid edge as evidence of seed rescue.

**Level 4: orientation/armature adaptation.** GloBI cannot establish that head orientation, true spine length/sharpness, phyllary posture or stickiness changed any interaction probability. That requires replicated same-head manipulations or similarly defensible exogenous variation with reproductive output.

### Mechanistic shortlisting: which Cirsium × guild cells matter first

1. Wild *Cirsium* with all three roles documented for the **same head and site**: legitimate visitor, head seed eater, parasitoid.
2. Species with contrasting nodding/erect head stage and accessible/inaccessible involucre, preferably in the **same population**, to separate ontogeny from actual trait variation.
3. Sites with two damaging enemy modes (external adult ovipositor vs pre-existing internal larvae/stem crawler) to test whether an outward spine barrier is selective.
4. Japanese *Cirsium* data are especially sparse in the current public literature collection; do not transfer northern temperate *C. arvense* mechanisms to Ryukyu *C. brevicaule/C. irumtiense* without a geographic source/functional check.

### Current decision

`PARTIAL_GLOBI_SOURCE_AUDIT_PASS__GENUS_WIDE_ENUMERATION_PENDING`

We have source-confirmed GloBI interaction examples and an executable source/provenance-aware genus query, but **no certified whole-genus interacting-species census** yet. A GitHub Actions queued/pending retrieval is not a scientific observation and should never be described as completed results. Do not declare organismal group absences or joint guild prevalence from this evidence layer.
