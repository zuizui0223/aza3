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
