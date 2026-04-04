# Recommended Pipeline

`scdown` is designed for **annotation-after downstream work**.

That means:

- do QC, normalization, dimensionality reduction, integration, and
  clustering upstream
- do label assignment upstream
- use `scdown` to review the annotated object quickly
- hand off the final formal analyses to specialist tools when needed

## Generic positioning

``` r
kable(scope_table("generic"))
```

| stage | use_scdown | scdown_role | specialist_tools | note | scenario_note |
|:---|:---|:---|:---|:---|:---|
| qc_and_cleanup | no | not designed for this stage | SoupX, scDblFinder, scater, scuttle | Handle droplets, ambient RNA, doublets, and low-quality cells upstream. | Generic single-cell workflow boundary. |
| normalize_cluster_integrate | no | not designed for this stage | Seurat, scran, batchelor, scater | Settle normalization, embedding, clustering, and integration before `scdown`. | Generic single-cell workflow boundary. |
| label_assignment | no | annotation sanity check only | SingleR, Azimuth, manual curation | Use `scdown` to review labels, not to create them from scratch. | Generic single-cell workflow boundary. |
| exploratory_review | yes | core use case | scdown | Maps, markers, genes, composition, signatures, communication, and reports are the main job of `scdown`. | Review any annotated object quickly. |
| lightweight_group_tests | yes | good for quick sample-aware screening | scdown | Useful for prioritizing endpoints, but not a substitute for formal specialist models. | Good for initial prioritization across groups. |
| publication_grade_state_models | handoff | handoff point | muscat, dreamlet | Use specialist models when differential state becomes a primary claim. | Promote key marker or state claims into specialist DE/DS models. |
| publication_grade_abundance | handoff | handoff point | miloR | Use neighborhood- or model-based abundance tools when abundance is central. | Promote abundance claims into specialist DA methods. |
| publication_grade_communication | handoff | handoff point | CellChat | Use a dedicated communication framework when communication is a paper-level result. | Promote communication claims into a specialist communication framework. |

``` r
kable(recommended_pipeline("generic"))
```

| stage | goal | recommended_tools | scdown_role |
|:---|:---|:---|:---|
| raw_qc | remove technical artifacts and low-quality droplets | SoupX, scDblFinder, scater, scuttle | no |
| normalize_cluster | normalize, reduce dimensions, cluster, and integrate batches | Seurat, scater, scran, batchelor | no |
| annotate | assign biologically meaningful cell labels | SingleR, Azimuth, manual curation | no |
| review_with_scdown | inspect maps, markers, genes, composition, signatures, and communication | scdown | core |
| rigorous_testing | run publication-grade differential state, abundance, or communication analyses | muscat, dreamlet, miloR, CellChat | handoff |

The shortest rule is:

- use `scdown` for review, prioritization, and compact reporting
- hand off to specialist tools once an endpoint becomes a primary claim

## What stays inside `scdown`

`scdown` is a good home for:

- maps and annotation review
- quick marker review
- user-selected genes
- sample-level composition summaries
- signature summaries
- exploratory communication
- compact HTML reports

## What should hand off

``` r
kable(scope_table("patient_comparison"))
```

| stage | use_scdown | scdown_role | specialist_tools | note | scenario_note |
|:---|:---|:---|:---|:---|:---|
| qc_and_cleanup | no | not designed for this stage | SoupX, scDblFinder, scater, scuttle | Handle droplets, ambient RNA, doublets, and low-quality cells upstream. | Patient-aware QC must be solved upstream. |
| normalize_cluster_integrate | no | not designed for this stage | Seurat, scran, batchelor, scater | Settle normalization, embedding, clustering, and integration before `scdown`. | Patient-aware normalization, integration, and embedding stay upstream. |
| label_assignment | no | annotation sanity check only | SingleR, Azimuth, manual curation | Use `scdown` to review labels, not to create them from scratch. | Patient labels should be stabilized before comparison. |
| exploratory_review | yes | core use case | scdown | Maps, markers, genes, composition, signatures, communication, and reports are the main job of `scdown`. | Use `scdown` to review and prioritize patient comparison endpoints. |
| lightweight_group_tests | yes | good for quick sample-aware screening | scdown | Useful for prioritizing endpoints, but not a substitute for formal specialist models. | Good for a first pass, not for the final claim set. |
| publication_grade_state_models | handoff | handoff point | muscat, dreamlet | Use specialist models when differential state becomes a primary claim. | Use specialist mixed or pseudobulk models for patient comparisons. |
| publication_grade_abundance | handoff | handoff point | miloR | Use neighborhood- or model-based abundance tools when abundance is central. | Use specialist DA models for patient abundance claims. |
| publication_grade_communication | handoff | handoff point | CellChat | Use a dedicated communication framework when communication is a paper-level result. | Use specialist communication models when patient-group signaling is central. |

When the question becomes formal, export a handoff object:

``` r
obj <- scdown(scdown_example())
export_for_handoff(obj, target = "muscat")
export_for_handoff(obj, target = "dreamlet")
export_for_handoff(obj, target = "miloR")
export_for_handoff(obj, target = "CellChat")
```

These exports do not replace the specialist package. They prepare the
object you would pass into the specialist workflow.

## PBMC

Use `scdown` after standard PBMC preprocessing and annotation.

``` r
kable(recommended_pipeline("pbmc"))
```

| stage | goal | recommended_tools | scdown_role |
|:---|:---|:---|:---|
| preprocess | standard PBMC QC, normalization, and clustering | Seurat or scater/scran | no |
| annotate | assign canonical immune populations | SingleR, Azimuth, manual marker review | no |
| review_with_scdown | rapidly review markers, cytotoxic programs, composition, and report | scdown | core |
| targeted_testing | test a few predefined signatures or compositions across groups | scdown, then muscat or dreamlet if the question becomes formal | optional_test |

This is usually enough when the goal is:

- checking whether annotation looks sensible
- comparing known immune signatures
- sharing a compact report with collaborators

## Tumor immune

Tumor-infiltrating immune datasets are a strong fit for `scdown` as a
review layer, especially for:

- markers
- user-selected genes
- composition
- immune signatures
- exploratory communication

``` r
kable(recommended_pipeline("tumor_immune"))
```

| stage | goal | recommended_tools | scdown_role |
|:---|:---|:---|:---|
| preprocess | clean tumor-infiltrating immune profiles and integrate samples | Seurat, scater/scran, batchelor | no |
| annotate | label T/NK, B/plasma, myeloid, and dendritic compartments | SingleR plus manual tumor-immune curation | no |
| review_with_scdown | review immune state markers, composition, signatures, and exploratory communication | scdown | core |
| rigorous_followup | promote key abundance, state, or communication claims into specialist methods | muscat, dreamlet, miloR, CellChat | handoff |

## Patient comparison

For donor-aware or patient-vs-patient comparisons, `scdown` is still
useful, but mainly as the layer that decides what should be promoted
into formal models.

``` r
kable(recommended_pipeline("patient_comparison"))
```

| stage | goal | recommended_tools | scdown_role |
|:---|:---|:---|:---|
| preprocess | prepare a donor-aware object with clean sample metadata | Seurat or Bioconductor upstream stack | no |
| annotate | stabilize labels before comparison | SingleR, Azimuth, manual review | no |
| review_with_scdown | use `scdown` to decide which endpoints are worth formal modeling | scdown | core |
| formal_models | fit sample-aware models for state, abundance, and pathway differences | dreamlet, muscat, miloR | handoff |

In this scenario the usual handoff is:

- `scdown` for review, prioritization, and compact reports
- `dreamlet` or `muscat` for differential state
- `miloR` for abundance-style shifts
- CellChat or another specialist framework when communication is a
  primary claim

## Minimal practical rule

If the question is:

- “What does this annotated dataset look like?”
- “Which cell types and signatures stand out?”
- “Can I quickly share a report?”

then `scdown` is a good home.

If the question is:

- “Is this the final model I will publish?”
- “Do I need covariates, random effects, or formal abundance modeling?”

then `scdown` should be the handoff point, not the final engine.
