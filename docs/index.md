# scdown

[![R-CMD-check](https://github.com/dai540/scdown/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dai540/scdown/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/dai540/scdown/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/dai540/scdown/actions/workflows/pkgdown.yaml)
[![Latest
release](https://img.shields.io/github/v/release/dai540/scdown?display_name=tag)](https://github.com/dai540/scdown/releases)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/dai540/scdown/blob/main/LICENSE.txt)

`scdown` is a simple downstream toolkit for annotated single-cell
RNA-seq data.

It is best used **after** QC, normalization, embedding, clustering, and
cell annotation. In other words, `scdown` is a downstream hub, not a
full upstream single-cell framework.

## What `scdown` is for

- make one downstream analysis object
- run fast exploratory summaries
- run sample-aware tests only when you need them
- export a compact HTML report for review and collaboration

## What `scdown` is not for

- droplet QC or ambient RNA correction
- normalization, HVG selection, integration, clustering, or trajectory
  fitting
- publication-grade differential state as a primary specialist method

Those parts are better handled by upstream or specialist tools such as
Seurat, `scater`, `scuttle`, `scran`, `SingleR`, `muscat`, `dreamlet`,
`miloR`, and CellChat.

Internally, `scdown` keeps data in a sparse-first matrix form, so it no
longer needs to expand everything into a huge long table before
analysis.

## Website

- Home: <https://dai540.github.io/scdown/>
- Getting started:
  <https://dai540.github.io/scdown/articles/getting-started.html>
- Recommended pipeline:
  <https://dai540.github.io/scdown/articles/recommended-pipeline.html>
- Reference: <https://dai540.github.io/scdown/reference/index.html>
- News: <https://dai540.github.io/scdown/news/index.html>

## Installation

``` r
install.packages("pak")
pak::pak("dai540/scdown")
```

## Citation

GitHub can read
[CITATION.cff](https://github.com/dai540/scdown/blob/main/CITATION.cff),
and R can read `citation("scdown")` from the package.

## Example

``` r
library(scdown)

obj <- scdown(scdown_example())

run_core(obj)
build_scdown_report(obj)
```

## Explore

``` r
plot_map(obj)
explore_markers(obj)
plot_gene(obj, genes = c("CD3D", "NKG7", "LYZ"))
plot_composition(obj)
plot_signature(obj, signatures = c("cytotoxic", "exhaustion"))
infer_communication(obj)
```

## Test

``` r
test_markers(obj)
test_composition(obj)
test_signature(obj, signatures = c("cytotoxic", "b_cell"))
test_communication(obj, n_perm = 25)
```

## Recommended positioning

``` r
scope_table("generic")
recommended_pipeline("generic")
recommended_pipeline("pbmc")
recommended_pipeline("tumor_immune")
recommended_pipeline("patient_comparison")
```

## Specialist handoff

If an endpoint becomes central enough for publication-grade modeling,
prepare a handoff object instead of stretching `scdown` beyond its
intended scope.

``` r
muscat_input <- export_for_handoff(obj, target = "muscat")
dreamlet_input <- export_for_handoff(obj, target = "dreamlet")
milor_input <- export_for_handoff(obj, target = "miloR")
cellchat_input <- export_for_handoff(obj, target = "CellChat")
```

[`build_scdown_report()`](https://dai540.github.io/scdown/reference/build_scdown_report.md)
also includes a scope and handoff section, plus marker and communication
tests when `run_tests = TRUE`. These built-in tests are for screening
and prioritization, not a replacement for specialist models.

## Custom resources

``` r
obj_custom <- scdown(
  scdown_example(),
  signatures = list(my_state = c("CD3D", "IL7R")),
  marker_panels = list(my_t_panel = c("CD3D", "TRAC")),
  lineage_mapping = data.frame(
    pattern = c("cd4", "cd8", "nk", "b", "mono", "macro"),
    lineage = c("T/NK", "T/NK", "T/NK", "B/Plasma", "Myeloid", "Myeloid")
  ),
  lr_pairs = data.frame(
    ligand = "CCL5",
    receptor = "CCR5",
    pathway = "chemokine"
  )
)

plot_signature(obj_custom, signatures = "my_state")
check_annotation(obj_custom, panels = "my_t_panel")
```

## Direct input

``` r
obj_sce <- scdown(
  sce,
  annotation_col = "cell_type",
  sample_col = "sample",
  reduced_dim = "UMAP",
  assay_name = "counts"
)

obj_seurat <- scdown(
  seu,
  annotation_col = "cell_type",
  sample_col = "sample",
  reduced_dim = "umap",
  assay_name = "RNA",
  expression_layer = "counts"
)
```

## Main functions

- [`scdown()`](https://dai540.github.io/scdown/reference/scdown.md)
- [`recommended_pipeline()`](https://dai540.github.io/scdown/reference/recommended_pipeline.md)
- [`scope_table()`](https://dai540.github.io/scdown/reference/scope_table.md)
- [`export_for_handoff()`](https://dai540.github.io/scdown/reference/export_for_handoff.md)
- [`summarize_dataset()`](https://dai540.github.io/scdown/reference/summarize_dataset.md)
- [`collapse_lineage()`](https://dai540.github.io/scdown/reference/collapse_lineage.md)
- [`plot_map()`](https://dai540.github.io/scdown/reference/plot_map.md)
- [`explore_markers()`](https://dai540.github.io/scdown/reference/explore_markers.md)
- [`find_markers()`](https://dai540.github.io/scdown/reference/find_markers.md)
- [`plot_gene()`](https://dai540.github.io/scdown/reference/plot_gene.md)
- [`plot_composition()`](https://dai540.github.io/scdown/reference/plot_composition.md)
- [`explore_composition()`](https://dai540.github.io/scdown/reference/explore_composition.md)
- [`check_annotation()`](https://dai540.github.io/scdown/reference/check_annotation.md)
- [`plot_signature()`](https://dai540.github.io/scdown/reference/plot_signature.md)
- [`explore_signature()`](https://dai540.github.io/scdown/reference/explore_signature.md)
- [`plot_heatmap()`](https://dai540.github.io/scdown/reference/plot_heatmap.md)
- [`infer_communication()`](https://dai540.github.io/scdown/reference/infer_communication.md)
- [`explore_communication()`](https://dai540.github.io/scdown/reference/explore_communication.md)
- [`test_markers()`](https://dai540.github.io/scdown/reference/test_markers.md)
- [`test_composition()`](https://dai540.github.io/scdown/reference/test_composition.md)
- [`test_signature()`](https://dai540.github.io/scdown/reference/test_signature.md)
- [`test_communication()`](https://dai540.github.io/scdown/reference/test_communication.md)
- [`run_core()`](https://dai540.github.io/scdown/reference/run_core.md)
- [`build_scdown_report()`](https://dai540.github.io/scdown/reference/build_scdown_report.md)
