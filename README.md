# scdown

[![R-CMD-check](https://github.com/dai540/scdown/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dai540/scdown/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/dai540/scdown/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/dai540/scdown/actions/workflows/pkgdown.yaml)
[![Latest release](https://img.shields.io/github/v/release/dai540/scdown?display_name=tag)](https://github.com/dai540/scdown/releases)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/dai540/scdown/blob/main/LICENSE.txt)

`scdown` is a minimal downstream toolkit for annotated single-cell RNA-seq.

It is built for the common case where you already have clusters or cell-type
labels and want to understand structure within one dataset without introducing
group- or replicate-aware inference.

## What `scdown` does

- show cluster or annotation maps
- extract marker genes
- check annotations with known marker panels
- summarize average expression by cluster
- inspect cluster similarity
- score signatures or pathways
- rank exploratory ligand-receptor interactions
- build a compact HTML report

## What `scdown` does not do

- QC, ambient RNA correction, doublet detection, normalization, clustering, or integration
- group-aware differential abundance or differential expression
- formal communication modeling for publication-grade claims

Use Seurat, `scater`, `scuttle`, `scran`, `SingleR`, `muscat`, `dreamlet`,
`miloR`, or CellChat for those stages.

## Website

- Home: [https://dai540.github.io/scdown/](https://dai540.github.io/scdown/)
- Getting started: [https://dai540.github.io/scdown/articles/getting-started.html](https://dai540.github.io/scdown/articles/getting-started.html)
- Reference: [https://dai540.github.io/scdown/reference/index.html](https://dai540.github.io/scdown/reference/index.html)
- News: [https://dai540.github.io/scdown/news/index.html](https://dai540.github.io/scdown/news/index.html)

## Installation

```r
install.packages("pak")
pak::pak("dai540/scdown")
```

## Inputs

`scdown()` accepts:

- a long-format table with `cell`, `gene`, `value`, and annotation columns
- a list with `expr`, `cells`, and optional `embedding`
- a `SingleCellExperiment`
- a `Seurat` object

The minimal metadata columns are:

- `cell`
- `annotation`
- `sample` if available

## Minimal workflow

```r
library(scdown)

obj <- scdown(
  scdown_example(),
  annotation_col = "annotation",
  sample_col = "sample"
)

plot_map(obj)
find_markers(obj)
check_annotation(obj)
plot_average_expression(obj)
plot_cluster_similarity(obj)
plot_signature(obj, signatures = c("cytotoxic", "antigen_presentation"))
infer_communication(obj)

build_scdown_report(obj)
```

## Main functions

- `scdown()`
- `summarize_dataset()`
- `plot_map()`
- `find_markers()`
- `check_annotation()`
- `plot_average_expression()`
- `plot_cluster_similarity()`
- `plot_signature()`
- `infer_communication()`
- `run_core()`
- `build_scdown_report()`

## Built-in resources

- `immune_signatures()`
- `immune_marker_panels()`
- `list_signatures()`
- `get_signature()`
