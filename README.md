# scdown

[![pkgdown](https://img.shields.io/badge/docs-pkgdown-315c86)](https://dai540.github.io/scdown/)
[![R-CMD-check](https://github.com/dai540/scdown/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dai540/scdown/actions/workflows/R-CMD-check.yaml)
[![GitHub release](https://img.shields.io/github/v/release/dai540/scdown)](https://github.com/dai540/scdown/releases)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/dai540/scdown/blob/main/LICENSE.txt)

`scdown` is a minimal downstream package for annotated single-cell RNA-seq.

It assumes that QC, normalization, clustering, and annotation are already done.
The package is intentionally narrow and focuses on five practical jobs inside
one dataset:

- cluster and annotation maps
- marker discovery
- annotation checks with known marker panels
- cluster-level expression and signature summaries
- exploratory ligand-receptor interaction ranking

Website: <https://dai540.github.io/scdown/>

## Installation

Install from GitHub with `pak`:

```r
install.packages("pak")
pak::pak("dai540/scdown")
```

or `remotes`:

```r
install.packages("remotes")
remotes::install_github("dai540/scdown")
```

Or install from a source tarball:

```r
install.packages("path/to/scdown_<version>.tar.gz", repos = NULL, type = "source")
```

Then load the package:

```r
library(scdown)
```

## Input contract

`scdown()` accepts:

- a long-format table with `cell`, `gene`, `value`, and annotation columns
- a list with `expr`, `cells`, and optional `embedding`
- a `SingleCellExperiment`
- a `Seurat` object

The minimal metadata columns are:

- `cell`
- `annotation`
- `sample` if available

## Minimal example

```r
library(scdown)

obj <- scdown(
  scdown_example(),
  annotation_col = "annotation",
  sample_col = "sample"
)

res <- run_core(
  obj,
  signatures = c("cytotoxic", "antigen_presentation")
)

build_scdown_report(obj)
```

## Core functions

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

## Main outputs

Typical outputs include:

- cluster map PNGs
- marker tables
- annotation check tables and heatmaps
- average-expression summaries
- cluster-similarity heatmaps
- signature score tables
- communication ranking tables
- one compact HTML report

## What scdown does not do

`scdown` does not try to replace upstream or specialist tools.

- It does not do QC, ambient RNA correction, doublet detection, normalization, integration, or clustering.
- It does not do group-aware differential abundance or differential expression.
- It does not claim publication-grade communication inference.

Use Seurat, `scater`, `scuttle`, `scran`, `SingleR`, `muscat`, `dreamlet`,
`miloR`, or CellChat for those stages.

## Documentation

- Home: <https://dai540.github.io/scdown/>
- Getting Started: <https://dai540.github.io/scdown/articles/getting-started.html>
- Reference: <https://dai540.github.io/scdown/reference/index.html>
- News: <https://dai540.github.io/scdown/news/index.html>

## Citation

```r
citation("scdown")
```

## Package layout

- `R/object.R`: input standardization
- `R/downstream.R`: minimal downstream workflow
- `R/resources.R`: built-in signatures, marker panels, and ligand-receptor pairs
- `R/report.R`: compact HTML reporting
- `vignettes/`: getting-started tutorial
