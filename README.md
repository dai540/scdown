# scdown

`scdown` is a simple R package for post-annotation downstream analysis of single-cell RNA-seq data.

The package is intentionally short and practical:

- build one analysis object
- inspect cell maps
- find positive markers
- visualize user-selected genes
- summarize composition
- check annotation with marker panels
- score signatures
- inspect exploratory communication
- export a compact HTML report

It accepts:

- a long-format expression table
- a list with `expr`, `cells`, and optional `embedding`
- `SingleCellExperiment`
- `Seurat`

## Installation

```r
install.packages("pak")
pak::pak("dai540/scdown")
```

## Example

```r
library(scdown)

obj <- scdown(scdown_example())

plot_map(obj)
find_markers(obj)
plot_gene(obj, genes = c("CD3D", "NKG7", "LYZ"))
plot_composition(obj)
check_annotation(obj)
plot_signature(obj, signatures = c("cytotoxic", "exhaustion"))
infer_communication(obj)
build_scdown_report(obj)
```

## Direct input

```r
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

- `scdown()`
- `summarize_dataset()`
- `collapse_lineage()`
- `plot_map()`
- `find_markers()`
- `plot_gene()`
- `plot_composition()`
- `check_annotation()`
- `plot_signature()`
- `plot_heatmap()`
- `infer_communication()`
- `run_core()`
- `build_scdown_report()`
