# Create a `scdown` object

`scdown()` standardizes annotation-after single-cell inputs into a
compact object for cluster interpretation within one dataset. Supported
inputs are:

## Usage

``` r
scdown(
  x,
  annotation_col = "annotation",
  sample_col = "sample",
  reduced_dim = "UMAP",
  embedding = NULL,
  assay_name = NULL,
  expression_layer = NULL,
  signatures = NULL,
  marker_panels = NULL,
  lr_pairs = NULL
)
```

## Arguments

- x:

  Input data.

- annotation_col:

  Column containing cluster or annotation labels.

- sample_col:

  Optional column containing sample identifiers.

- reduced_dim:

  Reduced-dimension name used when extracting embeddings from
  Bioconductor or Seurat objects.

- embedding:

  Optional data frame with `cell`, `UMAP1`, and `UMAP2`.

- assay_name:

  Optional assay name for `SingleCellExperiment` or `Seurat`.

- expression_layer:

  Optional Seurat layer or slot name.

- signatures:

  Optional named list of custom signatures.

- marker_panels:

  Optional named list of custom marker panels.

- lr_pairs:

  Optional data frame with `ligand`, `receptor`, and `pathway`.

## Value

An object of class `scdown_obj`.

## Details

- a long table with `cell`, `gene`, and `value`

- a list with `expr`, `cells`, and optional `embedding`

- a `SingleCellExperiment`

- a `Seurat` object
