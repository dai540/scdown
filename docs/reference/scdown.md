# Create a `scdown` object

Create a `scdown` object

## Usage

``` r
scdown(
  x,
  annotation_col = "cell_type",
  sample_col = "sample",
  group_col = "group",
  batch_col = NULL,
  reduced_dim = "UMAP",
  embedding = NULL,
  assay_name = NULL,
  expression_layer = NULL,
  signatures = NULL,
  marker_panels = NULL,
  lineage_mapping = NULL,
  lr_pairs = NULL
)
```

## Arguments

- x:

  A long-format expression table, a list with `expr`, optional `cells`,
  and optional `embedding`, or a `SingleCellExperiment` / `Seurat`
  object.

- annotation_col:

  Cell annotation column name.

- sample_col:

  Sample column name.

- group_col:

  Optional group column name.

- batch_col:

  Optional batch column name.

- reduced_dim:

  Name of the reduced dimension to use.

- embedding:

  Optional embedding data frame with `cell`, `UMAP1`, `UMAP2`.

- assay_name:

  Optional assay name for `SingleCellExperiment` or `Seurat`.

- expression_layer:

  Optional Seurat layer/slot name.

- signatures:

  Optional named list of custom signatures.

- marker_panels:

  Optional named list of custom marker panels.

- lineage_mapping:

  Optional data frame with `pattern` and `lineage`.

- lr_pairs:

  Optional data frame with `ligand`, `receptor`, and `pathway`.

## Value

An object of class `scdown_obj`.
