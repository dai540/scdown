# Plot expression heatmaps

Plot expression heatmaps

## Usage

``` r
plot_heatmap(
  obj,
  signatures = c("cytotoxic", "exhaustion", "antigen_presentation"),
  signature_sets = NULL,
  top_n_markers = 5L,
  outdir = NULL
)
```

## Arguments

- obj:

  A `scdown_obj`.

- signatures:

  Signatures to include in the signature heatmap.

- signature_sets:

  Optional custom signatures.

- top_n_markers:

  Number of markers per cell type to show.

- outdir:

  Optional output directory.

## Value

A list with heatmap tables and plots.
