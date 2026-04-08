# Run the minimal downstream workflow

Run the minimal downstream workflow

## Usage

``` r
run_core(
  obj,
  signatures = c("cytotoxic", "antigen_presentation"),
  panels = c("t_cell", "nk", "b_cell", "monocyte", "macrophage", "dendritic"),
  signature_sets = NULL,
  marker_panels = NULL,
  outdir = NULL
)
```

## Arguments

- obj:

  A `scdown_obj`.

- signatures:

  Character vector of signatures to score.

- panels:

  Character vector of marker panels to use for annotation checks.

- signature_sets:

  Optional named list of custom signatures.

- marker_panels:

  Optional named list of custom marker panels.

- outdir:

  Optional output directory.

## Value

A named list of downstream results.
