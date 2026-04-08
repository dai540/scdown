# Build a compact HTML report

Build a compact HTML report

## Usage

``` r
build_scdown_report(
  obj,
  signatures = c("cytotoxic", "antigen_presentation"),
  panels = c("t_cell", "nk", "b_cell", "monocyte", "macrophage", "dendritic"),
  signature_sets = NULL,
  marker_panels = NULL,
  outdir = tempfile("scdown-report-")
)
```

## Arguments

- obj:

  A `scdown_obj`.

- signatures:

  Character vector of built-in signature names or a named list of
  signatures.

- panels:

  Character vector of built-in marker panel names or a named list of
  marker panels.

- signature_sets:

  Optional named list of custom signatures.

- marker_panels:

  Optional named list of custom marker panels.

- outdir:

  Output directory.

## Value

Path to the generated HTML report.
