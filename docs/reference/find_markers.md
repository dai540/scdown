# Find cell-type markers

Find cell-type markers

## Usage

``` r
find_markers(
  obj,
  top_n = 10L,
  exclude_pattern = "^(RPL|RPS|MT-|MALAT1$|MTRNR|HSP)",
  min_pct_in_type = 0.1,
  min_pct_diff = 0.05,
  min_marker_score = 0.25,
  outdir = NULL
)
```

## Arguments

- obj:

  A `scdown_obj`.

- top_n:

  Number of markers per cell type.

- exclude_pattern:

  Pattern for low-information genes to exclude.

- min_pct_in_type:

  Minimum detection rate inside the target cell type.

- min_pct_diff:

  Minimum detection-rate improvement over the rest.

- min_marker_score:

  Minimum mean-expression improvement over the rest.

- outdir:

  Optional output directory.

## Value

A list with marker tables and plot paths.
