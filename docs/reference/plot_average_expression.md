# Plot average expression by annotation

Plot average expression by annotation

## Usage

``` r
plot_average_expression(
  obj,
  features = NULL,
  top_n_markers = 5L,
  outdir = NULL
)
```

## Arguments

- obj:

  A `scdown_obj`.

- features:

  Optional character vector of genes to show.

- top_n_markers:

  Number of markers per annotation used when `features` is `NULL`.

- outdir:

  Optional output directory.

## Value

A list with average-expression tables and plot paths.
