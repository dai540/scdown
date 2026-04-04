# Test markers with lightweight per-cell Wilcoxon statistics

Test markers with lightweight per-cell Wilcoxon statistics

## Usage

``` r
test_markers(
  obj,
  top_n = 10L,
  exclude_pattern = "^(RPL|RPS|MT-|MALAT1$|MTRNR|HSP)",
  min_pct_in_type = 0.1,
  max_features = 2000L,
  outdir = NULL
)
```

## Arguments

- obj:

  A `scdown_obj`.

- top_n:

  Number of marker tests to keep per cell type.

- exclude_pattern:

  Pattern for low-information genes to exclude.

- min_pct_in_type:

  Minimum detection rate inside the target cell type.

- max_features:

  Maximum features to test after prefiltering.

- outdir:

  Optional output directory.

## Value

A list with marker test tables and plot paths.

## Details

This is a quick review-layer test. It is useful for prioritizing marker
candidates, but it is not a replacement for pseudobulk, mixed-model, or
covariate-aware specialist workflows.
