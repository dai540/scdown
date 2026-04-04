# Test composition differences across groups

Test composition differences across groups

## Usage

``` r
test_composition(obj, level = c("cell_type", "lineage"), outdir = NULL)
```

## Arguments

- obj:

  A `scdown_obj`.

- level:

  Either `"cell_type"` or `"lineage"`.

- outdir:

  Optional output directory.

## Value

A list with composition test tables and plot paths.

## Details

This is a lightweight sample-level fraction test. It is useful for quick
screening, but specialist composition or abundance models should take
over when composition is a primary result.
