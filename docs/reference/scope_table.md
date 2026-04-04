# Scope table for `scdown`

`scdown` is strongest after annotation, when the goal is fast review,
lightweight testing, and compact reporting. This helper makes the
intended boundary explicit so users can see where `scdown` is enough and
where a specialist method should take over.

## Usage

``` r
scope_table(
  scenario = c("generic", "pbmc", "tumor_immune", "patient_comparison")
)
```

## Arguments

- scenario:

  One of `"generic"`, `"pbmc"`, `"tumor_immune"`, or
  `"patient_comparison"`.

## Value

A data frame with workflow stages, whether `scdown` is a good fit, and
which specialist tools are better suited when the analysis needs to go
further.
