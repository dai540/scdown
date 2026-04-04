# Run the core downstream workflow

Run the core downstream workflow

## Usage

``` r
run_core(
  obj,
  genes = c("CD3D", "NKG7", "LYZ", "MS4A1"),
  signatures = c("cytotoxic", "exhaustion", "antigen_presentation"),
  signature_sets = NULL,
  outdir = NULL
)
```

## Arguments

- obj:

  A `scdown_obj`.

- genes:

  User-selected genes.

- signatures:

  Signatures to score.

- signature_sets:

  Optional custom signatures.

- outdir:

  Optional output directory.

## Value

A named list of results.
