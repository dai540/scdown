# Test signature differences across groups

Test signature differences across groups

## Usage

``` r
test_signature(obj, signatures, signature_sets = NULL, outdir = NULL)
```

## Arguments

- obj:

  A `scdown_obj`.

- signatures:

  Character vector of signature names or a named list of signatures.

- signature_sets:

  Optional named list of custom signatures.

- outdir:

  Optional output directory.

## Value

A list with signature test tables and plot paths.

## Details

This is a compact group comparison over sample-level mean signature
scores. It is intended for review and prioritization, not as a
replacement for more flexible specialist modeling.
