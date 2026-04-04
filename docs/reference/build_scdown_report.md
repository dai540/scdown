# Build a compact HTML report

Build a compact HTML report

## Usage

``` r
build_scdown_report(
  obj,
  genes = c("CD3D", "NKG7", "LYZ", "MS4A1"),
  signatures = c("cytotoxic", "exhaustion", "antigen_presentation"),
  signature_sets = NULL,
  run_tests = TRUE,
  scenario = c("generic", "pbmc", "tumor_immune", "patient_comparison"),
  outdir = tempfile("scdown-report-")
)
```

## Arguments

- obj:

  A `scdown_obj`.

- genes:

  User-selected genes.

- signatures:

  User-selected signatures.

- signature_sets:

  Optional custom signatures.

- run_tests:

  Whether to include lightweight tests and, when possible, sample-aware
  tests.

- scenario:

  One of `"generic"`, `"pbmc"`, `"tumor_immune"`, or
  `"patient_comparison"`. Used for the scope and handoff section.

- outdir:

  Output directory.

## Value

Path to the generated HTML file.
