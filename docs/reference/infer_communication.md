# Infer exploratory cell-cell communication

Infer exploratory cell-cell communication

## Usage

``` r
infer_communication(
  obj,
  lr_pairs = NULL,
  min_detected = 0.1,
  min_specificity = 0.15,
  outdir = NULL
)
```

## Arguments

- obj:

  A `scdown_obj`.

- lr_pairs:

  Optional ligand-receptor table.

- min_detected:

  Minimum detection rate for ligand and receptor.

- min_specificity:

  Minimum specificity score.

- outdir:

  Optional output directory.

## Value

A list with communication tables and plot paths.
