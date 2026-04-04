# Explore cell-cell communication

Explore cell-cell communication

## Usage

``` r
explore_communication(
  obj,
  lr_pairs = NULL,
  min_detected = 0.1,
  min_specificity = 0.15,
  exclude_pattern = "(?i)unassigned|unknown|doublet",
  outdir = NULL
)
```

## Arguments

- obj:

  A `scdown_obj`.

- lr_pairs:

  Optional data frame with custom ligand-receptor pairs.

- min_detected:

  Minimum detection rate for ligand and receptor.

- min_specificity:

  Minimum positive enrichment over the rest.

- exclude_pattern:

  Pattern for labels to exclude from communication.

- outdir:

  Optional output directory.

## Value

A list with communication tables and plots.
