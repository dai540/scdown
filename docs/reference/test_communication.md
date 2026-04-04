# Test communication scores by permutation

Test communication scores by permutation

## Usage

``` r
test_communication(
  obj,
  lr_pairs = NULL,
  n_perm = 50L,
  seed = 1L,
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

- n_perm:

  Number of label permutations.

- seed:

  Random seed.

- min_detected:

  Minimum detection rate.

- min_specificity:

  Minimum specificity score.

- exclude_pattern:

  Pattern for labels to exclude.

- outdir:

  Optional output directory.

## Value

A list with communication test tables and plot paths.

## Details

This is an exploratory permutation layer over the built-in communication
score. Use a specialist framework when communication is a central claim.
