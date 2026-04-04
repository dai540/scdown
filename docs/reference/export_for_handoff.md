# Export `scdown` data for specialist handoff

This helper does not run `muscat`, `dreamlet`, `miloR`, or CellChat for
you. Instead, it prepares a compact handoff object with the matrices,
metadata, and embeddings that those downstream specialist tools usually
need.

## Usage

``` r
export_for_handoff(
  obj,
  target = c("muscat", "dreamlet", "miloR", "CellChat"),
  assay = c("counts", "analysis"),
  include_lineage = TRUE
)
```

## Arguments

- obj:

  A `scdown_obj`.

- target:

  One of `"muscat"`, `"dreamlet"`, `"miloR"`, or `"CellChat"`.

- assay:

  Which matrix to emphasize in the handoff: `"counts"` or `"analysis"`.

- include_lineage:

  Whether to ensure a lineage column is included in the exported cell
  metadata.

## Value

A named list with matrices, metadata, guidance, and, when the needed
packages are installed, a `SingleCellExperiment` object for Bioconductor
handoff targets.
