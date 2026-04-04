# scdown

`scdown` is a simple downstream toolkit for annotated single-cell RNA-seq data.

It is best used **after** QC, normalization, embedding, clustering, and cell
annotation. In other words, `scdown` is a downstream hub, not a full upstream
single-cell framework.

## What `scdown` is for

- make one downstream analysis object
- run fast exploratory summaries
- run sample-aware tests only when you need them
- export a compact HTML report for review and collaboration

## What `scdown` is not for

- droplet QC or ambient RNA correction
- normalization, HVG selection, integration, clustering, or trajectory fitting
- publication-grade differential state as a primary specialist method

Those parts are better handled by upstream or specialist tools such as Seurat,
`scater`, `scuttle`, `scran`, `SingleR`, `muscat`, `dreamlet`, `miloR`, and
CellChat.

Internally, `scdown` keeps data in a sparse-first matrix form, so it no longer
needs to expand everything into a huge long table before analysis.

## Website

- Home: [https://dai540.github.io/scdown/](https://dai540.github.io/scdown/)
- Getting started: [https://dai540.github.io/scdown/articles/getting-started.html](https://dai540.github.io/scdown/articles/getting-started.html)
- Recommended pipeline: [https://dai540.github.io/scdown/articles/recommended-pipeline.html](https://dai540.github.io/scdown/articles/recommended-pipeline.html)
- Reference: [https://dai540.github.io/scdown/reference/index.html](https://dai540.github.io/scdown/reference/index.html)

## Installation

```r
install.packages("pak")
pak::pak("dai540/scdown")
```

## Example

```r
library(scdown)

obj <- scdown(scdown_example())

run_core(obj)
build_scdown_report(obj)
```

## Explore

```r
plot_map(obj)
explore_markers(obj)
plot_gene(obj, genes = c("CD3D", "NKG7", "LYZ"))
plot_composition(obj)
plot_signature(obj, signatures = c("cytotoxic", "exhaustion"))
infer_communication(obj)
```

## Test

```r
test_markers(obj)
test_composition(obj)
test_signature(obj, signatures = c("cytotoxic", "b_cell"))
test_communication(obj, n_perm = 25)
```

## Recommended positioning

```r
scope_table("generic")
recommended_pipeline("generic")
recommended_pipeline("pbmc")
recommended_pipeline("tumor_immune")
recommended_pipeline("patient_comparison")
```

## Specialist handoff

If an endpoint becomes central enough for publication-grade modeling, prepare a
handoff object instead of stretching `scdown` beyond its intended scope.

```r
muscat_input <- export_for_handoff(obj, target = "muscat")
dreamlet_input <- export_for_handoff(obj, target = "dreamlet")
milor_input <- export_for_handoff(obj, target = "miloR")
cellchat_input <- export_for_handoff(obj, target = "CellChat")
```

`build_scdown_report()` also includes a scope and handoff section, plus marker
and communication tests when `run_tests = TRUE`.

## Custom resources

```r
obj_custom <- scdown(
  scdown_example(),
  signatures = list(my_state = c("CD3D", "IL7R")),
  marker_panels = list(my_t_panel = c("CD3D", "TRAC")),
  lineage_mapping = data.frame(
    pattern = c("cd4", "cd8", "nk", "b", "mono", "macro"),
    lineage = c("T/NK", "T/NK", "T/NK", "B/Plasma", "Myeloid", "Myeloid")
  ),
  lr_pairs = data.frame(
    ligand = "CCL5",
    receptor = "CCR5",
    pathway = "chemokine"
  )
)

plot_signature(obj_custom, signatures = "my_state")
check_annotation(obj_custom, panels = "my_t_panel")
```

## Direct input

```r
obj_sce <- scdown(
  sce,
  annotation_col = "cell_type",
  sample_col = "sample",
  reduced_dim = "UMAP",
  assay_name = "counts"
)

obj_seurat <- scdown(
  seu,
  annotation_col = "cell_type",
  sample_col = "sample",
  reduced_dim = "umap",
  assay_name = "RNA",
  expression_layer = "counts"
)
```

## Main functions

- `scdown()`
- `recommended_pipeline()`
- `scope_table()`
- `export_for_handoff()`
- `summarize_dataset()`
- `collapse_lineage()`
- `plot_map()`
- `explore_markers()`
- `find_markers()`
- `plot_gene()`
- `plot_composition()`
- `explore_composition()`
- `check_annotation()`
- `plot_signature()`
- `explore_signature()`
- `plot_heatmap()`
- `infer_communication()`
- `explore_communication()`
- `test_markers()`
- `test_composition()`
- `test_signature()`
- `test_communication()`
- `run_core()`
- `build_scdown_report()`
