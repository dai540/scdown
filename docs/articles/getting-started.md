# Getting Started

`scdown` is for annotation-after downstream analysis inside one dataset.

The package assumes that clustering and annotation are already done. The
core questions are:

- what clusters are present on the embedding?
- which genes define each cluster?
- do known marker panels support the annotation?
- which clusters are transcriptionally similar?
- which signatures or pathways are active?
- which ligand-receptor pairs are plausible?

## Example data

``` r
library(scdown)

example_data <- scdown_example()
str(example_data, max.level = 1)
#> List of 3
#>  $ expr     :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>  $ cells    :'data.frame':   108 obs. of  3 variables:
#>  $ embedding:'data.frame':   108 obs. of  3 variables:
```

The example object contains:

- `cells`: one row per cell with `cell`, `sample`, and `annotation`
- `embedding`: two-dimensional coordinates
- `expr`: a sparse count matrix

## Create a `scdown` object

``` r
obj <- scdown(
  example_data,
  annotation_col = "annotation",
  sample_col = "sample"
)

obj
#> <scdown_obj>
#>   cells: 108 
#>   genes: 21 
#>   annotations: 6 
#>   samples: 3 
#>   value type: counts
```

## Run the minimal downstream workflow

``` r
res <- run_core(
  obj,
  signatures = c("cytotoxic", "antigen_presentation")
)

names(res)
#> [1] "summary"            "map"                "markers"           
#> [4] "annotation_check"   "average_expression" "cluster_similarity"
#> [7] "signatures"         "communication"
```

## Marker discovery

``` r
head(res$markers$tables$top_markers)
#>   annotation gene  mean_in mean_rest pct_in pct_rest marker_score marker_rank
#> 1     B_cell CD74 7.601402  0.000000      1      0.0     7.601402    8.101402
#> 2      CD4_T CD3D 7.552637  1.510527      1      0.2     6.042110    6.442110
#> 3      CD4_T CCR5 6.860015  2.744006      1      0.4     4.116009    4.416009
#> 4      CD8_T PRF1 7.552637  0.000000      1      0.0     7.552637    8.052637
#> 5      CD8_T CD3D 7.552637  1.510527      1      0.2     6.042110    6.442110
#> 6      CD8_T CCL5 7.552637  1.510527      1      0.2     6.042110    6.442110
```

## Marker-based annotation check

``` r
head(res$annotation_check$tables$annotation_panel_scores)
#>   annotation score_name  score
#> 1     B_cell     b_cell 0.9750
#> 2      CD4_T     b_cell 0.5250
#> 3      CD8_T     b_cell 0.4750
#> 4 Macrophage     b_cell 0.5625
#> 5   Monocyte     b_cell 0.4500
#> 6         NK     b_cell 0.5250
```

## Average expression and similarity

``` r
head(res$average_expression$tables$average_expression)
#>   annotation gene mean_expression
#> 1      CD4_T CD74        0.000000
#> 2      CD4_T CD3D        7.552637
#> 3      CD4_T CCR5        6.860015
#> 4      CD4_T PRF1        0.000000
#> 5      CD4_T CCL5        0.000000
#> 6      CD4_T APOE        0.000000
head(res$cluster_similarity$tables$cluster_similarity)
#>   annotation_1 annotation_2 correlation
#> 1       B_cell       B_cell  1.00000000
#> 2        CD4_T       B_cell  0.41803263
#> 3        CD8_T       B_cell -0.03177355
#> 4   Macrophage       B_cell  0.60050246
#> 5     Monocyte       B_cell  0.57845072
#> 6           NK       B_cell  0.41803263
```

## Signature scoring

``` r
head(res$signatures$tables$signature_score_by_annotation)
#>   annotation           score_name  score
#> 1     B_cell antigen_presentation 0.8750
#> 2      CD4_T antigen_presentation 0.2375
#> 3      CD8_T antigen_presentation 0.1750
#> 4 Macrophage antigen_presentation 0.5625
#> 5   Monocyte antigen_presentation 0.4125
#> 6         NK antigen_presentation 0.2375
```

## Exploratory communication

``` r
head(res$communication$tables$communication_pair_summary)
#>     sender receiver     score
#> 1    CD8_T    CD4_T 24.869376
#> 2 Monocyte    CD4_T  8.870093
#> 3       NK    CD4_T 24.869376
#> 4    CD8_T    CD8_T 24.869376
#> 5 Monocyte    CD8_T  8.870093
#> 6       NK    CD8_T 24.869376
```

## Build a compact report

``` r
report_path <- build_scdown_report(obj)
report_path
#> [1] "C:\\Users\\daiki\\AppData\\Local\\Temp\\RtmpoVBCfw\\scdown-report-3d8c96931c9/scdown_report.html"
```

The report writes one HTML file plus CSV and PNG artifacts for each
section.
