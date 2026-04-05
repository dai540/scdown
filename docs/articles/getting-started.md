# Getting Started

`scdown` stays intentionally short:

1.  build one `scdown` object
2.  run quick exploratory summaries
3.  run sample-aware tests only when needed

## Build an object

``` r
obj <- scdown(scdown_example())
obj
#> <scdown_obj>
#>   cells: 48 
#>   genes: 15 
#>   samples: 4 
#>   annotations: 6 
#>   value type: counts
```

## Explore

``` r
map_res <- plot_map(obj)
marker_res <- explore_markers(obj, top_n = 3)
gene_res <- plot_gene(obj, genes = c("CD3D", "NKG7", "LYZ"))
sig_res <- plot_signature(obj, signatures = c("cytotoxic", "b_cell"))

kable(head(map_res$tables$embedding_cells, 6))
```

| cell   | UMAP1 | UMAP2 | sample | cell_type | group   | lineage |
|:-------|------:|------:|:-------|:----------|:--------|:--------|
| cell_1 | -2.43 |   2.1 | S1     | CD4_T     | control | T/NK    |
| cell_2 | -2.07 |   1.9 | S2     | CD4_T     | control | T/NK    |
| cell_3 | -2.16 |   2.1 | S3     | CD4_T     | case    | T/NK    |
| cell_4 | -1.90 |   1.9 | S4     | CD4_T     | case    | T/NK    |
| cell_5 | -1.73 |   1.5 | S1     | CD8_T     | control | T/NK    |
| cell_6 | -1.37 |   1.3 | S2     | CD8_T     | control | T/NK    |

``` r
kable(marker_res$tables$top_marker_per_celltype)
```

| cell_type | gene | mean_in_type | mean_rest | pct_in_type | pct_rest | marker_score | marker_table_rank |
|:---|:---|---:|---:|---:|---:|---:|---:|
| B_cell | MS4A1 | 7.964999 | 5.527988 | 1 | 1 | 2.437011 | 2.437011 |
| B_cell | CD79A | 7.893474 | 5.528006 | 1 | 1 | 2.365468 | 2.365468 |
| CD4_T | IL7R | 7.765374 | 5.553538 | 1 | 1 | 2.211836 | 2.211836 |
| CD4_T | CD3D | 7.765375 | 5.987036 | 1 | 1 | 1.778339 | 1.778339 |
| CD4_T | CCR5 | 7.420032 | 6.282412 | 1 | 1 | 1.137620 | 1.137620 |
| CD8_T | CCL5 | 7.862159 | 5.695226 | 1 | 1 | 2.166933 | 2.166933 |
| CD8_T | CD3D | 7.732316 | 5.993648 | 1 | 1 | 1.738669 | 1.738669 |
| CD8_T | CCR5 | 7.386979 | 6.289023 | 1 | 1 | 1.097956 | 1.097956 |
| Macro | C1QC | 7.675544 | 5.571430 | 1 | 1 | 2.104114 | 2.104114 |
| Macro | APOE | 7.675544 | 5.571436 | 1 | 1 | 2.104108 | 2.104108 |
| Macro | HLA-DRA | 7.779318 | 5.709578 | 1 | 1 | 2.069739 | 2.069739 |
| Mono | CXCL8 | 7.571326 | 5.592191 | 1 | 1 | 1.979135 | 1.979135 |
| Mono | FCN1 | 7.571327 | 5.592244 | 1 | 1 | 1.979083 | 1.979083 |
| Mono | LYZ | 7.571327 | 5.592245 | 1 | 1 | 1.979081 | 1.979081 |
| NK | NKG7 | 7.732318 | 5.560177 | 1 | 1 | 2.172141 | 2.172141 |
| NK | GNLY | 7.732318 | 5.560178 | 1 | 1 | 2.172140 | 2.172140 |
| NK | CCR5 | 7.386979 | 6.289023 | 1 | 1 | 1.097956 | 1.097956 |

``` r
kable(gene_res$tables$gene_average_by_celltype)
```

| gene | cell_type | mean_expression |
|:-----|:----------|----------------:|
| CD3D | B_cell    |        5.725435 |
| LYZ  | B_cell    |        5.725385 |
| NKG7 | B_cell    |        5.725385 |
| CD3D | CD4_T     |        7.765375 |
| LYZ  | CD4_T     |        5.597757 |
| NKG7 | CD4_T     |        5.597757 |
| CD3D | CD8_T     |        7.732316 |
| LYZ  | CD8_T     |        5.564884 |
| NKG7 | CD8_T     |        5.564884 |
| CD3D | Macro     |        5.508223 |
| LYZ  | Macro     |        5.508318 |
| NKG7 | Macro     |        5.508318 |
| CD3D | Mono      |        5.404446 |
| LYZ  | Mono      |        7.571327 |
| NKG7 | Mono      |        5.404541 |
| CD3D | NK        |        5.564760 |
| LYZ  | NK        |        5.564884 |
| NKG7 | NK        |        7.732318 |

``` r
kable(sig_res$tables$signature_score_by_celltype)
```

| signature | cell_type |     score |
|:----------|:----------|----------:|
| b_cell    | B_cell    | 0.9642857 |
| cytotoxic | B_cell    | 0.4330357 |
| b_cell    | CD4_T     | 0.3839286 |
| cytotoxic | CD4_T     | 0.4017857 |
| b_cell    | CD8_T     | 0.3750000 |
| cytotoxic | CD8_T     | 0.4107143 |
| b_cell    | Macro     | 0.3928571 |
| cytotoxic | Macro     | 0.3928571 |
| b_cell    | Mono      | 0.3303571 |
| cytotoxic | Mono      | 0.3303571 |
| b_cell    | NK        | 0.3571429 |
| cytotoxic | NK        | 0.9642857 |

## Test

``` r
comp_test <- test_composition(obj)
sig_test <- test_signature(obj, signatures = c("cytotoxic", "b_cell"))

kable(comp_test$tables$composition_test_table)
```

|        | feature | method | p_value | effect | note                        | p_adj |
|:-------|:--------|:-------|--------:|-------:|:----------------------------|------:|
| B_cell | B_cell  | wilcox |       1 |      0 | no variation between groups |     1 |
| CD4_T  | CD4_T   | wilcox |       1 |      0 | no variation between groups |     1 |
| CD8_T  | CD8_T   | wilcox |       1 |      0 | no variation between groups |     1 |
| Macro  | Macro   | wilcox |       1 |      0 | no variation between groups |     1 |
| Mono   | Mono    | wilcox |       1 |      0 | no variation between groups |     1 |
| NK     | NK      | wilcox |       1 |      0 | no variation between groups |     1 |

``` r
kable(sig_test$tables$signature_test_table)
```

|  | signature | cell_type | method | p_value | effect | note | p_adj |
|:---|:---|:---|:---|---:|---:|:---|---:|
| b_cell.B_cell | b_cell | B_cell | wilcox | 1.0000000 | 0.0000000 | no variation between groups | 1.0000000 |
| cytotoxic.B_cell | cytotoxic | B_cell | wilcox | 0.2452781 | 0.2589286 |  | 0.3679172 |
| b_cell.CD4_T | b_cell | CD4_T | wilcox | 0.2452781 | 0.3392857 |  | 0.3679172 |
| cytotoxic.CD4_T | cytotoxic | CD4_T | wilcox | 0.2452781 | 0.1964286 |  | 0.3679172 |
| b_cell.CD8_T | b_cell | CD8_T | wilcox | 0.2452781 | 0.3214286 |  | 0.3679172 |
| cytotoxic.CD8_T | cytotoxic | CD8_T | wilcox | 0.6985354 | 0.2142857 |  | 0.8382424 |
| b_cell.Macro | b_cell | Macro | wilcox | 0.2452781 | 0.2678571 |  | 0.3679172 |
| cytotoxic.Macro | cytotoxic | Macro | wilcox | 0.2452781 | 0.2678571 |  | 0.3679172 |
| b_cell.Mono | b_cell | Mono | wilcox | 0.6985354 | 0.1071429 |  | 0.8382424 |
| cytotoxic.Mono | cytotoxic | Mono | wilcox | 0.2452781 | 0.3750000 |  | 0.3679172 |
| b_cell.NK | b_cell | NK | wilcox | 0.2452781 | 0.1964286 |  | 0.3679172 |
| cytotoxic.NK | cytotoxic | NK | wilcox | 1.0000000 | 0.0000000 | no variation between groups | 1.0000000 |

## Custom resources

``` r
obj_custom <- scdown(
  scdown_example(),
  signatures = list(my_state = c("CD3D", "IL7R")),
  marker_panels = list(my_t_panel = c("CD3D", "TRAC"))
)

kable(plot_signature(obj_custom, signatures = "my_state")$tables$signature_score_by_celltype)
```

| signature | cell_type |     score |
|:----------|:----------|----------:|
| my_state  | B_cell    | 0.4241071 |
| my_state  | CD4_T     | 0.9642857 |
| my_state  | CD8_T     | 0.6562500 |
| my_state  | Macro     | 0.3928571 |
| my_state  | Mono      | 0.3303571 |
| my_state  | NK        | 0.3571429 |

## Report

``` r
report_path <- build_scdown_report(obj)
report_path
#> [1] "C:\\Users\\daiki\\AppData\\Local\\Temp\\RtmpSg9BgS\\scdown-report-269c26082377/scdown_report.html"
```
