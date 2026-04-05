test_that("lineage collapse uses defaults and custom mappings", {
  obj <- scdown(scdown_example())
  obj_default <- collapse_lineage(obj)
  expect_true("lineage" %in% names(obj_default$cells))

  custom_map <- data.frame(
    pattern = c("cd4", "cd8", "nk", "b", "mono", "macro"),
    lineage = c("L1", "L1", "L1", "L2", "L3", "L3"),
    stringsAsFactors = FALSE
  )
  obj_custom <- collapse_lineage(obj, mapping = custom_map)
  expect_true(all(obj_custom$cells$lineage %in% c("L1", "L2", "L3")))
})

test_that("explore functions return expected tables", {
  obj <- collapse_lineage(scdown(scdown_example()))
  expect_true("embedding_cells" %in% names(plot_map(obj)$tables))
  expect_true("top_marker_per_celltype" %in% names(explore_markers(obj)$tables))
  expect_true("gene_average_by_celltype" %in% names(plot_gene(obj, genes = c("CD3D", "LYZ"))$tables))
  expect_true("composition_by_sample" %in% names(plot_composition(obj)$tables))
  expect_true("annotation_check_score" %in% names(check_annotation(obj)$tables))
  expect_true("signature_score_by_celltype" %in% names(plot_signature(obj, signatures = "cytotoxic")$tables))
  expect_true("communication_pair_summary" %in% names(infer_communication(obj)$tables))
})

test_that("markers stay positive and specific in example data", {
  obj <- collapse_lineage(scdown(scdown_example()))
  top_markers <- find_markers(obj, top_n = 3)$tables$top_marker_per_celltype
  expect_true(all(top_markers$marker_score > 0))
  expect_true("MS4A1" %in% top_markers$gene[top_markers$cell_type == "B_cell"])
  expect_true("NKG7" %in% top_markers$gene[top_markers$cell_type == "NK"])
})

test_that("plot_gene records missing genes instead of failing", {
  obj <- scdown(scdown_example())
  res <- plot_gene(obj, genes = c("CD3D", "NOT_A_GENE"))
  expect_true("missing_genes" %in% names(res$tables))
  expect_true("NOT_A_GENE" %in% res$tables$missing_genes$gene)
})
