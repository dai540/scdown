test_that("marker discovery recovers expected lineage genes", {
  obj <- scdown(scdown_example(), annotation_col = "annotation", sample_col = "sample")
  res <- find_markers(obj, top_n = 5)
  top <- res$tables$top_markers

  nk_markers <- top$gene[top$annotation == "NK"]
  b_markers <- top$gene[top$annotation == "B_cell"]
  macro_markers <- top$gene[top$annotation == "Macrophage"]

  expect_true(any(nk_markers %in% c("NKG7", "GNLY", "CCL5")))
  expect_true(any(b_markers %in% c("MS4A1", "CD79A", "CD74")))
  expect_true(any(macro_markers %in% c("C1QC", "APOE", "HLA-DRA", "FCER1A")))
})

test_that("core downstream functions return expected tables", {
  obj <- scdown(
    scdown_example(),
    annotation_col = "annotation",
    sample_col = "sample",
    signatures = list(my_state = c("CD3D", "IL7R"))
  )

  map_res <- plot_map(obj)
  expect_true(all(c("cell", "UMAP1", "UMAP2") %in% names(map_res$tables$embedding_cells)))

  check_res <- check_annotation(obj)
  expect_true(nrow(check_res$tables$annotation_panel_scores) > 0)

  avg_res <- plot_average_expression(obj)
  expect_true(nrow(avg_res$tables$average_expression) > 0)

  sim_res <- plot_cluster_similarity(obj)
  expect_true(nrow(sim_res$tables$cluster_similarity) > 0)

  sig_res <- plot_signature(obj, signatures = c("cytotoxic", "my_state"))
  expect_true(all(c("annotation", "score_name", "score") %in% names(sig_res$tables$signature_score_by_annotation)))

  comm_res <- infer_communication(obj)
  expect_true(all(c("sender", "receiver", "score") %in% names(comm_res$tables$communication_pair_summary)))
})

test_that("run_core returns the minimal downstream sections", {
  obj <- scdown(scdown_example(), annotation_col = "annotation", sample_col = "sample")
  res <- run_core(obj)
  expect_named(
    res,
    c("summary", "map", "markers", "annotation_check", "average_expression", "cluster_similarity", "signatures", "communication")
  )
})
