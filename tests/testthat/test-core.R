test_that("constructor returns sparse-first scdown object", {
  obj <- scdown(scdown_example())
  expect_s3_class(obj, "scdown_obj")
  expect_true(all(c("expr", "cells", "embedding", "resources") %in% names(obj)))
  expect_equal(ncol(obj$expr), nrow(obj$cells))
  expect_equal(rownames(obj$expr), obj$features$gene)
})

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

test_that("built-ins and custom resources are available", {
  obj <- scdown(
    scdown_example(),
    signatures = list(my_state = c("CD3D", "IL7R")),
    marker_panels = list(my_panel = c("CD3D", "TRAC"))
  )
  expect_true("cytotoxic" %in% list_signatures())
  expect_true("my_state" %in% list_signatures(obj))
  expect_true("NKG7" %in% get_signature("cytotoxic"))
  expect_true("CD3D" %in% get_signature("my_state", obj))
})

test_that("recommended pipeline helper returns expected scenarios", {
  tbl <- recommended_pipeline("generic")
  expect_true(all(c("stage", "goal", "recommended_tools", "scdown_role") %in% names(tbl)))
  expect_true(any(tbl$scdown_role == "core"))

  pbmc <- recommended_pipeline("pbmc")
  expect_true(any(grepl("scdown", pbmc$recommended_tools, fixed = TRUE)))

  scope <- scope_table("generic")
  expect_true(all(c("stage", "use_scdown", "specialist_tools", "note") %in% names(scope)))
  expect_true(any(scope$use_scdown == "handoff"))
})

test_that("handoff export returns target-ready structures", {
  obj <- collapse_lineage(scdown(scdown_example()))

  muscat_export <- export_for_handoff(obj, target = "muscat")
  expect_true(all(c("expr", "cells", "sample_metadata", "embedding", "guidance") %in% names(muscat_export)))
  expect_equal(muscat_export$cluster_col, obj$annotation_col)

  cellchat_export <- export_for_handoff(obj, target = "CellChat")
  expect_true("cellchat_input" %in% names(cellchat_export))
  expect_equal(cellchat_export$cellchat_input$group.by, obj$annotation_col)
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

test_that("sample-aware test functions return expected tables", {
  obj <- collapse_lineage(scdown(scdown_example()))
  expect_true("top_marker_tests" %in% names(test_markers(obj, top_n = 3)$tables))
  expect_true("composition_test_table" %in% names(test_composition(obj)$tables))
  expect_true("signature_test_table" %in% names(test_signature(obj, signatures = c("cytotoxic", "b_cell"))$tables))
  expect_true("communication_test_table" %in% names(test_communication(obj, n_perm = 5)$tables))
})

test_that("markers stay positive and specific in example data", {
  obj <- collapse_lineage(scdown(scdown_example()))
  top_markers <- find_markers(obj, top_n = 3)$tables$top_marker_per_celltype
  expect_true(all(top_markers$marker_score > 0))
  expect_true("MS4A1" %in% top_markers$gene[top_markers$cell_type == "B_cell"])
  expect_true("NKG7" %in% top_markers$gene[top_markers$cell_type == "NK"])
})

test_that("report builder creates html file", {
  obj <- scdown(scdown_example())
  path <- build_scdown_report(obj, outdir = tempfile("scdown-report-"))
  expect_true(file.exists(path))
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_match(text, "Scope and handoff")
  expect_match(text, "Markers")
  expect_match(text, "Communication")
})

test_that("constructor can accept long-format tables", {
  demo <- scdown_example()
  long_expr <- .matrix_to_long_subset(
    raw_mat = demo$expr,
    analysis_mat = demo$expr,
    genes = rownames(demo$expr),
    cells = demo$cells,
    sample_col = "sample",
    annotation_col = "cell_type"
  )
  long_expr$group <- rep(demo$cells$group, times = nrow(demo$expr))
  obj <- scdown(long_expr, annotation_col = "cell_type", sample_col = "sample", group_col = "group", embedding = demo$embedding)
  expect_s3_class(obj, "scdown_obj")
})

test_that("constructor can accept SingleCellExperiment when available", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("S4Vectors")

  mat <- matrix(
    c(5, 1, 0, 2, 3, 4),
    nrow = 3,
    dimnames = list(c("CD3D", "NKG7", "LYZ"), c("cell1", "cell2"))
  )
  coldata <- S4Vectors::DataFrame(
    sample = c("S1", "S1"),
    cell_type = c("T_cell", "Myeloid"),
    row.names = c("cell1", "cell2")
  )
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = mat),
    colData = coldata
  )
  SingleCellExperiment::reducedDim(sce, "UMAP") <- matrix(
    c(0, 0, 1, 1),
    ncol = 2,
    dimnames = list(c("cell1", "cell2"), c("UMAP1", "UMAP2"))
  )
  obj <- scdown(sce, annotation_col = "cell_type", sample_col = "sample")
  expect_s3_class(obj, "scdown_obj")
})

test_that("constructor can accept Seurat when available", {
  skip_if_not_installed("Matrix")
  skip_if_not_installed("SeuratObject")

  mat <- Matrix::Matrix(
    c(5, 1, 0, 2, 3, 4),
    nrow = 3,
    dimnames = list(c("CD3D", "NKG7", "LYZ"), c("cell1", "cell2")),
    sparse = TRUE
  )
  meta <- data.frame(
    sample = c("S1", "S1"),
    cell_type = c("T_cell", "Myeloid"),
    row.names = c("cell1", "cell2"),
    stringsAsFactors = FALSE
  )
  seu <- SeuratObject::CreateSeuratObject(counts = mat, assay = "RNA", meta.data = meta)
  emb <- matrix(c(0, 0, 1, 1), ncol = 2, dimnames = list(c("cell1", "cell2"), c("UMAP_1", "UMAP_2")))
  seu[["umap"]] <- SeuratObject::CreateDimReducObject(embeddings = emb, key = "UMAP_", assay = "RNA")
  obj <- scdown(seu, annotation_col = "cell_type", sample_col = "sample", reduced_dim = "umap", assay_name = "RNA", expression_layer = "counts")
  expect_s3_class(obj, "scdown_obj")
})
