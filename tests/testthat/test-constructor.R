test_that("constructor returns sparse-first scdown object", {
  obj <- scdown(scdown_example())
  expect_s3_class(obj, "scdown_obj")
  expect_true(all(c("expr", "cells", "embedding", "resources") %in% names(obj)))
  expect_equal(ncol(obj$expr), nrow(obj$cells))
  expect_equal(rownames(obj$expr), obj$features$gene)
  expect_true(inherits(obj$expr, "sparseMatrix"))
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

test_that("sparse large input stays sparse through constructor and handoff", {
  demo <- example_sparse_large()
  obj <- scdown(demo)
  expect_true(inherits(obj$expr, "sparseMatrix"))

  muscat_export <- export_for_handoff(obj, target = "muscat")
  expect_true(inherits(muscat_export$expr, "sparseMatrix"))
  expect_equal(ncol(muscat_export$expr), nrow(obj$cells))
})
