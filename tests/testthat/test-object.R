to_long_expr <- function(x) {
  expr <- x$expr
  if (inherits(expr, "sparseMatrix")) {
    triplet <- Matrix::summary(expr)
    data.frame(
      cell = colnames(expr)[triplet$j],
      gene = rownames(expr)[triplet$i],
      value = triplet$x,
      annotation = x$cells$annotation[match(colnames(expr)[triplet$j], x$cells$cell)],
      sample = x$cells$sample[match(colnames(expr)[triplet$j], x$cells$cell)],
      stringsAsFactors = FALSE
    )
  } else {
    triplet <- as.data.frame(as.table(as.matrix(expr)), stringsAsFactors = FALSE)
    names(triplet) <- c("gene", "cell", "value")
    triplet$annotation <- x$cells$annotation[match(triplet$cell, x$cells$cell)]
    triplet$sample <- x$cells$sample[match(triplet$cell, x$cells$cell)]
    triplet
  }
}

test_that("constructor accepts example list input", {
  obj <- scdown(scdown_example(), annotation_col = "annotation", sample_col = "sample")
  expect_s3_class(obj, "scdown_obj")
  expect_equal(ncol(obj$expr), nrow(obj$cells))
  expect_equal(nrow(obj$embedding), nrow(obj$cells))
  expect_identical(obj$value_type, "counts")
})

test_that("constructor accepts long-format input and fills sample when absent", {
  x <- scdown_example()
  long_df <- to_long_expr(x)
  long_df$sample <- NULL
  obj <- scdown(long_df, annotation_col = "annotation", sample_col = "sample")
  expect_s3_class(obj, "scdown_obj")
  expect_true(all(obj$cells$sample == "sample_1"))
})

test_that("constructor accepts SingleCellExperiment when available", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")
  x <- scdown_example()
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = x$expr),
    colData = S4Vectors::DataFrame(x$cells)
  )
  SingleCellExperiment::reducedDim(sce, "UMAP") <- as.matrix(x$embedding[, c("UMAP1", "UMAP2")])
  obj <- scdown(sce, annotation_col = "annotation", sample_col = "sample", reduced_dim = "UMAP")
  expect_s3_class(obj, "scdown_obj")
  expect_equal(nrow(obj$embedding), ncol(x$expr))
})

test_that("constructor accepts Seurat when available", {
  skip_if_not_installed("SeuratObject")
  x <- scdown_example()
  seu <- SeuratObject::CreateSeuratObject(counts = x$expr, meta.data = x$cells)
  emb <- as.matrix(x$embedding[, c("UMAP1", "UMAP2")])
  rownames(emb) <- x$embedding$cell
  seu[["UMAP"]] <- SeuratObject::CreateDimReducObject(
    embeddings = emb,
    key = "UMAP_",
    assay = SeuratObject::DefaultAssay(seu)
  )
  obj <- scdown(seu, annotation_col = "annotation", sample_col = "sample", reduced_dim = "UMAP")
  expect_s3_class(obj, "scdown_obj")
  expect_equal(nrow(obj$embedding), ncol(x$expr))
})
