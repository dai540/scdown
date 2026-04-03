test_that("constructor returns scdown object", {
  obj <- scdown(scdown_example())
  expect_s3_class(obj, "scdown_obj")
  expect_true(all(c("expr", "cells", "embedding") %in% names(obj)))
})

test_that("lineage collapse adds lineage", {
  obj <- collapse_lineage(scdown(scdown_example()))
  expect_true("lineage" %in% names(obj$cells))
})

test_that("built-ins are available", {
  expect_true("cytotoxic" %in% list_signatures())
  expect_true("NKG7" %in% get_signature("cytotoxic"))
})

test_that("core analysis functions return expected tables", {
  obj <- collapse_lineage(scdown(scdown_example()))
  expect_true("embedding_cells" %in% names(plot_map(obj)$tables))
  expect_true("top_marker_per_celltype" %in% names(find_markers(obj)$tables))
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

test_that("report builder creates html file", {
  obj <- scdown(scdown_example())
  path <- build_scdown_report(obj, outdir = tempfile("scdown-report-"))
  expect_true(file.exists(path))
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
    assays = list(logcounts = mat),
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
