#' Create a `scdown` object
#'
#' @param x A long-format expression table, a list with `expr`, optional
#'   `cells`, and optional `embedding`, or a `SingleCellExperiment` / `Seurat`
#'   object.
#' @param annotation_col Cell annotation column name.
#' @param sample_col Sample column name.
#' @param group_col Optional group column name.
#' @param batch_col Optional batch column name.
#' @param reduced_dim Name of the reduced dimension to use.
#' @param embedding Optional embedding data frame with `cell`, `UMAP1`, `UMAP2`.
#' @param assay_name Optional assay name for `SingleCellExperiment` or `Seurat`.
#' @param expression_layer Optional Seurat layer/slot name.
#'
#' @return An object of class `scdown_obj`.
#' @export
scdown <- function(x,
                   annotation_col = "cell_type",
                   sample_col = "sample",
                   group_col = "group",
                   batch_col = NULL,
                   reduced_dim = "UMAP",
                   embedding = NULL,
                   assay_name = NULL,
                   expression_layer = NULL) {
  parsed <- .coerce_scdown_input(
    x = x,
    annotation_col = annotation_col,
    sample_col = sample_col,
    group_col = group_col,
    batch_col = batch_col,
    reduced_dim = reduced_dim,
    embedding = embedding,
    assay_name = assay_name,
    expression_layer = expression_layer
  )
  out <- list(
    expr = parsed$expr,
    cells = parsed$cells,
    embedding = parsed$embedding,
    annotation_col = annotation_col,
    sample_col = sample_col,
    group_col = group_col,
    batch_col = batch_col,
    reduced_dim = reduced_dim,
    lineage_col = "lineage"
  )
  class(out) <- "scdown_obj"
  out
}

#' @export
print.scdown_obj <- function(x, ...) {
  cat("<scdown_obj>\n")
  cat("  cells:", length(unique(x$cells$cell)), "\n")
  cat("  samples:", length(unique(x$cells[[x$sample_col]])), "\n")
  cat("  annotations:", length(unique(x$cells[[x$annotation_col]])), "\n")
  invisible(x)
}

.coerce_scdown_input <- function(x,
                                 annotation_col,
                                 sample_col,
                                 group_col,
                                 batch_col,
                                 reduced_dim,
                                 embedding,
                                 assay_name,
                                 expression_layer) {
  if (is.data.frame(x)) {
    expr <- x
    .assert_columns(expr, c(sample_col, "cell", annotation_col, "gene", "value"), "x")
    cells <- unique(expr[, .existing_columns(expr, c("cell", sample_col, annotation_col, group_col, batch_col)), drop = FALSE])
  } else if (is.list(x) && !is.null(x$expr)) {
    expr <- x$expr
    .assert_columns(expr, c(sample_col, "cell", annotation_col, "gene", "value"), "x$expr")
    cells <- if (!is.null(x$cells)) x$cells else unique(expr[, .existing_columns(expr, c("cell", sample_col, annotation_col, group_col, batch_col)), drop = FALSE])
    if (is.null(embedding) && !is.null(x$embedding)) {
      embedding <- x$embedding
    }
  } else if (inherits(x, "SingleCellExperiment")) {
    assay_name <- .resolve_sce_assay(x, assay_name)
    cells <- as.data.frame(SummarizedExperiment::colData(x), stringsAsFactors = FALSE)
    cells$cell <- colnames(x)
    .assert_columns(cells, c("cell", sample_col, annotation_col), "colData(x)")
    expr <- .matrix_to_long(
      mat = SummarizedExperiment::assay(x, assay_name),
      cells = cells,
      sample_col = sample_col,
      annotation_col = annotation_col,
      group_col = group_col,
      batch_col = batch_col
    )
    if (is.null(embedding)) {
      embedding <- .extract_sce_embedding(x, reduced_dim)
    }
  } else if (inherits(x, "Seurat")) {
    cells <- as.data.frame(x[[]], stringsAsFactors = FALSE)
    cells$cell <- rownames(cells)
    .assert_columns(cells, c("cell", sample_col, annotation_col), "x[[]]")
    expr <- .matrix_to_long(
      mat = .extract_seurat_matrix(x, assay_name = assay_name, expression_layer = expression_layer),
      cells = cells,
      sample_col = sample_col,
      annotation_col = annotation_col,
      group_col = group_col,
      batch_col = batch_col
    )
    if (is.null(embedding)) {
      embedding <- .extract_seurat_embedding(x, reduced_dim)
    }
  } else {
    stop("`x` must be a long-format data.frame, a list with `expr`, or a supported `SingleCellExperiment` / `Seurat` object.", call. = FALSE)
  }

  cells <- as.data.frame(cells, stringsAsFactors = FALSE)
  expr <- as.data.frame(expr, stringsAsFactors = FALSE)
  if (!is.null(embedding)) {
    .assert_columns(embedding, c("cell", "UMAP1", "UMAP2"), "embedding")
    embedding <- as.data.frame(embedding, stringsAsFactors = FALSE)
  } else {
    embedding <- .deterministic_embedding(unique(cells$cell))
  }
  list(expr = expr, cells = cells, embedding = embedding)
}

.resolve_sce_assay <- function(x, assay_name = NULL) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("`SummarizedExperiment` is required for `SingleCellExperiment` input.", call. = FALSE)
  }
  available <- names(SummarizedExperiment::assays(x))
  if (!length(available)) {
    stop("No assays found in `SingleCellExperiment` input.", call. = FALSE)
  }
  if (!is.null(assay_name)) {
    if (!assay_name %in% available) {
      stop("Requested assay was not found in `SingleCellExperiment` input.", call. = FALSE)
    }
    return(assay_name)
  }
  preferred <- c("logcounts", "counts")
  hit <- preferred[preferred %in% available]
  if (length(hit)) hit[[1]] else available[[1]]
}

.matrix_to_long <- function(mat, cells, sample_col, annotation_col, group_col, batch_col) {
  mat <- as.matrix(mat)
  if (is.null(rownames(mat))) {
    rownames(mat) <- paste0("gene_", seq_len(nrow(mat)))
  }
  if (is.null(colnames(mat))) {
    colnames(mat) <- cells$cell
  }
  long <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  names(long) <- c("gene", "cell", "value")
  meta_cols <- .existing_columns(cells, c("cell", sample_col, annotation_col, group_col, batch_col))
  out <- merge(cells[, meta_cols, drop = FALSE], long, by = "cell", all.y = TRUE, sort = FALSE)
  out[, c(.existing_columns(out, c(sample_col, "cell", annotation_col, group_col, batch_col)), "gene", "value")]
}

.extract_sce_embedding <- function(x, reduced_dim) {
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    return(.deterministic_embedding(colnames(x)))
  }
  emb <- tryCatch(SingleCellExperiment::reducedDim(x, reduced_dim), error = function(e) NULL)
  if (is.null(emb) || ncol(as.matrix(emb)) < 2L) {
    return(.deterministic_embedding(colnames(x)))
  }
  emb <- as.matrix(emb)
  data.frame(cell = rownames(emb) %||% colnames(x), UMAP1 = emb[, 1], UMAP2 = emb[, 2], stringsAsFactors = FALSE)
}

.extract_seurat_matrix <- function(x, assay_name = NULL, expression_layer = NULL) {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("`SeuratObject` is required for `Seurat` input.", call. = FALSE)
  }
  assay_name <- assay_name %||% SeuratObject::DefaultAssay(x)
  expression_layer <- expression_layer %||% "data"
  if (exists("LayerData", where = asNamespace("SeuratObject"), mode = "function")) {
    mat <- tryCatch(SeuratObject::LayerData(x, assay = assay_name, layer = expression_layer), error = function(e) NULL)
    if (!is.null(mat)) {
      return(mat)
    }
  }
  SeuratObject::GetAssayData(x, assay = assay_name, slot = expression_layer)
}

.extract_seurat_embedding <- function(x, reduced_dim) {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    return(.deterministic_embedding(rownames(x[[]])))
  }
  emb <- tryCatch(SeuratObject::Embeddings(x[[reduced_dim]]), error = function(e) NULL)
  if (is.null(emb) || ncol(as.matrix(emb)) < 2L) {
    return(.deterministic_embedding(rownames(x[[]])))
  }
  emb <- as.matrix(emb)
  data.frame(cell = rownames(emb), UMAP1 = emb[, 1], UMAP2 = emb[, 2], stringsAsFactors = FALSE)
}
