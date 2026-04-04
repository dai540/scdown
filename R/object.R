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
#' @param signatures Optional named list of custom signatures.
#' @param marker_panels Optional named list of custom marker panels.
#' @param lineage_mapping Optional data frame with `pattern` and `lineage`.
#' @param lr_pairs Optional data frame with `ligand`, `receptor`, and `pathway`.
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
                   expression_layer = NULL,
                   signatures = NULL,
                   marker_panels = NULL,
                   lineage_mapping = NULL,
                   lr_pairs = NULL) {
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
    features = data.frame(gene = rownames(parsed$expr), stringsAsFactors = FALSE),
    embedding = parsed$embedding,
    annotation_col = annotation_col,
    sample_col = sample_col,
    group_col = group_col,
    batch_col = batch_col,
    reduced_dim = reduced_dim,
    lineage_col = "lineage",
    value_type = if (.is_count_like_matrix(parsed$expr)) "counts" else "continuous",
    resources = list(
      signatures = .validate_named_gene_sets(signatures, "signatures"),
      marker_panels = .validate_named_gene_sets(marker_panels, "marker_panels"),
      lineage_mapping = .validate_lineage_mapping(lineage_mapping),
      lr_pairs = .validate_lr_pairs(lr_pairs)
    )
  )
  class(out) <- "scdown_obj"
  out
}

#' @export
print.scdown_obj <- function(x, ...) {
  cat("<scdown_obj>\n")
  cat("  cells:", ncol(x$expr), "\n")
  cat("  genes:", nrow(x$expr), "\n")
  cat("  samples:", length(unique(x$cells[[x$sample_col]])), "\n")
  cat("  annotations:", length(unique(x$cells[[x$annotation_col]])), "\n")
  cat("  value type:", x$value_type, "\n")
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
    .assert_columns(x, c(sample_col, "cell", annotation_col, "gene", "value"), "x")
    cells <- unique(x[, .existing_columns(x, c("cell", sample_col, annotation_col, group_col, batch_col)), drop = FALSE])
    expr <- .long_to_matrix(x, cells$cell)
  } else if (is.list(x) && !is.null(x$expr)) {
    if (is.data.frame(x$expr)) {
      .assert_columns(x$expr, c("cell", "gene", "value"), "x$expr")
      cells <- if (!is.null(x$cells)) {
        x$cells
      } else {
        unique(x$expr[, .existing_columns(x$expr, c("cell", sample_col, annotation_col, group_col, batch_col)), drop = FALSE])
      }
      cells <- .standardize_cells(cells, unique(x$expr$cell), sample_col, annotation_col, group_col, batch_col)
      expr <- .long_to_matrix(x$expr, cells$cell)
    } else if (is.matrix(x$expr) || .is_sparse_matrix(x$expr)) {
      if (is.null(x$cells)) {
        stop("List inputs with matrix-like `expr` must also supply `cells` metadata.", call. = FALSE)
      }
      expr <- x$expr
      if (is.null(rownames(expr))) {
        rownames(expr) <- paste0("gene_", seq_len(nrow(expr)))
      }
      if (is.null(colnames(expr))) {
        stop("Matrix-like `expr` inputs must have column names.", call. = FALSE)
      }
      cells <- .standardize_cells(x$cells, colnames(expr), sample_col, annotation_col, group_col, batch_col)
    } else {
      stop("Unsupported `x$expr` format.", call. = FALSE)
    }
    if (is.null(embedding) && !is.null(x$embedding)) {
      embedding <- x$embedding
    }
  } else if (inherits(x, "SingleCellExperiment")) {
    assay_name <- .resolve_sce_assay(x, assay_name)
    expr <- SummarizedExperiment::assay(x, assay_name)
    if (is.null(rownames(expr))) {
      rownames(expr) <- rownames(x)
    }
    cells <- as.data.frame(SummarizedExperiment::colData(x), stringsAsFactors = FALSE)
    cells$cell <- colnames(x)
    cells <- .standardize_cells(cells, colnames(expr), sample_col, annotation_col, group_col, batch_col)
    if (is.null(embedding)) {
      embedding <- .extract_sce_embedding(x, reduced_dim)
    }
  } else if (inherits(x, "Seurat")) {
    expr <- .extract_seurat_matrix(x, assay_name = assay_name, expression_layer = expression_layer)
    cells <- as.data.frame(x[[]], stringsAsFactors = FALSE)
    cells$cell <- rownames(cells)
    cells <- .standardize_cells(cells, colnames(expr), sample_col, annotation_col, group_col, batch_col)
    if (is.null(embedding)) {
      embedding <- .extract_seurat_embedding(x, reduced_dim)
    }
  } else {
    stop("`x` must be a long-format data.frame, a list with `expr`, or a supported `SingleCellExperiment` / `Seurat` object.", call. = FALSE)
  }

  expr <- if (.is_sparse_matrix(expr)) expr else as.matrix(expr)
  if (!.is_sparse_matrix(expr)) {
    storage.mode(expr) <- "numeric"
  }
  if (is.null(embedding)) {
    embedding <- .compute_embedding_from_matrix(expr)
  } else {
    .assert_columns(embedding, c("cell", "UMAP1", "UMAP2"), "embedding")
    embedding <- as.data.frame(embedding, stringsAsFactors = FALSE)
    embedding <- embedding[match(colnames(expr), embedding$cell), c("cell", "UMAP1", "UMAP2"), drop = FALSE]
  }
  list(expr = expr, cells = cells, embedding = embedding)
}

.long_to_matrix <- function(expr, cell_order = NULL) {
  expr <- as.data.frame(expr, stringsAsFactors = FALSE)
  cell_levels <- cell_order %||% unique(expr$cell)
  gene_levels <- unique(expr$gene)
  i <- match(expr$gene, gene_levels)
  j <- match(expr$cell, cell_levels)
  if (requireNamespace("Matrix", quietly = TRUE)) {
    mat <- Matrix::sparseMatrix(
      i = i,
      j = j,
      x = as.numeric(expr$value),
      dims = c(length(gene_levels), length(cell_levels)),
      dimnames = list(gene_levels, cell_levels)
    )
    return(mat)
  }
  mat <- stats::xtabs(value ~ gene + cell, data = data.frame(gene = expr$gene, cell = expr$cell, value = expr$value, stringsAsFactors = FALSE))
  mat[, cell_levels, drop = FALSE]
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
  preferred <- c("counts", "logcounts")
  hit <- preferred[preferred %in% available]
  if (length(hit)) hit[[1]] else available[[1]]
}

.extract_sce_embedding <- function(x, reduced_dim) {
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    return(NULL)
  }
  emb <- tryCatch(SingleCellExperiment::reducedDim(x, reduced_dim), error = function(e) NULL)
  if (is.null(emb) || ncol(as.matrix(emb)) < 2L) {
    return(NULL)
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
    return(NULL)
  }
  emb <- tryCatch(SeuratObject::Embeddings(x[[reduced_dim]]), error = function(e) NULL)
  if (is.null(emb) || ncol(as.matrix(emb)) < 2L) {
    return(NULL)
  }
  emb <- as.matrix(emb)
  data.frame(cell = rownames(emb), UMAP1 = emb[, 1], UMAP2 = emb[, 2], stringsAsFactors = FALSE)
}
