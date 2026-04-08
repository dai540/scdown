#' Create a `scdown` object
#'
#' `scdown()` standardizes annotation-after single-cell inputs into a compact
#' object for cluster interpretation within one dataset. Supported inputs are:
#'
#' - a long table with `cell`, `gene`, and `value`
#' - a list with `expr`, `cells`, and optional `embedding`
#' - a `SingleCellExperiment`
#' - a `Seurat` object
#'
#' @param x Input data.
#' @param annotation_col Column containing cluster or annotation labels.
#' @param sample_col Optional column containing sample identifiers.
#' @param reduced_dim Reduced-dimension name used when extracting embeddings
#'   from Bioconductor or Seurat objects.
#' @param embedding Optional data frame with `cell`, `UMAP1`, and `UMAP2`.
#' @param assay_name Optional assay name for `SingleCellExperiment` or `Seurat`.
#' @param expression_layer Optional Seurat layer or slot name.
#' @param signatures Optional named list of custom signatures.
#' @param marker_panels Optional named list of custom marker panels.
#' @param lr_pairs Optional data frame with `ligand`, `receptor`, and `pathway`.
#'
#' @return An object of class `scdown_obj`.
#' @export
scdown <- function(x,
                   annotation_col = "annotation",
                   sample_col = "sample",
                   reduced_dim = "UMAP",
                   embedding = NULL,
                   assay_name = NULL,
                   expression_layer = NULL,
                   signatures = NULL,
                   marker_panels = NULL,
                   lr_pairs = NULL) {
  parsed <- .coerce_scdown_input(
    x = x,
    annotation_col = annotation_col,
    sample_col = sample_col,
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
    reduced_dim = reduced_dim,
    value_type = if (.is_count_like_matrix(parsed$expr)) "counts" else "continuous",
    resources = list(
      signatures = .validate_named_gene_sets(signatures, "signatures"),
      marker_panels = .validate_named_gene_sets(marker_panels, "marker_panels"),
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
  cat("  annotations:", length(unique(x$cells[[x$annotation_col]])), "\n")
  cat("  samples:", length(unique(x$cells[[x$sample_col]])), "\n")
  cat("  value type:", x$value_type, "\n")
  invisible(x)
}

.coerce_scdown_input <- function(x,
                                 annotation_col,
                                 sample_col,
                                 reduced_dim,
                                 embedding,
                                 assay_name,
                                 expression_layer) {
  if (is.data.frame(x)) {
    .assert_columns(x, c("cell", "gene", "value", annotation_col), "x")
    if (!sample_col %in% names(x)) {
      x[[sample_col]] <- "sample_1"
    }
    cells <- unique(x[, c("cell", sample_col, annotation_col), drop = FALSE])
    cells <- .standardize_cells(cells, unique(x$cell), sample_col, annotation_col)
    expr <- .long_to_matrix(x, cell_order = cells$cell)
  } else if (is.list(x) && !is.null(x$expr)) {
    expr <- .coerce_list_expr(
      x = x,
      annotation_col = annotation_col,
      sample_col = sample_col
    )
    cells <- attr(expr, "cells")
    if (is.null(embedding) && !is.null(x$embedding)) {
      embedding <- x$embedding
    }
  } else if (inherits(x, "SingleCellExperiment")) {
    assay_name <- .resolve_sce_assay(x, assay_name)
    expr <- SummarizedExperiment::assay(x, assay_name)
    cells <- as.data.frame(SummarizedExperiment::colData(x), stringsAsFactors = FALSE)
    cells$cell <- colnames(x)
    cells <- .standardize_cells(cells, colnames(expr), sample_col, annotation_col)
    if (is.null(embedding)) {
      embedding <- .extract_sce_embedding(x, reduced_dim)
    }
  } else if (inherits(x, "Seurat")) {
    expr <- .extract_seurat_matrix(x, assay_name = assay_name, expression_layer = expression_layer)
    cells <- as.data.frame(x[[]], stringsAsFactors = FALSE)
    cells$cell <- rownames(cells)
    cells <- .standardize_cells(cells, colnames(expr), sample_col, annotation_col)
    if (is.null(embedding)) {
      embedding <- .extract_seurat_embedding(x, reduced_dim)
    }
  } else {
    stop(
      "`x` must be a supported long table, list, `SingleCellExperiment`, or `Seurat` object.",
      call. = FALSE
    )
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

.coerce_list_expr <- function(x, annotation_col, sample_col) {
  if (is.data.frame(x$expr)) {
    expr_df <- as.data.frame(x$expr, stringsAsFactors = FALSE)
    .assert_columns(expr_df, c("cell", "gene", "value", annotation_col), "x$expr")
    cells <- if (!is.null(x$cells)) {
      x$cells
    } else {
      if (!sample_col %in% names(expr_df)) {
        expr_df[[sample_col]] <- "sample_1"
      }
      unique(expr_df[, c("cell", sample_col, annotation_col), drop = FALSE])
    }
    cells <- .standardize_cells(cells, unique(expr_df$cell), sample_col, annotation_col)
    expr <- .long_to_matrix(expr_df, cell_order = cells$cell)
    attr(expr, "cells") <- cells
    return(expr)
  }

  if (!(is.matrix(x$expr) || .is_sparse_matrix(x$expr))) {
    stop("List inputs must contain a matrix-like or long-table `expr` field.", call. = FALSE)
  }
  if (is.null(x$cells)) {
    stop("List inputs with matrix-like `expr` must supply `cells` metadata.", call. = FALSE)
  }
  expr <- x$expr
  if (is.null(rownames(expr))) {
    rownames(expr) <- paste0("gene_", seq_len(nrow(expr)))
  }
  if (is.null(colnames(expr))) {
    stop("Matrix-like `expr` inputs must have column names.", call. = FALSE)
  }
  cells <- .standardize_cells(x$cells, colnames(expr), sample_col, annotation_col)
  attr(expr, "cells") <- cells
  expr
}

.long_to_matrix <- function(expr, cell_order) {
  expr <- as.data.frame(expr, stringsAsFactors = FALSE)
  gene_levels <- unique(as.character(expr$gene))
  cell_levels <- unique(as.character(cell_order))
  i <- match(expr$gene, gene_levels)
  j <- match(expr$cell, cell_levels)
  if (requireNamespace("Matrix", quietly = TRUE)) {
    return(Matrix::sparseMatrix(
      i = i,
      j = j,
      x = as.numeric(expr$value),
      dims = c(length(gene_levels), length(cell_levels)),
      dimnames = list(gene_levels, cell_levels)
    ))
  }
  mat <- stats::xtabs(value ~ gene + cell, data = data.frame(
    gene = expr$gene,
    cell = expr$cell,
    value = expr$value,
    stringsAsFactors = FALSE
  ))
  mat[, cell_levels, drop = FALSE]
}

.resolve_sce_assay <- function(x, assay_name = NULL) {
  available <- names(SummarizedExperiment::assays(x))
  if (!length(available)) {
    stop("No assays were found in the `SingleCellExperiment` input.", call. = FALSE)
  }
  if (!is.null(assay_name)) {
    if (!assay_name %in% available) {
      stop("Requested assay was not found in the `SingleCellExperiment` input.", call. = FALSE)
    }
    return(assay_name)
  }
  preferred <- c("counts", "logcounts")
  hit <- preferred[preferred %in% available]
  if (length(hit)) hit[[1]] else available[[1]]
}

.extract_sce_embedding <- function(x, reduced_dim) {
  emb <- tryCatch(SingleCellExperiment::reducedDim(x, reduced_dim), error = function(e) NULL)
  if (is.null(emb) || ncol(as.matrix(emb)) < 2L) {
    return(NULL)
  }
  emb <- as.matrix(emb)
  data.frame(
    cell = rownames(emb) %||% colnames(x),
    UMAP1 = emb[, 1],
    UMAP2 = emb[, 2],
    stringsAsFactors = FALSE
  )
}

.extract_seurat_matrix <- function(x, assay_name = NULL, expression_layer = NULL) {
  assay_name <- assay_name %||% SeuratObject::DefaultAssay(x)
  preferred <- unique(c(expression_layer %||% "data", "counts", "data", "scale.data"))
  if (exists("LayerData", where = asNamespace("SeuratObject"), mode = "function")) {
    for (layer_name in preferred) {
      mat <- tryCatch(
        suppressWarnings(SeuratObject::LayerData(x, assay = assay_name, layer = layer_name)),
        error = function(e) NULL
      )
      if (!is.null(mat) && nrow(mat) > 0L && ncol(mat) > 0L) {
        return(mat)
      }
    }
  }
  for (slot_name in preferred) {
    mat <- tryCatch(
      suppressWarnings(SeuratObject::GetAssayData(x, assay = assay_name, slot = slot_name)),
      error = function(e) NULL
    )
    if (!is.null(mat) && nrow(mat) > 0L && ncol(mat) > 0L) {
      return(mat)
    }
  }
  stop("No usable expression layer was found in the Seurat input.", call. = FALSE)
}

.extract_seurat_embedding <- function(x, reduced_dim) {
  emb <- tryCatch(SeuratObject::Embeddings(x[[reduced_dim]]), error = function(e) NULL)
  if (is.null(emb) || ncol(as.matrix(emb)) < 2L) {
    return(NULL)
  }
  emb <- as.matrix(emb)
  data.frame(
    cell = rownames(emb),
    UMAP1 = emb[, 1],
    UMAP2 = emb[, 2],
    stringsAsFactors = FALSE
  )
}
