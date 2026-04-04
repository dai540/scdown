.assert_columns <- function(x, required, what = "data") {
  miss <- setdiff(required, names(x))
  if (length(miss) > 0L) {
    stop(sprintf("`%s` is missing required columns: %s.", what, paste(miss, collapse = ", ")), call. = FALSE)
  }
}

.existing_columns <- function(x, cols) {
  cols <- cols[!is.na(cols) & nzchar(cols)]
  cols[cols %in% names(x)]
}

.factor_info <- function(x) {
  vals <- unique(as.character(x))
  fac <- factor(as.character(x), levels = vals)
  list(index = as.integer(fac), levels = vals)
}

.safe_mean <- function(x) {
  if (!length(x)) {
    return(NA_real_)
  }
  mean(x, na.rm = TRUE)
}

.dir_create <- function(path) {
  if (!is.null(path) && !dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}

.write_tables <- function(tables, outdir) {
  if (is.null(outdir)) {
    return(invisible(NULL))
  }
  .dir_create(outdir)
  for (nm in names(tables)) {
    utils::write.csv(tables[[nm]], file.path(outdir, paste0(nm, ".csv")), row.names = FALSE)
  }
  invisible(NULL)
}

.with_png <- function(path, expr, width = 1200, height = 900) {
  grDevices::png(path, width = width, height = height, res = 120)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
  invisible(path)
}

.save_plot_paths <- function(plotters, outdir) {
  if (is.null(outdir)) {
    return(character())
  }
  .dir_create(outdir)
  out <- character()
  for (nm in names(plotters)) {
    path <- file.path(outdir, paste0(nm, ".png"))
    .with_png(path, plotters[[nm]]())
    out[nm] <- path
  }
  out
}

.plot_message <- function(label) {
  graphics::plot.new()
  graphics::text(0.5, 0.5, label)
}

.top_n_per_group <- function(df, group_col, score_col, n = 10L) {
  split_df <- split(df, df[[group_col]])
  out <- lapply(split_df, function(x) {
    x <- x[order(x[[score_col]], decreasing = TRUE), , drop = FALSE]
    utils::head(x, n)
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

.is_sparse_matrix <- function(x) {
  inherits(x, "sparseMatrix")
}

.as_numeric_probe <- function(x, max_n = 10000L) {
  if (.is_sparse_matrix(x)) {
    probe <- x@x
  } else {
    probe <- as.numeric(x)
  }
  probe <- probe[is.finite(probe)]
  if (length(probe) > max_n) {
    probe <- probe[seq_len(max_n)]
  }
  probe
}

.is_count_like_matrix <- function(x) {
  probe <- .as_numeric_probe(x)
  if (!length(probe)) {
    return(TRUE)
  }
  if (any(probe < 0)) {
    return(FALSE)
  }
  all(abs(probe - round(probe)) < 1e-8)
}

.matrix_col_sums <- function(x) {
  if (.is_sparse_matrix(x)) {
    Matrix::colSums(x)
  } else {
    colSums(x)
  }
}

.matrix_row_means <- function(x) {
  if (.is_sparse_matrix(x)) {
    as.numeric(Matrix::rowMeans(x))
  } else {
    rowMeans(x)
  }
}

.matrix_row_vars <- function(x) {
  means <- .matrix_row_means(x)
  sq_means <- .matrix_row_means(x ^ 2)
  pmax(sq_means - means^2, 0)
}

.normalize_matrix_counts <- function(mat) {
  totals <- .matrix_col_sums(mat)
  totals[totals <= 0] <- 1
  if (.is_sparse_matrix(mat)) {
    norm <- Matrix::t(Matrix::t(mat) / totals)
    norm <- norm * 1e4
    norm@x <- log1p(norm@x)
    return(norm)
  }
  log1p(t(t(as.matrix(mat)) / totals) * 1e4)
}

.analysis_matrix <- function(obj) {
  if (identical(obj$value_type, "counts")) {
    return(.normalize_matrix_counts(obj$expr))
  }
  obj$expr
}

.detected_matrix <- function(obj) {
  obj$expr > 0
}

.deterministic_embedding <- function(cells) {
  n <- length(cells)
  if (!n) {
    return(data.frame(cell = character(), UMAP1 = numeric(), UMAP2 = numeric(), stringsAsFactors = FALSE))
  }
  theta <- seq(0, 2 * pi, length.out = n + 1L)[seq_len(n)]
  radius <- 1 + ((seq_len(n) - 1L) %% 7L) / 10
  data.frame(
    cell = cells,
    UMAP1 = cos(theta) * radius,
    UMAP2 = sin(theta) * radius,
    stringsAsFactors = FALSE
  )
}

.compute_embedding_from_matrix <- function(mat, n_features = 800L) {
  if (!ncol(mat)) {
    return(data.frame(cell = character(), UMAP1 = numeric(), UMAP2 = numeric(), stringsAsFactors = FALSE))
  }
  if (nrow(mat) < 2L || ncol(mat) < 3L) {
    return(.deterministic_embedding(colnames(mat)))
  }
  norm <- .normalize_matrix_counts(mat)
  vars <- .matrix_row_vars(norm)
  vars[is.na(vars)] <- 0
  keep <- order(vars, decreasing = TRUE)[seq_len(min(length(vars), n_features))]
  emb_input <- t(as.matrix(norm[keep, , drop = FALSE]))
  if (requireNamespace("uwot", quietly = TRUE) && nrow(emb_input) >= 10L) {
    set.seed(1)
    emb <- uwot::umap(
      emb_input,
      n_components = 2,
      n_neighbors = min(15, max(2, nrow(emb_input) - 1L)),
      min_dist = 0.3,
      metric = "cosine",
      verbose = FALSE
    )
  } else {
    emb <- stats::prcomp(emb_input, center = TRUE, scale. = TRUE)$x[, 1:2, drop = FALSE]
  }
  data.frame(cell = colnames(mat), UMAP1 = emb[, 1], UMAP2 = emb[, 2], stringsAsFactors = FALSE)
}

.rank_signature_score <- function(mat, genes) {
  hit <- intersect(genes, rownames(mat))
  if (!length(hit)) {
    return(stats::setNames(rep(0, ncol(mat)), colnames(mat)))
  }
  ranks <- apply(as.matrix(mat), 2, function(x) rank(-x, ties.method = "average"))
  if (is.null(dim(ranks))) {
    ranks <- matrix(ranks, ncol = 1L, dimnames = list(rownames(mat), colnames(mat)))
  }
  denom <- max(nrow(ranks) - 1L, 1L)
  sig_ranks <- ranks[hit, , drop = FALSE]
  score <- colMeans(1 - ((sig_ranks - 1) / denom))
  stats::setNames(as.numeric(score), colnames(mat))
}

.group_test <- function(values, groups) {
  keep <- is.finite(values) & !is.na(groups)
  values <- values[keep]
  groups <- as.character(groups[keep])
  lev <- unique(groups)
  if (length(lev) < 2L) {
    return(list(method = "insufficient_groups", p_value = NA_real_, effect = NA_real_, note = "need at least two groups"))
  }
  sample_counts <- table(groups)
  if (any(sample_counts < 2L)) {
    return(list(method = "insufficient_replicates", p_value = NA_real_, effect = NA_real_, note = "need at least two samples per group"))
  }
  group_means <- tapply(values, groups, mean, na.rm = TRUE)
  effect <- max(group_means) - min(group_means)
  if (!is.finite(effect) || effect == 0 || stats::var(values) == 0) {
    method <- if (length(lev) == 2L) "wilcox" else "kruskal"
    return(list(method = method, p_value = 1, effect = 0, note = "no variation between groups"))
  }
  if (length(lev) == 2L) {
    p_value <- tryCatch(stats::wilcox.test(values ~ factor(groups), exact = FALSE)$p.value, error = function(e) NA_real_)
    return(list(method = "wilcox", p_value = p_value, effect = effect, note = ""))
  }
  p_value <- tryCatch(stats::kruskal.test(values ~ factor(groups))$p.value, error = function(e) NA_real_)
  list(method = "kruskal", p_value = p_value, effect = effect, note = "")
}

.standardize_cells <- function(cells, cell_names, sample_col, annotation_col, group_col, batch_col) {
  cells <- as.data.frame(cells, stringsAsFactors = FALSE)
  if (!"cell" %in% names(cells)) {
    cells$cell <- rownames(cells) %||% cell_names
  }
  cells$cell <- as.character(cells$cell)
  .assert_columns(cells, c("cell", sample_col, annotation_col), "cells")
  keep <- .existing_columns(cells, c("cell", sample_col, annotation_col, group_col, batch_col))
  cells <- cells[, unique(keep), drop = FALSE]
  miss <- setdiff(cell_names, cells$cell)
  if (length(miss)) {
    stop("Cell metadata is missing rows for some expression columns.", call. = FALSE)
  }
  cells <- cells[match(cell_names, cells$cell), , drop = FALSE]
  rownames(cells) <- NULL
  cells
}

.matrix_to_long_subset <- function(raw_mat, analysis_mat, genes, cells, sample_col, annotation_col) {
  hit <- intersect(genes, rownames(raw_mat))
  if (!length(hit)) {
    return(data.frame(
      cell = character(),
      sample = character(),
      cell_type = character(),
      gene = character(),
      value = numeric(),
      analysis_value = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  out <- lapply(hit, function(gene) {
    data.frame(
      cell = colnames(raw_mat),
      sample = cells[[sample_col]],
      cell_type = cells[[annotation_col]],
      gene = gene,
      value = as.numeric(raw_mat[gene, ]),
      analysis_value = as.numeric(analysis_mat[gene, ]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
