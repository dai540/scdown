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
  do.call(rbind, out)
}

.is_count_like <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x) || any(x < 0)) {
    return(FALSE)
  }
  probe <- x[seq_len(min(length(x), 10000L))]
  all(abs(probe - round(probe)) < 1e-8)
}

.normalize_count_values <- function(values, cells) {
  totals <- stats::aggregate(values, by = list(cell = cells), FUN = sum)
  names(totals)[2] <- "library_size"
  lib_size <- totals$library_size[match(cells, totals$cell)]
  lib_size[is.na(lib_size) | lib_size <= 0] <- 1
  log1p((values / lib_size) * 1e4)
}

.analysis_expr <- function(obj) {
  expr <- obj$expr
  out <- expr
  count_like <- .is_count_like(expr$value)
  if (count_like) {
    out$analysis_value <- .normalize_count_values(expr$value, expr$cell)
    out$detected <- expr$value > 0
  } else {
    out$analysis_value <- expr$value
    out$detected <- expr$value > 0
  }
  out$value_type <- if (count_like) "counts" else "continuous"
  out
}

.expr_matrix <- function(expr, value_col = "analysis_value") {
  stats::xtabs(
    value ~ gene + cell,
    data = data.frame(
      gene = expr$gene,
      cell = expr$cell,
      value = expr[[value_col]],
      stringsAsFactors = FALSE
    )
  )
}

.cell_type_sizes <- function(expr, annotation_col) {
  cell_meta <- unique(expr[, c("cell", annotation_col), drop = FALSE])
  stats::aggregate(
    cell ~ cell_type,
    data = setNames(cell_meta, c("cell", "cell_type")),
    FUN = length
  )
}

.normalize_matrix_counts <- function(mat) {
  mat <- as.matrix(mat)
  totals <- colSums(mat)
  totals[totals <= 0] <- 1
  log1p(t(t(mat) / totals) * 1e4)
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
  vars <- apply(norm, 1, stats::var)
  vars[is.na(vars)] <- 0
  keep <- order(vars, decreasing = TRUE)[seq_len(min(length(vars), n_features))]
  emb_input <- t(norm[keep, , drop = FALSE])
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
  ranks <- apply(mat, 2, function(x) rank(-x, ties.method = "average"))
  if (is.null(dim(ranks))) {
    ranks <- matrix(ranks, ncol = 1L, dimnames = list(rownames(mat), colnames(mat)))
  }
  denom <- max(nrow(ranks) - 1L, 1L)
  sig_ranks <- ranks[hit, , drop = FALSE]
  score <- colMeans(1 - ((sig_ranks - 1) / denom))
  stats::setNames(as.numeric(score), colnames(mat))
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
