.marker_test_table <- function(obj,
                               exclude_pattern = "^(RPL|RPS|MT-|MALAT1$|MTRNR|HSP)",
                               min_pct_in_type = 0.1,
                               max_features = 2000L) {
  stats_table <- .marker_statistics(obj, exclude_pattern = exclude_pattern)
  stats_table <- stats_table[stats_table$pct_in_type >= min_pct_in_type | stats_table$pct_rest >= min_pct_in_type / 2, , drop = FALSE]
  if (nrow(stats_table) > max_features * length(unique(stats_table$cell_type))) {
    feature_order <- unique(stats_table[order(-stats_table$marker_table_rank), "gene"])
    keep <- utils::head(feature_order, max_features)
    stats_table <- stats_table[stats_table$gene %in% keep, , drop = FALSE]
  }
  mat <- .analysis_matrix(obj)
  ann <- as.character(obj$cells[[obj$annotation_col]])
  out <- lapply(split(stats_table, stats_table$cell_type), function(df) {
    cell_type <- df$cell_type[[1]]
    idx <- ann == cell_type
    rest <- !idx
    p_values <- vapply(df$gene, function(gene) {
      values <- as.numeric(mat[gene, ])
      tryCatch(stats::wilcox.test(values[idx], values[rest], exact = FALSE)$p.value, error = function(e) NA_real_)
    }, numeric(1))
    df$p_value <- p_values
    df$p_adj <- stats::p.adjust(df$p_value, method = "BH")
    df
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out[order(out$cell_type, out$p_adj, -out$marker_score), , drop = FALSE]
}

#' Test markers with lightweight per-cell Wilcoxon statistics
#'
#' @param obj A `scdown_obj`.
#' @param top_n Number of marker tests to keep per cell type.
#' @param exclude_pattern Pattern for low-information genes to exclude.
#' @param min_pct_in_type Minimum detection rate inside the target cell type.
#' @param max_features Maximum features to test after prefiltering.
#' @param outdir Optional output directory.
#'
#' @return A list with marker test tables and plot paths.
#'
#' @details
#' This is a quick review-layer test. It is useful for prioritizing marker
#' candidates, but it is not a replacement for pseudobulk, mixed-model, or
#' covariate-aware specialist workflows.
#' @export
test_markers <- function(obj,
                         top_n = 10L,
                         exclude_pattern = "^(RPL|RPS|MT-|MALAT1$|MTRNR|HSP)",
                         min_pct_in_type = 0.1,
                         max_features = 2000L,
                         outdir = NULL) {
  stopifnot(inherits(obj, "scdown_obj"))
  marker_test_table <- .marker_test_table(
    obj = obj,
    exclude_pattern = exclude_pattern,
    min_pct_in_type = min_pct_in_type,
    max_features = max_features
  )
  top_marker_tests <- .top_n_per_group(marker_test_table[order(-marker_test_table$marker_score), , drop = FALSE], "cell_type", "marker_score", n = top_n)
  top_marker_tests <- top_marker_tests[order(top_marker_tests$cell_type, top_marker_tests$p_adj, -top_marker_tests$marker_score), , drop = FALSE]
  tables <- list(
    marker_test_table = marker_test_table,
    top_marker_tests = top_marker_tests
  )
  .write_tables(tables, outdir)
  plotters <- list(
    marker_test_heatmap = function() {
      if (!nrow(top_marker_tests)) {
        return(.plot_message("No marker tests"))
      }
      mat <- stats::xtabs(marker_score ~ cell_type + gene, data = top_marker_tests)
      graphics::image(seq_len(ncol(mat)), seq_len(nrow(mat)), t(mat[nrow(mat):1, , drop = FALSE]), axes = FALSE, col = grDevices::hcl.colors(20, "Plasma"), main = "Tested markers")
      graphics::axis(1, at = seq_len(ncol(mat)), labels = colnames(mat), las = 2, cex.axis = 0.6)
      graphics::axis(2, at = seq_len(nrow(mat)), labels = rev(rownames(mat)), las = 2, cex.axis = 0.7)
    }
  )
  list(tables = tables, plots = .save_plot_paths(plotters, outdir))
}

#' Test composition differences across groups
#'
#' @param obj A `scdown_obj`.
#' @param level Either `"cell_type"` or `"lineage"`.
#' @param outdir Optional output directory.
#'
#' @return A list with composition test tables and plot paths.
#'
#' @details
#' This is a lightweight sample-level fraction test. It is useful for quick
#' screening, but specialist composition or abundance models should take over
#' when composition is a primary result.
#' @export
test_composition <- function(obj, level = c("cell_type", "lineage"), outdir = NULL) {
  stopifnot(inherits(obj, "scdown_obj"))
  level <- match.arg(level)
  if (!obj$lineage_col %in% names(obj$cells)) {
    obj <- collapse_lineage(obj)
  }
  tables_comp <- .composition_tables(obj)
  df <- if (identical(level, "cell_type")) tables_comp$composition_by_sample else tables_comp$composition_major_lineage
  if (!"group" %in% names(df)) {
    stop("`test_composition()` requires `group_col` with at least two groups.", call. = FALSE)
  }
  feature_col <- if (identical(level, "cell_type")) "cell_type" else "lineage"
  test_table <- do.call(rbind, lapply(split(df, df[[feature_col]]), function(x) {
    test_res <- .group_test(x$fraction, x$group)
    data.frame(
      feature = x[[feature_col]][[1]],
      method = test_res$method,
      p_value = test_res$p_value,
      effect = test_res$effect,
      note = test_res$note,
      stringsAsFactors = FALSE
    )
  }))
  test_table$p_adj <- stats::p.adjust(test_table$p_value, method = "BH")
  tables <- list(composition_test_table = test_table)
  .write_tables(tables, outdir)
  plotters <- list(
    composition_test_overview = function() {
      df_plot <- test_table
      df_plot$score <- -log10(df_plot$p_adj)
      df_plot$score[!is.finite(df_plot$score)] <- 0
      graphics::barplot(df_plot$score, names.arg = df_plot$feature, las = 2, col = "slateblue", main = paste("Composition test:", level), ylab = "-log10 adjusted p")
    }
  )
  list(tables = tables, plots = .save_plot_paths(plotters, outdir))
}

#' Test signature differences across groups
#'
#' @param obj A `scdown_obj`.
#' @param signatures Character vector of signature names or a named list of signatures.
#' @param signature_sets Optional named list of custom signatures.
#' @param outdir Optional output directory.
#'
#' @return A list with signature test tables and plot paths.
#'
#' @details
#' This is a compact group comparison over sample-level mean signature scores.
#' It is intended for review and prioritization, not as a replacement for more
#' flexible specialist modeling.
#' @export
test_signature <- function(obj, signatures, signature_sets = NULL, outdir = NULL) {
  stopifnot(inherits(obj, "scdown_obj"))
  if (!obj$group_col %in% names(obj$cells)) {
    stop("`test_signature()` requires `group_col` with at least two groups.", call. = FALSE)
  }
  sig_plot <- plot_signature(obj, signatures = signatures, signature_sets = signature_sets)
  score_long <- sig_plot$tables$signature_score_by_cell
  sample_info <- unique(obj$cells[, c("cell", obj$sample_col, obj$group_col, obj$annotation_col), drop = FALSE])
  names(sample_info) <- c("cell", "sample", "group", "cell_type")
  score_long <- merge(score_long, sample_info, by = c("cell", "cell_type"), all.x = TRUE)
  sample_scores <- stats::aggregate(score ~ signature + sample + group + cell_type, data = score_long, FUN = mean)
  test_table <- do.call(rbind, lapply(split(sample_scores, interaction(sample_scores$signature, sample_scores$cell_type, drop = TRUE)), function(x) {
    test_res <- .group_test(x$score, x$group)
    data.frame(
      signature = x$signature[[1]],
      cell_type = x$cell_type[[1]],
      method = test_res$method,
      p_value = test_res$p_value,
      effect = test_res$effect,
      note = test_res$note,
      stringsAsFactors = FALSE
    )
  }))
  test_table$p_adj <- stats::p.adjust(test_table$p_value, method = "BH")
  tables <- list(
    signature_sample_scores = sample_scores,
    signature_test_table = test_table
  )
  .write_tables(tables, outdir)
  plotters <- list(
    signature_test_overview = function() {
      df_plot <- test_table
      if (!nrow(df_plot)) {
        return(.plot_message("No signature tests"))
      }
      df_plot$score <- -log10(df_plot$p_adj)
      df_plot$score[!is.finite(df_plot$score)] <- 0
      graphics::barplot(df_plot$score, names.arg = paste(df_plot$signature, df_plot$cell_type, sep = ":"), las = 2, col = "darkgreen", main = "Signature tests", ylab = "-log10 adjusted p")
    }
  )
  list(tables = tables, plots = .save_plot_paths(plotters, outdir))
}

#' Test communication scores by permutation
#'
#' @param obj A `scdown_obj`.
#' @param lr_pairs Optional data frame with custom ligand-receptor pairs.
#' @param n_perm Number of label permutations.
#' @param seed Random seed.
#' @param min_detected Minimum detection rate.
#' @param min_specificity Minimum specificity score.
#' @param exclude_pattern Pattern for labels to exclude.
#' @param outdir Optional output directory.
#'
#' @return A list with communication test tables and plot paths.
#'
#' @details
#' This is an exploratory permutation layer over the built-in communication
#' score. Use a specialist framework when communication is a central claim.
#' @export
test_communication <- function(obj,
                               lr_pairs = NULL,
                               n_perm = 50L,
                               seed = 1L,
                               min_detected = 0.1,
                               min_specificity = 0.15,
                               exclude_pattern = "(?i)unassigned|unknown|doublet",
                               outdir = NULL) {
  stopifnot(inherits(obj, "scdown_obj"))
  observed <- .communication_tables(
    obj = obj,
    lr_pairs = lr_pairs,
    min_detected = min_detected,
    min_specificity = min_specificity,
    exclude_pattern = exclude_pattern
  )
  observed_pairs <- observed$communication_pair_summary
  if (!nrow(observed_pairs)) {
    tables <- list(
      communication_test_table = data.frame(sender = character(), receiver = character(), score = numeric(), p_empirical = numeric(), stringsAsFactors = FALSE)
    )
    .write_tables(tables, outdir)
    return(list(tables = tables, plots = .save_plot_paths(list(communication_test_overview = function() .plot_message("No communication pairs")), outdir)))
  }
  set.seed(seed)
  labels <- obj$cells[[obj$annotation_col]]
  perm_scores <- replicate(n_perm, {
    perm <- sample(labels)
    perm_table <- .communication_tables(
      obj = obj,
      labels = perm,
      lr_pairs = lr_pairs,
      min_detected = min_detected,
      min_specificity = min_specificity,
      exclude_pattern = exclude_pattern
    )$communication_pair_summary
    merged <- merge(observed_pairs[, c("sender", "receiver")], perm_table, by = c("sender", "receiver"), all.x = TRUE)
    merged$score[is.na(merged$score)] <- 0
    merged$score
  })
  if (is.null(dim(perm_scores))) {
    perm_scores <- matrix(perm_scores, ncol = 1L)
  }
  communication_test_table <- observed_pairs
  communication_test_table$p_empirical <- (1 + rowSums(perm_scores >= communication_test_table$score, na.rm = TRUE)) / (ncol(perm_scores) + 1)
  communication_test_table$p_adj <- stats::p.adjust(communication_test_table$p_empirical, method = "BH")
  tables <- list(communication_test_table = communication_test_table)
  .write_tables(tables, outdir)
  plotters <- list(
    communication_test_overview = function() {
      df <- communication_test_table
      df$score_plot <- -log10(df$p_adj)
      df$score_plot[!is.finite(df$score_plot)] <- 0
      graphics::barplot(df$score_plot, names.arg = paste(df$sender, df$receiver, sep = "->"), las = 2, col = "orange", main = "Communication permutation test", ylab = "-log10 adjusted p")
    }
  )
  list(tables = tables, plots = .save_plot_paths(plotters, outdir))
}
