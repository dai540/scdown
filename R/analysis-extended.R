.score_signature <- function(obj, genes, score_name) {
  expr_all <- .analysis_expr(obj)
  mat <- .expr_matrix(expr_all, value_col = "analysis_value")
  score <- data.frame(
    cell = colnames(mat),
    score = .rank_signature_score(mat, genes),
    stringsAsFactors = FALSE
  )
  names(score)[2] <- score_name
  score
}

#' Check annotation with marker panels
#'
#' @param obj A `scdown_obj`.
#' @param outdir Optional output directory.
#'
#' @return A list with annotation-check tables and plots.
#' @export
check_annotation <- function(obj, outdir = NULL) {
  stopifnot(inherits(obj, "scdown_obj"))
  panels <- immune_marker_panels()
  panel_names <- names(panels)
  cell_scores <- unique(obj$cells["cell"])
  for (nm in panel_names) {
    cell_scores <- merge(cell_scores, .score_signature(obj, panels[[nm]], nm), by = "cell", all.x = TRUE)
  }
  cell_scores[is.na(cell_scores)] <- 0
  long <- reshape(
    merge(cell_scores, obj$cells[, c("cell", obj$annotation_col), drop = FALSE], by = "cell", all.x = TRUE),
    varying = panel_names,
    v.names = "score",
    timevar = "panel",
    times = panel_names,
    direction = "long"
  )
  rownames(long) <- NULL
  annotation_check_score <- stats::aggregate(score ~ panel + long[[obj$annotation_col]], data = long, FUN = mean)
  names(annotation_check_score)[2] <- "cell_type"
  annotation_marker_panel <- data.frame(
    panel = rep(names(panels), lengths(panels)),
    gene = unlist(panels, use.names = FALSE),
    stringsAsFactors = FALSE
  )
  tables <- list(
    annotation_marker_panel = annotation_marker_panel,
    annotation_check_score = annotation_check_score
  )
  .write_tables(tables, outdir)
  plotters <- list(
    annotation_check_heatmap = function() {
      mat <- stats::xtabs(score ~ cell_type + panel, data = annotation_check_score)
      graphics::image(seq_len(ncol(mat)), seq_len(nrow(mat)), t(mat[nrow(mat):1, , drop = FALSE]), axes = FALSE, col = grDevices::hcl.colors(20, "Temps"), main = "Annotation check")
      graphics::axis(1, at = seq_len(ncol(mat)), labels = colnames(mat), las = 2, cex.axis = 0.7)
      graphics::axis(2, at = seq_len(nrow(mat)), labels = rev(rownames(mat)), las = 2, cex.axis = 0.7)
    }
  )
  list(tables = tables, plots = .save_plot_paths(plotters, outdir))
}

#' Plot signature scores
#'
#' @param obj A `scdown_obj`.
#' @param signatures Character vector of built-in signature names.
#' @param outdir Optional output directory.
#'
#' @return A list with signature tables and plots.
#' @export
plot_signature <- function(obj, signatures, outdir = NULL) {
  stopifnot(inherits(obj, "scdown_obj"))
  signatures <- unique(as.character(signatures))
  cell_scores <- unique(obj$cells["cell"])
  for (sig in signatures) {
    cell_scores <- merge(cell_scores, .score_signature(obj, get_signature(sig), sig), by = "cell", all.x = TRUE)
  }
  cell_scores[is.na(cell_scores)] <- 0
  score_long <- reshape(
    merge(cell_scores, obj$cells[, c("cell", obj$annotation_col), drop = FALSE], by = "cell", all.x = TRUE),
    varying = signatures,
    v.names = "score",
    timevar = "signature",
    times = signatures,
    direction = "long"
  )
  rownames(score_long) <- NULL
  score_by_celltype <- stats::aggregate(score ~ signature + score_long[[obj$annotation_col]], data = score_long, FUN = mean)
  names(score_by_celltype)[2] <- "cell_type"
  tables <- list(
    signature_score_by_cell = score_long,
    signature_score_by_celltype = score_by_celltype
  )
  .write_tables(tables, outdir)
  plotters <- list(
    signature_dotplot = function() {
      df <- score_by_celltype
      cell_info <- .factor_info(df$cell_type)
      sig_info <- .factor_info(df$signature)
      graphics::plot(cell_info$index, sig_info$index, cex = pmax(df$score, 0.3), pch = 19, col = "purple", xlab = "Cell type", ylab = "Signature", xaxt = "n", yaxt = "n", main = "Signatures")
      graphics::axis(1, at = seq_along(cell_info$levels), labels = cell_info$levels, las = 2, cex.axis = 0.7)
      graphics::axis(2, at = seq_along(sig_info$levels), labels = sig_info$levels, las = 2, cex.axis = 0.8)
    },
    signature_heatmap = function() {
      mat <- stats::xtabs(score ~ cell_type + signature, data = score_by_celltype)
      graphics::image(seq_len(ncol(mat)), seq_len(nrow(mat)), t(mat[nrow(mat):1, , drop = FALSE]), axes = FALSE, col = grDevices::hcl.colors(20, "BluYl"), main = "Signature heatmap")
      graphics::axis(1, at = seq_len(ncol(mat)), labels = colnames(mat), las = 2)
      graphics::axis(2, at = seq_len(nrow(mat)), labels = rev(rownames(mat)), las = 2)
    }
  )
  list(tables = tables, plots = .save_plot_paths(plotters, outdir))
}

#' Plot expression heatmaps
#'
#' @param obj A `scdown_obj`.
#' @param signatures Signatures to include in the signature heatmap.
#' @param top_n_markers Number of markers per cell type to show.
#' @param outdir Optional output directory.
#'
#' @return A list with heatmap tables and plots.
#' @export
plot_heatmap <- function(obj,
                         signatures = c("cytotoxic", "exhaustion", "antigen_presentation"),
                         top_n_markers = 5L,
                         outdir = NULL) {
  stopifnot(inherits(obj, "scdown_obj"))
  expr <- .analysis_expr(obj)
  average_expression <- stats::aggregate(analysis_value ~ expr[[obj$annotation_col]] + gene, data = expr, FUN = mean)
  names(average_expression)[1:3] <- c("cell_type", "gene", "mean_expression")
  marker_res <- find_markers(obj, top_n = top_n_markers)
  marker_heatmap_matrix <- marker_res$tables$top_marker_per_celltype
  signature_res <- plot_signature(obj, signatures = signatures)
  signature_heatmap_matrix <- signature_res$tables$signature_score_by_celltype
  tables <- list(
    average_expression = average_expression,
    marker_heatmap_matrix = marker_heatmap_matrix,
    signature_heatmap_matrix = signature_heatmap_matrix
  )
  .write_tables(tables, outdir)
  plotters <- list(
    average_expression_heatmap = function() {
      top_genes <- unique(marker_heatmap_matrix$gene)
      df <- average_expression[average_expression$gene %in% top_genes, , drop = FALSE]
      if (!nrow(df)) {
        return(.plot_message("No genes available"))
      }
      mat <- stats::xtabs(mean_expression ~ cell_type + gene, data = df)
      graphics::image(seq_len(ncol(mat)), seq_len(nrow(mat)), t(mat[nrow(mat):1, , drop = FALSE]), axes = FALSE, col = grDevices::hcl.colors(20, "Teal"), main = "Average expression")
      graphics::axis(1, at = seq_len(ncol(mat)), labels = colnames(mat), las = 2, cex.axis = 0.6)
      graphics::axis(2, at = seq_len(nrow(mat)), labels = rev(rownames(mat)), las = 2, cex.axis = 0.7)
    },
    marker_heatmap = function() {
      mat <- stats::xtabs(marker_score ~ cell_type + gene, data = marker_heatmap_matrix)
      graphics::image(seq_len(ncol(mat)), seq_len(nrow(mat)), t(mat[nrow(mat):1, , drop = FALSE]), axes = FALSE, col = grDevices::hcl.colors(20, "Inferno"), main = "Marker heatmap")
      graphics::axis(1, at = seq_len(ncol(mat)), labels = colnames(mat), las = 2, cex.axis = 0.6)
      graphics::axis(2, at = seq_len(nrow(mat)), labels = rev(rownames(mat)), las = 2, cex.axis = 0.7)
    },
    signature_heatmap = function() {
      mat <- stats::xtabs(score ~ cell_type + signature, data = signature_heatmap_matrix)
      graphics::image(seq_len(ncol(mat)), seq_len(nrow(mat)), t(mat[nrow(mat):1, , drop = FALSE]), axes = FALSE, col = grDevices::hcl.colors(20, "Viridis"), main = "Signature heatmap")
      graphics::axis(1, at = seq_len(ncol(mat)), labels = colnames(mat), las = 2, cex.axis = 0.7)
      graphics::axis(2, at = seq_len(nrow(mat)), labels = rev(rownames(mat)), las = 2, cex.axis = 0.7)
    }
  )
  list(tables = tables, plots = .save_plot_paths(plotters, outdir))
}

#' Infer exploratory cell-cell communication
#'
#' @param obj A `scdown_obj`.
#' @param min_detected Minimum detection rate for ligand and receptor.
#' @param min_specificity Minimum positive enrichment over the rest.
#' @param exclude_pattern Pattern for labels to exclude from communication.
#' @param outdir Optional output directory.
#'
#' @return A list with communication tables and plots.
#' @export
infer_communication <- function(obj,
                                min_detected = 0.1,
                                min_specificity = 0.15,
                                exclude_pattern = "(?i)unassigned|unknown|doublet",
                                outdir = NULL) {
  stopifnot(inherits(obj, "scdown_obj"))
  if (!obj$lineage_col %in% names(obj$cells)) {
    obj <- collapse_lineage(obj)
  }
  expr <- .analysis_expr(obj)
  ann <- obj$annotation_col
  if (!is.null(exclude_pattern) && nzchar(exclude_pattern)) {
    keep <- !grepl(exclude_pattern, expr[[ann]], perl = TRUE)
    expr <- expr[keep, , drop = FALSE]
  }
  pairs <- .lr_pairs()
  cells_by_type <- .cell_type_sizes(expr, ann)
  names(cells_by_type)[2] <- "n_cells_type"
  n_cells_total <- sum(cells_by_type$n_cells_type)
  sum_type <- stats::aggregate(analysis_value ~ cell_type + gene, data = setNames(expr[, c(ann, "gene", "analysis_value")], c("cell_type", "gene", "analysis_value")), FUN = sum)
  detect_type <- stats::aggregate(detected ~ cell_type + gene, data = setNames(expr[, c(ann, "gene", "detected")], c("cell_type", "gene", "detected")), FUN = sum)
  total_sum <- stats::aggregate(analysis_value ~ gene, data = expr, FUN = sum)
  total_detect <- stats::aggregate(detected ~ gene, data = expr, FUN = sum)
  names(sum_type)[3] <- "sum_in_type"
  names(detect_type)[3] <- "detected_in_type"
  names(total_sum)[2] <- "sum_total"
  names(total_detect)[2] <- "detected_total"
  avg <- merge(sum_type, detect_type, by = c("cell_type", "gene"))
  avg <- merge(avg, cells_by_type, by = "cell_type")
  avg <- merge(avg, total_sum, by = "gene")
  avg <- merge(avg, total_detect, by = "gene")
  avg$mean_expression <- avg$sum_in_type / avg$n_cells_type
  avg$mean_rest <- (avg$sum_total - avg$sum_in_type) / pmax(n_cells_total - avg$n_cells_type, 1)
  avg$pct_detected <- avg$detected_in_type / avg$n_cells_type
  avg$pct_rest <- (avg$detected_total - avg$detected_in_type) / pmax(n_cells_total - avg$n_cells_type, 1)
  avg$specificity <- pmax(avg$mean_expression - avg$mean_rest, 0)
  sender <- merge(avg, pairs, by.x = "gene", by.y = "ligand")
  receiver <- merge(avg, pairs, by.x = "gene", by.y = "receptor")
  communication_table <- merge(
    sender[, c("cell_type", "pathway", "gene", "receptor", "mean_expression", "pct_detected", "specificity")],
    receiver[, c("cell_type", "pathway", "ligand", "gene", "mean_expression", "pct_detected", "specificity")],
    by = "pathway",
    suffixes = c("_sender", "_receiver")
  )
  names(communication_table)[names(communication_table) == "cell_type_sender"] <- "sender"
  names(communication_table)[names(communication_table) == "cell_type_receiver"] <- "receiver"
  communication_table <- communication_table[
    communication_table$pct_detected_sender >= min_detected &
      communication_table$pct_detected_receiver >= min_detected &
      communication_table$specificity_sender >= min_specificity &
      communication_table$specificity_receiver >= min_specificity,
    ,
    drop = FALSE
  ]
  communication_table$score <- communication_table$specificity_sender *
    communication_table$specificity_receiver *
    communication_table$pct_detected_sender *
    communication_table$pct_detected_receiver
  communication_table <- communication_table[order(-communication_table$score), c("pathway", "sender", "receiver", "gene_sender", "gene_receiver", "mean_expression_sender", "mean_expression_receiver", "pct_detected_sender", "pct_detected_receiver", "specificity_sender", "specificity_receiver", "score")]
  names(communication_table)[4:5] <- c("ligand", "receptor")
  if (nrow(communication_table)) {
    communication_pair_summary <- stats::aggregate(score ~ sender + receiver, data = communication_table, FUN = sum)
    lineage_lookup <- unique(obj$cells[, c(obj$annotation_col, obj$lineage_col), drop = FALSE])
    names(lineage_lookup) <- c("cell_type", "lineage")
    communication_lineage_summary <- merge(communication_pair_summary, lineage_lookup, by.x = "sender", by.y = "cell_type", all.x = TRUE)
    names(communication_lineage_summary)[names(communication_lineage_summary) == "lineage"] <- "sender_lineage"
    communication_lineage_summary <- merge(communication_lineage_summary, lineage_lookup, by.x = "receiver", by.y = "cell_type", all.x = TRUE)
    names(communication_lineage_summary)[names(communication_lineage_summary) == "lineage"] <- "receiver_lineage"
    communication_lineage_summary <- stats::aggregate(score ~ sender_lineage + receiver_lineage, data = communication_lineage_summary, FUN = sum)
  } else {
    communication_pair_summary <- data.frame(sender = character(), receiver = character(), score = numeric(), stringsAsFactors = FALSE)
    communication_lineage_summary <- data.frame(sender_lineage = character(), receiver_lineage = character(), score = numeric(), stringsAsFactors = FALSE)
  }
  tables <- list(
    communication_table = communication_table,
    communication_pair_summary = communication_pair_summary,
    communication_lineage_summary = communication_lineage_summary
  )
  .write_tables(tables, outdir)
  plotters <- list(
    communication_network = function() {
      df <- communication_pair_summary
      if (!nrow(df)) {
        return(.plot_message("No communication pairs"))
      }
      graphics::barplot(df$score, names.arg = paste(df$sender, df$receiver, sep = "->"), las = 2, col = "goldenrod", main = "Communication pairs", ylab = "Score")
    }
  )
  list(tables = tables, plots = .save_plot_paths(plotters, outdir))
}

#' Run the core downstream workflow
#'
#' @param obj A `scdown_obj`.
#' @param genes User-selected genes.
#' @param signatures Signatures to score.
#' @param outdir Optional output directory.
#'
#' @return A named list of results.
#' @export
run_core <- function(obj,
                     genes = c("CD3D", "NKG7", "LYZ", "MS4A1"),
                     signatures = c("cytotoxic", "exhaustion", "antigen_presentation"),
                     outdir = NULL) {
  stopifnot(inherits(obj, "scdown_obj"))
  if (!is.null(outdir)) .dir_create(outdir)
  obj <- collapse_lineage(obj)
  list(
    summary = summarize_dataset(obj, outdir = if (is.null(outdir)) NULL else file.path(outdir, "summary")),
    map = plot_map(obj, outdir = if (is.null(outdir)) NULL else file.path(outdir, "map")),
    markers = find_markers(obj, outdir = if (is.null(outdir)) NULL else file.path(outdir, "markers")),
    genes = plot_gene(obj, genes = genes, outdir = if (is.null(outdir)) NULL else file.path(outdir, "genes")),
    composition = plot_composition(obj, outdir = if (is.null(outdir)) NULL else file.path(outdir, "composition")),
    annotation = check_annotation(obj, outdir = if (is.null(outdir)) NULL else file.path(outdir, "annotation_check")),
    signatures = plot_signature(obj, signatures = signatures, outdir = if (is.null(outdir)) NULL else file.path(outdir, "signatures")),
    heatmap = plot_heatmap(obj, signatures = signatures, outdir = if (is.null(outdir)) NULL else file.path(outdir, "heatmap")),
    communication = infer_communication(obj, outdir = if (is.null(outdir)) NULL else file.path(outdir, "communication"))
  )
}
