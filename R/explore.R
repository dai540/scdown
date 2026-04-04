.marker_statistics <- function(obj, features = rownames(obj$expr), exclude_pattern = NULL) {
  ann <- as.character(obj$cells[[obj$annotation_col]])
  mat <- .analysis_matrix(obj)
  det <- .detected_matrix(obj)
  features <- intersect(features, rownames(mat))
  if (!is.null(exclude_pattern) && nzchar(exclude_pattern)) {
    features <- features[!grepl(exclude_pattern, features, ignore.case = TRUE)]
  }
  mat <- mat[features, , drop = FALSE]
  det <- det[features, , drop = FALSE]
  cell_types <- unique(ann)
  out <- lapply(cell_types, function(cell_type) {
    idx <- ann == cell_type
    rest <- !idx
    mean_in_type <- .matrix_row_means(mat[, idx, drop = FALSE])
    mean_rest <- if (any(rest)) .matrix_row_means(mat[, rest, drop = FALSE]) else rep(0, length(features))
    pct_in_type <- .matrix_row_means(det[, idx, drop = FALSE])
    pct_rest <- if (any(rest)) .matrix_row_means(det[, rest, drop = FALSE]) else rep(0, length(features))
    data.frame(
      cell_type = cell_type,
      gene = features,
      mean_in_type = mean_in_type,
      mean_rest = mean_rest,
      pct_in_type = pct_in_type,
      pct_rest = pct_rest,
      marker_score = mean_in_type - mean_rest,
      marker_table_rank = (mean_in_type - mean_rest) + 0.5 * pmax(pct_in_type - pct_rest, 0),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out[order(out$cell_type, -out$marker_table_rank, -out$pct_in_type), , drop = FALSE]
}

#' Summarize dataset structure
#'
#' @param obj A `scdown_obj`.
#' @param outdir Optional output directory.
#'
#' @return A list with summary tables.
#' @export
summarize_dataset <- function(obj, outdir = NULL) {
  stopifnot(inherits(obj, "scdown_obj"))
  if (!obj$lineage_col %in% names(obj$cells)) {
    obj <- collapse_lineage(obj)
  }
  dataset_summary <- data.frame(
    metric = c("cells", "samples", "cell_types", "genes", "value_type"),
    value = c(
      ncol(obj$expr),
      length(unique(obj$cells[[obj$sample_col]])),
      length(unique(obj$cells[[obj$annotation_col]])),
      nrow(obj$expr),
      obj$value_type
    ),
    stringsAsFactors = FALSE
  )
  celltype_count <- as.data.frame(table(obj$cells[[obj$annotation_col]]), stringsAsFactors = FALSE)
  names(celltype_count) <- c("cell_type", "n_cells")
  lineage_count <- as.data.frame(table(obj$cells[[obj$lineage_col]]), stringsAsFactors = FALSE)
  names(lineage_count) <- c("lineage", "n_cells")
  sample_summary <- as.data.frame(table(obj$cells[[obj$sample_col]]), stringsAsFactors = FALSE)
  names(sample_summary) <- c("sample", "n_cells")
  tables <- list(
    dataset_summary = dataset_summary,
    celltype_count = celltype_count,
    lineage_count = lineage_count,
    sample_summary = sample_summary
  )
  if (obj$group_col %in% names(obj$cells)) {
    group_summary <- as.data.frame(table(obj$cells[[obj$group_col]]), stringsAsFactors = FALSE)
    names(group_summary) <- c("group", "n_cells")
    tables$group_summary <- group_summary
  }
  .write_tables(tables, outdir)
  list(tables = tables, plots = character())
}

#' Collapse fine labels into major lineages
#'
#' @param obj A `scdown_obj`.
#' @param mapping Optional data frame with `pattern` and `lineage`.
#'
#' @return The updated `scdown_obj`.
#' @export
collapse_lineage <- function(obj, mapping = NULL) {
  stopifnot(inherits(obj, "scdown_obj"))
  mapping <- .resource_lineage_mapping(obj, override = mapping)
  labels <- as.character(obj$cells[[obj$annotation_col]])
  lineage <- rep("Other", length(labels))
  for (i in seq_len(nrow(mapping))) {
    hit <- grepl(mapping$pattern[[i]], labels, ignore.case = TRUE)
    lineage[hit] <- mapping$lineage[[i]]
  }
  obj$cells[[obj$lineage_col]] <- lineage
  obj
}

#' Plot cell map
#'
#' @param obj A `scdown_obj`.
#' @param outdir Optional output directory.
#'
#' @return A list with tables and plot paths.
#' @export
plot_map <- function(obj, outdir = NULL) {
  stopifnot(inherits(obj, "scdown_obj"))
  if (!obj$lineage_col %in% names(obj$cells)) {
    obj <- collapse_lineage(obj)
  }
  embedding_cells <- merge(obj$embedding, obj$cells, by = "cell", all.x = TRUE, sort = FALSE)
  tables <- list(embedding_cells = embedding_cells)
  .write_tables(tables, outdir)
  plotters <- list(
    map_by_celltype = function() {
      info <- .factor_info(embedding_cells[[obj$annotation_col]])
      cols <- grDevices::hcl.colors(length(info$levels), "Set 2")
      graphics::plot(embedding_cells$UMAP1, embedding_cells$UMAP2, col = cols[info$index], pch = 19, xlab = "Dim 1", ylab = "Dim 2", main = "Cells by annotation")
      graphics::legend("topright", legend = info$levels, col = cols, pch = 19, cex = 0.8)
    },
    map_by_lineage = function() {
      info <- .factor_info(embedding_cells[[obj$lineage_col]])
      cols <- grDevices::hcl.colors(length(info$levels), "Dark 3")
      graphics::plot(embedding_cells$UMAP1, embedding_cells$UMAP2, col = cols[info$index], pch = 19, xlab = "Dim 1", ylab = "Dim 2", main = "Cells by lineage")
      graphics::legend("topright", legend = info$levels, col = cols, pch = 19, cex = 0.8)
    }
  )
  if (obj$group_col %in% names(embedding_cells)) {
    plotters$map_by_group <- function() {
      info <- .factor_info(embedding_cells[[obj$group_col]])
      cols <- grDevices::hcl.colors(length(info$levels), "Dynamic")
      graphics::plot(embedding_cells$UMAP1, embedding_cells$UMAP2, col = cols[info$index], pch = 19, xlab = "Dim 1", ylab = "Dim 2", main = "Cells by group")
      graphics::legend("topright", legend = info$levels, col = cols, pch = 19, cex = 0.8)
    }
  }
  list(tables = tables, plots = .save_plot_paths(plotters, outdir))
}

#' Explore cell-type markers
#'
#' @param obj A `scdown_obj`.
#' @param top_n Number of markers per cell type.
#' @param exclude_pattern Pattern for low-information genes to exclude.
#' @param min_pct_in_type Minimum detection rate inside the target cell type.
#' @param min_pct_diff Minimum detection-rate improvement over the rest.
#' @param min_marker_score Minimum mean-expression improvement over the rest.
#' @param outdir Optional output directory.
#'
#' @return A list with marker tables and plot paths.
#' @export
explore_markers <- function(obj,
                            top_n = 10L,
                            exclude_pattern = "^(RPL|RPS|MT-|MALAT1$|MTRNR|HSP)",
                            min_pct_in_type = 0.1,
                            min_pct_diff = 0.05,
                            min_marker_score = 0.25,
                            outdir = NULL) {
  stopifnot(inherits(obj, "scdown_obj"))
  marker_table <- .marker_statistics(obj, exclude_pattern = exclude_pattern)
  keep <- marker_table$pct_in_type >= min_pct_in_type &
    (marker_table$pct_in_type - marker_table$pct_rest) >= min_pct_diff &
    marker_table$marker_score >= min_marker_score
  marker_hits <- marker_table[keep, , drop = FALSE]
  if (!nrow(marker_hits)) marker_hits <- marker_table[marker_table$marker_score > 0, , drop = FALSE]
  if (!nrow(marker_hits)) marker_hits <- marker_table
  top_marker_per_celltype <- .top_n_per_group(marker_hits, "cell_type", "marker_table_rank", n = top_n)
  tables <- list(marker_table = marker_table, top_marker_per_celltype = top_marker_per_celltype)
  .write_tables(tables, outdir)
  plotters <- list(
    marker_dotplot = function() {
      top <- utils::head(top_marker_per_celltype, 30)
      cell_info <- .factor_info(top$cell_type)
      gene_info <- .factor_info(top$gene)
      graphics::plot(cell_info$index, gene_info$index, cex = pmax(top$marker_score, 0.2), pch = 19, col = "steelblue", xlab = "Cell type", ylab = "Gene", xaxt = "n", yaxt = "n", main = "Top exploratory markers")
      graphics::axis(1, at = seq_along(cell_info$levels), labels = cell_info$levels, las = 2, cex.axis = 0.7)
      graphics::axis(2, at = seq_along(gene_info$levels), labels = gene_info$levels, las = 2, cex.axis = 0.7)
    },
    marker_heatmap = function() {
      mat <- stats::xtabs(marker_score ~ cell_type + gene, data = top_marker_per_celltype)
      graphics::image(seq_len(ncol(mat)), seq_len(nrow(mat)), t(mat[nrow(mat):1, , drop = FALSE]), axes = FALSE, col = grDevices::hcl.colors(20, "YlOrRd"), main = "Exploratory marker heatmap")
      graphics::axis(1, at = seq_len(ncol(mat)), labels = colnames(mat), las = 2, cex.axis = 0.6)
      graphics::axis(2, at = seq_len(nrow(mat)), labels = rev(rownames(mat)), las = 2, cex.axis = 0.7)
    }
  )
  list(tables = tables, plots = .save_plot_paths(plotters, outdir))
}

#' Find cell-type markers
#'
#' @inheritParams explore_markers
#'
#' @return A list with marker tables and plot paths.
#' @export
find_markers <- function(obj,
                         top_n = 10L,
                         exclude_pattern = "^(RPL|RPS|MT-|MALAT1$|MTRNR|HSP)",
                         min_pct_in_type = 0.1,
                         min_pct_diff = 0.05,
                         min_marker_score = 0.25,
                         outdir = NULL) {
  explore_markers(
    obj = obj,
    top_n = top_n,
    exclude_pattern = exclude_pattern,
    min_pct_in_type = min_pct_in_type,
    min_pct_diff = min_pct_diff,
    min_marker_score = min_marker_score,
    outdir = outdir
  )
}

#' Plot user-selected genes
#'
#' @param obj A `scdown_obj`.
#' @param genes Character vector of genes to visualize.
#' @param outdir Optional output directory.
#'
#' @return A list with gene tables and plot paths.
#' @export
plot_gene <- function(obj, genes, outdir = NULL) {
  stopifnot(inherits(obj, "scdown_obj"))
  genes <- unique(as.character(genes))
  raw_mat <- obj$expr
  analysis_mat <- .analysis_matrix(obj)
  gene_expression_long <- .matrix_to_long_subset(raw_mat, analysis_mat, genes, obj$cells, obj$sample_col, obj$annotation_col)
  missing_genes <- data.frame(gene = setdiff(genes, rownames(raw_mat)), stringsAsFactors = FALSE)
  gene_average_by_celltype <- if (nrow(gene_expression_long)) {
    stats::aggregate(analysis_value ~ gene + cell_type, data = gene_expression_long, FUN = mean)
  } else {
    data.frame(gene = character(), cell_type = character(), analysis_value = numeric(), stringsAsFactors = FALSE)
  }
  names(gene_average_by_celltype)[3] <- "mean_expression"
  gene_average_by_sample <- if (nrow(gene_expression_long)) {
    stats::aggregate(analysis_value ~ gene + sample, data = gene_expression_long, FUN = mean)
  } else {
    data.frame(gene = character(), sample = character(), analysis_value = numeric(), stringsAsFactors = FALSE)
  }
  names(gene_average_by_sample)[3] <- "mean_expression"
  tables <- list(
    gene_expression_long = gene_expression_long,
    gene_average_by_celltype = gene_average_by_celltype,
    gene_average_by_sample = gene_average_by_sample,
    missing_genes = missing_genes
  )
  .write_tables(tables, outdir)
  plotters <- list()
  for (gene in intersect(genes, rownames(analysis_mat))) {
    plotters[[paste0("featureplot_", gene)]] <- local({
      g <- gene
      function() {
        values <- as.numeric(analysis_mat[g, ])
        df <- cbind(obj$embedding, value = values)
        pal <- grDevices::hcl.colors(20, "YlGnBu")
        bins <- cut(df$value, breaks = 20, include.lowest = TRUE)
        graphics::plot(df$UMAP1, df$UMAP2, col = pal[bins], pch = 19, xlab = "Dim 1", ylab = "Dim 2", main = paste("Expression:", g))
      }
    })
  }
  if (nrow(gene_average_by_celltype)) {
    plotters$dotplot_selected_genes <- function() {
      df <- gene_average_by_celltype
      cell_info <- .factor_info(df$cell_type)
      gene_info <- .factor_info(df$gene)
      graphics::plot(cell_info$index, gene_info$index, cex = pmax(df$mean_expression, 0.3), pch = 19, col = "navy", xlab = "Cell type", ylab = "Gene", xaxt = "n", yaxt = "n", main = "Selected genes")
      graphics::axis(1, at = seq_along(cell_info$levels), labels = cell_info$levels, las = 2, cex.axis = 0.7)
      graphics::axis(2, at = seq_along(gene_info$levels), labels = gene_info$levels, las = 2, cex.axis = 0.8)
    }
  }
  list(tables = tables, plots = .save_plot_paths(plotters, outdir))
}

.composition_tables <- function(obj) {
  comp <- as.data.frame(table(obj$cells[[obj$sample_col]], obj$cells[[obj$annotation_col]]), stringsAsFactors = FALSE)
  names(comp) <- c("sample", "cell_type", "n_cells")
  totals <- stats::aggregate(n_cells ~ sample, data = comp, FUN = sum)
  names(totals)[2] <- "n_total"
  composition_by_sample <- merge(comp, totals, by = "sample")
  composition_by_sample$fraction <- composition_by_sample$n_cells / composition_by_sample$n_total

  lineage_comp <- as.data.frame(table(obj$cells[[obj$sample_col]], obj$cells[[obj$lineage_col]]), stringsAsFactors = FALSE)
  names(lineage_comp) <- c("sample", "lineage", "n_cells")
  lineage_totals <- stats::aggregate(n_cells ~ sample, data = lineage_comp, FUN = sum)
  names(lineage_totals)[2] <- "n_total"
  composition_major_lineage <- merge(lineage_comp, lineage_totals, by = "sample")
  composition_major_lineage$fraction <- composition_major_lineage$n_cells / composition_major_lineage$n_total

  sample_meta <- unique(obj$cells[, .existing_columns(obj$cells, c(obj$sample_col, obj$group_col)), drop = FALSE])
  names(sample_meta)[names(sample_meta) == obj$sample_col] <- "sample"
  if ("sample" %in% names(sample_meta) && obj$group_col %in% names(obj$cells)) {
    names(sample_meta)[names(sample_meta) == obj$group_col] <- "group"
    composition_by_sample <- merge(composition_by_sample, sample_meta, by = "sample", all.x = TRUE)
    composition_major_lineage <- merge(composition_major_lineage, sample_meta, by = "sample", all.x = TRUE)
    composition_by_group <- stats::aggregate(fraction ~ group + cell_type, data = composition_by_sample, FUN = mean)
  } else {
    composition_by_group <- data.frame(group = character(), cell_type = character(), fraction = numeric(), stringsAsFactors = FALSE)
  }

  list(
    composition_by_sample = composition_by_sample,
    composition_major_lineage = composition_major_lineage,
    composition_by_group = composition_by_group
  )
}

#' Plot composition summaries
#'
#' @param obj A `scdown_obj`.
#' @param outdir Optional output directory.
#'
#' @return A list with composition tables and plot paths.
#' @export
plot_composition <- function(obj, outdir = NULL) {
  stopifnot(inherits(obj, "scdown_obj"))
  if (!obj$lineage_col %in% names(obj$cells)) {
    obj <- collapse_lineage(obj)
  }
  tables <- .composition_tables(obj)
  .write_tables(tables, outdir)
  plotters <- list(
    composition_stacked = function() {
      mat <- stats::xtabs(fraction ~ sample + cell_type, data = tables$composition_by_sample)
      graphics::barplot(t(mat), beside = FALSE, col = grDevices::hcl.colors(ncol(mat), "Set 3"), legend.text = colnames(mat), args.legend = list(x = "topright", cex = 0.8), main = "Cell-type composition", ylab = "Fraction")
    },
    composition_major_lineage = function() {
      mat <- stats::xtabs(fraction ~ sample + lineage, data = tables$composition_major_lineage)
      graphics::barplot(t(mat), beside = FALSE, col = grDevices::hcl.colors(ncol(mat), "Dark 2"), legend.text = colnames(mat), args.legend = list(x = "topright", cex = 0.8), main = "Major lineage composition", ylab = "Fraction")
    }
  )
  if (nrow(tables$composition_by_group)) {
    plotters$composition_boxplot <- function() {
      df <- tables$composition_by_sample
      graphics::boxplot(fraction ~ interaction(group, cell_type, lex.order = TRUE), data = df, las = 2, col = "grey80", main = "Sample-level composition by group", ylab = "Fraction")
    }
  }
  list(tables = tables, plots = .save_plot_paths(plotters, outdir))
}

#' Explore composition summaries
#'
#' @inheritParams plot_composition
#'
#' @return A list with composition tables and plot paths.
#' @export
explore_composition <- function(obj, outdir = NULL) {
  plot_composition(obj = obj, outdir = outdir)
}

.score_signature <- function(obj, genes, score_name) {
  mat <- .analysis_matrix(obj)
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
#' @param panels Character vector of panel names or a named list of custom panels.
#' @param marker_panels Optional named list of custom marker panels to merge with built-ins.
#' @param outdir Optional output directory.
#'
#' @return A list with annotation-check tables and plots.
#' @export
check_annotation <- function(obj, panels = NULL, marker_panels = NULL, outdir = NULL) {
  stopifnot(inherits(obj, "scdown_obj"))
  panels <- .resolve_gene_sets(obj, panels, resource = "marker_panels", override = marker_panels)
  panel_names <- names(panels)
  cell_scores <- unique(obj$cells["cell"])
  for (nm in panel_names) {
    cell_scores <- merge(cell_scores, .score_signature(obj, panels[[nm]], nm), by = "cell", all.x = TRUE)
  }
  cell_scores[is.na(cell_scores)] <- 0
  long <- stats::reshape(
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
    },
    annotation_check_dotplot = function() {
      df <- annotation_check_score
      cell_info <- .factor_info(df$cell_type)
      panel_info <- .factor_info(df$panel)
      graphics::plot(cell_info$index, panel_info$index, cex = pmax(df$score, 0.3), pch = 19, col = "firebrick", xlab = "Cell type", ylab = "Panel", xaxt = "n", yaxt = "n", main = "Annotation check")
      graphics::axis(1, at = seq_along(cell_info$levels), labels = cell_info$levels, las = 2, cex.axis = 0.7)
      graphics::axis(2, at = seq_along(panel_info$levels), labels = panel_info$levels, las = 2, cex.axis = 0.8)
    }
  )
  list(tables = tables, plots = .save_plot_paths(plotters, outdir))
}

#' Plot signature scores
#'
#' @param obj A `scdown_obj`.
#' @param signatures Character vector of signature names or a named list of signatures.
#' @param signature_sets Optional named list of custom signatures to merge with built-ins.
#' @param outdir Optional output directory.
#'
#' @return A list with signature tables and plots.
#' @export
plot_signature <- function(obj, signatures, signature_sets = NULL, outdir = NULL) {
  stopifnot(inherits(obj, "scdown_obj"))
  signatures <- .resolve_gene_sets(obj, signatures, resource = "signatures", override = signature_sets)
  sig_names <- names(signatures)
  cell_scores <- unique(obj$cells["cell"])
  for (sig in sig_names) {
    cell_scores <- merge(cell_scores, .score_signature(obj, signatures[[sig]], sig), by = "cell", all.x = TRUE)
  }
  cell_scores[is.na(cell_scores)] <- 0
  score_long <- stats::reshape(
    merge(cell_scores, obj$cells[, c("cell", obj$annotation_col), drop = FALSE], by = "cell", all.x = TRUE),
    varying = sig_names,
    v.names = "score",
    timevar = "signature",
    times = sig_names,
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

#' Explore signature scores
#'
#' @inheritParams plot_signature
#'
#' @return A list with signature tables and plots.
#' @export
explore_signature <- function(obj, signatures, signature_sets = NULL, outdir = NULL) {
  plot_signature(obj = obj, signatures = signatures, signature_sets = signature_sets, outdir = outdir)
}

#' Plot expression heatmaps
#'
#' @param obj A `scdown_obj`.
#' @param signatures Signatures to include in the signature heatmap.
#' @param signature_sets Optional custom signatures.
#' @param top_n_markers Number of markers per cell type to show.
#' @param outdir Optional output directory.
#'
#' @return A list with heatmap tables and plots.
#' @export
plot_heatmap <- function(obj,
                         signatures = c("cytotoxic", "exhaustion", "antigen_presentation"),
                         signature_sets = NULL,
                         top_n_markers = 5L,
                         outdir = NULL) {
  stopifnot(inherits(obj, "scdown_obj"))
  mat <- .analysis_matrix(obj)
  ann <- obj$cells[[obj$annotation_col]]
  average_expression <- do.call(rbind, lapply(unique(ann), function(cell_type) {
    idx <- ann == cell_type
    data.frame(
      cell_type = cell_type,
      gene = rownames(mat),
      mean_expression = .matrix_row_means(mat[, idx, drop = FALSE]),
      stringsAsFactors = FALSE
    )
  }))
  marker_res <- explore_markers(obj, top_n = top_n_markers)
  marker_heatmap_matrix <- marker_res$tables$top_marker_per_celltype
  signature_res <- plot_signature(obj, signatures = signatures, signature_sets = signature_sets)
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
      mat_plot <- stats::xtabs(mean_expression ~ cell_type + gene, data = df)
      graphics::image(seq_len(ncol(mat_plot)), seq_len(nrow(mat_plot)), t(mat_plot[nrow(mat_plot):1, , drop = FALSE]), axes = FALSE, col = grDevices::hcl.colors(20, "Teal"), main = "Average expression")
      graphics::axis(1, at = seq_len(ncol(mat_plot)), labels = colnames(mat_plot), las = 2, cex.axis = 0.6)
      graphics::axis(2, at = seq_len(nrow(mat_plot)), labels = rev(rownames(mat_plot)), las = 2, cex.axis = 0.7)
    },
    marker_heatmap = function() {
      mat_plot <- stats::xtabs(marker_score ~ cell_type + gene, data = marker_heatmap_matrix)
      graphics::image(seq_len(ncol(mat_plot)), seq_len(nrow(mat_plot)), t(mat_plot[nrow(mat_plot):1, , drop = FALSE]), axes = FALSE, col = grDevices::hcl.colors(20, "Inferno"), main = "Marker heatmap")
      graphics::axis(1, at = seq_len(ncol(mat_plot)), labels = colnames(mat_plot), las = 2, cex.axis = 0.6)
      graphics::axis(2, at = seq_len(nrow(mat_plot)), labels = rev(rownames(mat_plot)), las = 2, cex.axis = 0.7)
    },
    signature_heatmap = function() {
      mat_plot <- stats::xtabs(score ~ cell_type + signature, data = signature_heatmap_matrix)
      graphics::image(seq_len(ncol(mat_plot)), seq_len(nrow(mat_plot)), t(mat_plot[nrow(mat_plot):1, , drop = FALSE]), axes = FALSE, col = grDevices::hcl.colors(20, "Viridis"), main = "Signature heatmap")
      graphics::axis(1, at = seq_len(ncol(mat_plot)), labels = colnames(mat_plot), las = 2, cex.axis = 0.7)
      graphics::axis(2, at = seq_len(nrow(mat_plot)), labels = rev(rownames(mat_plot)), las = 2, cex.axis = 0.7)
    }
  )
  list(tables = tables, plots = .save_plot_paths(plotters, outdir))
}

.communication_tables <- function(obj,
                                  labels = NULL,
                                  lr_pairs = NULL,
                                  min_detected = 0.1,
                                  min_specificity = 0.15,
                                  exclude_pattern = "(?i)unassigned|unknown|doublet") {
  if (!obj$lineage_col %in% names(obj$cells)) {
    obj <- collapse_lineage(obj)
  }
  labels <- as.character(labels %||% obj$cells[[obj$annotation_col]])
  keep_cells <- rep(TRUE, length(labels))
  if (!is.null(exclude_pattern) && nzchar(exclude_pattern)) {
    keep_cells <- !grepl(exclude_pattern, labels, perl = TRUE)
  }
  labels <- labels[keep_cells]
  empty_table <- data.frame(
    pathway = character(), sender = character(), receiver = character(),
    ligand = character(), receptor = character(),
    mean_expression_sender = numeric(), mean_expression_receiver = numeric(),
    pct_detected_sender = numeric(), pct_detected_receiver = numeric(),
    specificity_sender = numeric(), specificity_receiver = numeric(),
    score = numeric(), stringsAsFactors = FALSE
  )
  if (!length(labels)) {
    return(list(
      communication_table = empty_table,
      communication_pair_summary = data.frame(sender = character(), receiver = character(), score = numeric(), stringsAsFactors = FALSE),
      communication_lineage_summary = data.frame(sender_lineage = character(), receiver_lineage = character(), score = numeric(), stringsAsFactors = FALSE)
    ))
  }
  pairs <- .resource_lr_pairs(obj, override = lr_pairs)
  genes <- intersect(unique(c(pairs$ligand, pairs$receptor)), rownames(obj$expr))
  if (!length(genes)) {
    return(list(
      communication_table = empty_table,
      communication_pair_summary = data.frame(sender = character(), receiver = character(), score = numeric(), stringsAsFactors = FALSE),
      communication_lineage_summary = data.frame(sender_lineage = character(), receiver_lineage = character(), score = numeric(), stringsAsFactors = FALSE)
    ))
  }
  mat <- .analysis_matrix(obj)[genes, keep_cells, drop = FALSE]
  det <- .detected_matrix(obj)[genes, keep_cells, drop = FALSE]
  cell_types <- unique(labels)
  avg <- do.call(rbind, lapply(cell_types, function(cell_type) {
    idx <- labels == cell_type
    rest <- !idx
    mean_expression <- .matrix_row_means(mat[, idx, drop = FALSE])
    mean_rest <- if (any(rest)) .matrix_row_means(mat[, rest, drop = FALSE]) else rep(0, length(genes))
    pct_detected <- .matrix_row_means(det[, idx, drop = FALSE])
    data.frame(
      cell_type = cell_type,
      gene = genes,
      mean_expression = mean_expression,
      pct_detected = pct_detected,
      specificity = pmax(mean_expression - mean_rest, 0),
      stringsAsFactors = FALSE
    )
  }))
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
  rownames(communication_table) <- NULL
  if (!nrow(communication_table)) {
    communication_pair_summary <- data.frame(sender = character(), receiver = character(), score = numeric(), stringsAsFactors = FALSE)
    communication_lineage_summary <- data.frame(sender_lineage = character(), receiver_lineage = character(), score = numeric(), stringsAsFactors = FALSE)
  } else {
    communication_pair_summary <- stats::aggregate(score ~ sender + receiver, data = communication_table, FUN = sum)
    lookup <- unique(obj$cells[, c(obj$annotation_col, obj$lineage_col), drop = FALSE])
    names(lookup) <- c("cell_type", "lineage")
    communication_lineage_summary <- merge(communication_pair_summary, lookup, by.x = "sender", by.y = "cell_type", all.x = TRUE)
    names(communication_lineage_summary)[names(communication_lineage_summary) == "lineage"] <- "sender_lineage"
    communication_lineage_summary <- merge(communication_lineage_summary, lookup, by.x = "receiver", by.y = "cell_type", all.x = TRUE)
    names(communication_lineage_summary)[names(communication_lineage_summary) == "lineage"] <- "receiver_lineage"
    communication_lineage_summary <- stats::aggregate(score ~ sender_lineage + receiver_lineage, data = communication_lineage_summary, FUN = sum)
  }
  list(
    communication_table = communication_table,
    communication_pair_summary = communication_pair_summary,
    communication_lineage_summary = communication_lineage_summary
  )
}

#' Infer exploratory cell-cell communication
#'
#' @param obj A `scdown_obj`.
#' @param lr_pairs Optional data frame with custom ligand-receptor pairs.
#' @param min_detected Minimum detection rate for ligand and receptor.
#' @param min_specificity Minimum positive enrichment over the rest.
#' @param exclude_pattern Pattern for labels to exclude from communication.
#' @param outdir Optional output directory.
#'
#' @return A list with communication tables and plots.
#' @export
infer_communication <- function(obj,
                                lr_pairs = NULL,
                                min_detected = 0.1,
                                min_specificity = 0.15,
                                exclude_pattern = "(?i)unassigned|unknown|doublet",
                                outdir = NULL) {
  stopifnot(inherits(obj, "scdown_obj"))
  tables <- .communication_tables(
    obj = obj,
    lr_pairs = lr_pairs,
    min_detected = min_detected,
    min_specificity = min_specificity,
    exclude_pattern = exclude_pattern
  )
  .write_tables(tables, outdir)
  plotters <- list(
    communication_network = function() {
      df <- tables$communication_pair_summary
      if (!nrow(df)) {
        return(.plot_message("No communication pairs"))
      }
      graphics::barplot(df$score, names.arg = paste(df$sender, df$receiver, sep = "->"), las = 2, col = "goldenrod", main = "Communication pairs", ylab = "Score")
    },
    communication_heatmap = function() {
      df <- tables$communication_lineage_summary
      if (!nrow(df)) {
        return(.plot_message("No lineage summary"))
      }
      mat <- stats::xtabs(score ~ sender_lineage + receiver_lineage, data = df)
      graphics::image(seq_len(ncol(mat)), seq_len(nrow(mat)), t(mat[nrow(mat):1, , drop = FALSE]), axes = FALSE, col = grDevices::hcl.colors(20, "Sunset"), main = "Communication by lineage")
      graphics::axis(1, at = seq_len(ncol(mat)), labels = colnames(mat), las = 2)
      graphics::axis(2, at = seq_len(nrow(mat)), labels = rev(rownames(mat)), las = 2)
    }
  )
  list(tables = tables, plots = .save_plot_paths(plotters, outdir))
}

#' Explore cell-cell communication
#'
#' @inheritParams infer_communication
#'
#' @return A list with communication tables and plots.
#' @export
explore_communication <- function(obj,
                                  lr_pairs = NULL,
                                  min_detected = 0.1,
                                  min_specificity = 0.15,
                                  exclude_pattern = "(?i)unassigned|unknown|doublet",
                                  outdir = NULL) {
  infer_communication(
    obj = obj,
    lr_pairs = lr_pairs,
    min_detected = min_detected,
    min_specificity = min_specificity,
    exclude_pattern = exclude_pattern,
    outdir = outdir
  )
}

#' Run the core downstream workflow
#'
#' @param obj A `scdown_obj`.
#' @param genes User-selected genes.
#' @param signatures Signatures to score.
#' @param signature_sets Optional custom signatures.
#' @param outdir Optional output directory.
#'
#' @return A named list of results.
#' @export
run_core <- function(obj,
                     genes = c("CD3D", "NKG7", "LYZ", "MS4A1"),
                     signatures = c("cytotoxic", "exhaustion", "antigen_presentation"),
                     signature_sets = NULL,
                     outdir = NULL) {
  stopifnot(inherits(obj, "scdown_obj"))
  if (!is.null(outdir)) .dir_create(outdir)
  obj <- collapse_lineage(obj)
  list(
    summary = summarize_dataset(obj, outdir = if (is.null(outdir)) NULL else file.path(outdir, "summary")),
    map = plot_map(obj, outdir = if (is.null(outdir)) NULL else file.path(outdir, "map")),
    markers = explore_markers(obj, outdir = if (is.null(outdir)) NULL else file.path(outdir, "markers")),
    genes = plot_gene(obj, genes = genes, outdir = if (is.null(outdir)) NULL else file.path(outdir, "genes")),
    composition = plot_composition(obj, outdir = if (is.null(outdir)) NULL else file.path(outdir, "composition")),
    annotation = check_annotation(obj, outdir = if (is.null(outdir)) NULL else file.path(outdir, "annotation_check")),
    signatures = plot_signature(obj, signatures = signatures, signature_sets = signature_sets, outdir = if (is.null(outdir)) NULL else file.path(outdir, "signatures")),
    heatmap = plot_heatmap(obj, signatures = signatures, signature_sets = signature_sets, outdir = if (is.null(outdir)) NULL else file.path(outdir, "heatmap")),
    communication = infer_communication(obj, outdir = if (is.null(outdir)) NULL else file.path(outdir, "communication"))
  )
}
