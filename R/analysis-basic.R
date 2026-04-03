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
    metric = c("cells", "samples", "cell_types", "genes"),
    value = c(
      length(unique(obj$cells$cell)),
      length(unique(obj$cells[[obj$sample_col]])),
      length(unique(obj$cells[[obj$annotation_col]])),
      length(unique(obj$expr$gene))
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
  mapping <- mapping %||% default_lineage_mapping()
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
  list(tables = tables, plots = .save_plot_paths(plotters, outdir))
}

#' Find cell-type markers
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
find_markers <- function(obj,
                         top_n = 10L,
                         exclude_pattern = "^(RPL|RPS|MT-|MALAT1$|MTRNR|HSP)",
                         min_pct_in_type = 0.1,
                         min_pct_diff = 0.05,
                         min_marker_score = 0.25,
                         outdir = NULL) {
  stopifnot(inherits(obj, "scdown_obj"))
  expr <- .analysis_expr(obj)
  ann <- obj$annotation_col
  expr <- expr[!grepl(exclude_pattern, expr$gene, ignore.case = TRUE), , drop = FALSE]
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
  marker_table <- merge(sum_type, detect_type, by = c("cell_type", "gene"))
  marker_table <- merge(marker_table, cells_by_type, by = "cell_type")
  marker_table <- merge(marker_table, total_sum, by = "gene")
  marker_table <- merge(marker_table, total_detect, by = "gene")
  marker_table$mean_in_type <- marker_table$sum_in_type / marker_table$n_cells_type
  marker_table$mean_rest <- (marker_table$sum_total - marker_table$sum_in_type) / pmax(n_cells_total - marker_table$n_cells_type, 1)
  marker_table$pct_in_type <- marker_table$detected_in_type / marker_table$n_cells_type
  marker_table$pct_rest <- (marker_table$detected_total - marker_table$detected_in_type) / pmax(n_cells_total - marker_table$n_cells_type, 1)
  marker_table$marker_score <- marker_table$mean_in_type - marker_table$mean_rest
  marker_table$marker_table_rank <- marker_table$marker_score + 0.5 * pmax(marker_table$pct_in_type - marker_table$pct_rest, 0)
  marker_table <- marker_table[order(marker_table$cell_type, -marker_table$marker_table_rank, -marker_table$pct_in_type), c("cell_type", "gene", "mean_in_type", "mean_rest", "pct_in_type", "pct_rest", "marker_score", "marker_table_rank")]
  rownames(marker_table) <- NULL
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
      graphics::plot(cell_info$index, gene_info$index, cex = pmax(top$marker_score, 0.2), pch = 19, col = "steelblue", xlab = "Cell type", ylab = "Gene", xaxt = "n", yaxt = "n", main = "Top markers")
      graphics::axis(1, at = seq_along(cell_info$levels), labels = cell_info$levels, las = 2, cex.axis = 0.7)
      graphics::axis(2, at = seq_along(gene_info$levels), labels = gene_info$levels, las = 2, cex.axis = 0.7)
    },
    marker_heatmap = function() {
      top <- utils::head(top_marker_per_celltype, top_n * length(unique(top_marker_per_celltype$cell_type)))
      mat <- stats::xtabs(marker_score ~ cell_type + gene, data = top)
      graphics::image(seq_len(ncol(mat)), seq_len(nrow(mat)), t(mat[nrow(mat):1, , drop = FALSE]), axes = FALSE, col = grDevices::hcl.colors(20, "YlOrRd"), main = "Marker heatmap")
      graphics::axis(1, at = seq_len(ncol(mat)), labels = colnames(mat), las = 2, cex.axis = 0.6)
      graphics::axis(2, at = seq_len(nrow(mat)), labels = rev(rownames(mat)), las = 2, cex.axis = 0.7)
    }
  )
  list(tables = tables, plots = .save_plot_paths(plotters, outdir))
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
  expr <- .analysis_expr(obj)
  genes <- unique(as.character(genes))
  expr_gene <- expr[expr$gene %in% genes, , drop = FALSE]
  missing_genes <- data.frame(gene = setdiff(genes, unique(expr$gene)), stringsAsFactors = FALSE)
  gene_average_by_celltype <- stats::aggregate(analysis_value ~ gene + expr_gene[[obj$annotation_col]], data = expr_gene, FUN = mean)
  names(gene_average_by_celltype)[2:3] <- c("cell_type", "mean_expression")
  gene_average_by_sample <- stats::aggregate(analysis_value ~ gene + expr_gene[[obj$sample_col]], data = expr_gene, FUN = mean)
  names(gene_average_by_sample)[2:3] <- c("sample", "mean_expression")
  gene_expression_long <- expr_gene[, c("cell", obj$sample_col, obj$annotation_col, "gene", "value", "analysis_value"), drop = FALSE]
  names(gene_expression_long)[2:3] <- c("sample", "cell_type")
  tables <- list(
    gene_expression_long = gene_expression_long,
    gene_average_by_celltype = gene_average_by_celltype,
    gene_average_by_sample = gene_average_by_sample,
    missing_genes = missing_genes
  )
  .write_tables(tables, outdir)
  plotters <- list()
  for (gene in intersect(genes, unique(expr$gene))) {
    plotters[[paste0("featureplot_", gene)]] <- local({
      g <- gene
      function() {
        df <- merge(obj$embedding, expr[expr$gene == g, c("cell", "analysis_value"), drop = FALSE], by = "cell", all.x = TRUE, sort = FALSE)
        graphics::plot(df$UMAP1, df$UMAP2, col = grDevices::hcl.colors(20, "YlGnBu")[cut(df$analysis_value, breaks = 20, include.lowest = TRUE)], pch = 19, xlab = "Dim 1", ylab = "Dim 2", main = paste("Expression:", g))
      }
    })
  }
  plotters$dotplot_selected_genes <- function() {
    df <- gene_average_by_celltype
    cell_info <- .factor_info(df$cell_type)
    gene_info <- .factor_info(df$gene)
    graphics::plot(cell_info$index, gene_info$index, cex = pmax(df$mean_expression, 0.3), pch = 19, col = "navy", xlab = "Cell type", ylab = "Gene", xaxt = "n", yaxt = "n", main = "Selected genes")
    graphics::axis(1, at = seq_along(cell_info$levels), labels = cell_info$levels, las = 2, cex.axis = 0.7)
    graphics::axis(2, at = seq_along(gene_info$levels), labels = gene_info$levels, las = 2, cex.axis = 0.8)
  }
  list(tables = tables, plots = .save_plot_paths(plotters, outdir))
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
  tables <- list(
    composition_by_sample = composition_by_sample,
    composition_major_lineage = composition_major_lineage
  )
  .write_tables(tables, outdir)
  plotters <- list(
    composition_stacked = function() {
      mat <- stats::xtabs(fraction ~ sample + cell_type, data = composition_by_sample)
      graphics::barplot(t(mat), beside = FALSE, col = grDevices::hcl.colors(ncol(mat), "Set 3"), legend.text = colnames(mat), args.legend = list(x = "topright", cex = 0.8), main = "Cell-type composition", ylab = "Fraction")
    },
    composition_major_lineage = function() {
      mat <- stats::xtabs(fraction ~ sample + lineage, data = composition_major_lineage)
      graphics::barplot(t(mat), beside = FALSE, col = grDevices::hcl.colors(ncol(mat), "Dark 2"), legend.text = colnames(mat), args.legend = list(x = "topright", cex = 0.8), main = "Major lineage composition", ylab = "Fraction")
    }
  )
  list(tables = tables, plots = .save_plot_paths(plotters, outdir))
}
