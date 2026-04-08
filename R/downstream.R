.marker_statistics <- function(obj, features = rownames(obj$expr), exclude_pattern = NULL) {
  labels <- as.character(obj$cells[[obj$annotation_col]])
  mat <- .analysis_matrix(obj)
  det <- .detected_matrix(obj)
  features <- intersect(features, rownames(mat))
  if (!is.null(exclude_pattern) && nzchar(exclude_pattern)) {
    features <- features[!grepl(exclude_pattern, features, ignore.case = TRUE)]
  }
  mat <- mat[features, , drop = FALSE]
  det <- det[features, , drop = FALSE]
  out <- lapply(unique(labels), function(label) {
    idx <- labels == label
    rest <- !idx
    mean_in <- .matrix_row_means(mat[, idx, drop = FALSE])
    mean_rest <- if (any(rest)) .matrix_row_means(mat[, rest, drop = FALSE]) else rep(0, length(features))
    pct_in <- .matrix_row_means(det[, idx, drop = FALSE])
    pct_rest <- if (any(rest)) .matrix_row_means(det[, rest, drop = FALSE]) else rep(0, length(features))
    data.frame(
      annotation = label,
      gene = features,
      mean_in = mean_in,
      mean_rest = mean_rest,
      pct_in = pct_in,
      pct_rest = pct_rest,
      marker_score = mean_in - mean_rest,
      marker_rank = (mean_in - mean_rest) + 0.5 * pmax(pct_in - pct_rest, 0),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out[order(out$annotation, -out$marker_rank), , drop = FALSE]
}

.average_expression_table <- function(obj, features = rownames(obj$expr)) {
  features <- intersect(features, rownames(obj$expr))
  labels <- as.character(obj$cells[[obj$annotation_col]])
  mat <- .analysis_matrix(obj)[features, , drop = FALSE]
  out <- do.call(rbind, lapply(unique(labels), function(label) {
    idx <- labels == label
    data.frame(
      annotation = label,
      gene = features,
      mean_expression = .matrix_row_means(mat[, idx, drop = FALSE]),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

.score_gene_set <- function(obj, genes, score_name) {
  mat <- .analysis_matrix(obj)
  data.frame(
    cell = colnames(mat),
    score = .rank_signature_score(mat, genes),
    score_name = score_name,
    stringsAsFactors = FALSE
  )
}

.communication_tables <- function(obj,
                                  lr_pairs = NULL,
                                  min_detected = 0.1,
                                  min_specificity = 0.15) {
  labels <- as.character(obj$cells[[obj$annotation_col]])
  pairs <- .resource_lr_pairs(obj, override = lr_pairs)
  genes <- intersect(unique(c(pairs$ligand, pairs$receptor)), rownames(obj$expr))
  empty_table <- data.frame(
    pathway = character(),
    sender = character(),
    receiver = character(),
    ligand = character(),
    receptor = character(),
    score = numeric(),
    stringsAsFactors = FALSE
  )
  if (!length(genes)) {
    return(list(
      communication_table = empty_table,
      communication_pair_summary = data.frame(sender = character(), receiver = character(), score = numeric(), stringsAsFactors = FALSE)
    ))
  }

  mat <- .analysis_matrix(obj)[genes, , drop = FALSE]
  det <- .detected_matrix(obj)[genes, , drop = FALSE]
  avg <- do.call(rbind, lapply(unique(labels), function(label) {
    idx <- labels == label
    rest <- !idx
    mean_expression <- .matrix_row_means(mat[, idx, drop = FALSE])
    mean_rest <- if (any(rest)) .matrix_row_means(mat[, rest, drop = FALSE]) else rep(0, length(genes))
    pct_detected <- .matrix_row_means(det[, idx, drop = FALSE])
    data.frame(
      annotation = label,
      gene = genes,
      mean_expression = mean_expression,
      pct_detected = pct_detected,
      specificity = pmax(mean_expression - mean_rest, 0),
      stringsAsFactors = FALSE
    )
  }))

  sender <- merge(avg, pairs, by.x = "gene", by.y = "ligand")
  receiver <- merge(avg, pairs, by.x = "gene", by.y = "receptor")
  comm <- merge(
    sender[, c("annotation", "pathway", "gene", "receptor", "pct_detected", "specificity")],
    receiver[, c("annotation", "pathway", "ligand", "gene", "pct_detected", "specificity")],
    by = "pathway",
    suffixes = c("_sender", "_receiver")
  )
  names(comm)[names(comm) == "annotation_sender"] <- "sender"
  names(comm)[names(comm) == "annotation_receiver"] <- "receiver"

  comm <- comm[
    comm$pct_detected_sender >= min_detected &
      comm$pct_detected_receiver >= min_detected &
      comm$specificity_sender >= min_specificity &
      comm$specificity_receiver >= min_specificity,
    ,
    drop = FALSE
  ]
  if (!nrow(comm)) {
    return(list(
      communication_table = empty_table,
      communication_pair_summary = data.frame(sender = character(), receiver = character(), score = numeric(), stringsAsFactors = FALSE)
    ))
  }

  comm$score <- comm$specificity_sender *
    comm$specificity_receiver *
    comm$pct_detected_sender *
    comm$pct_detected_receiver
  comm <- comm[order(-comm$score), c("pathway", "sender", "receiver", "gene_sender", "gene_receiver", "score")]
  names(comm)[4:5] <- c("ligand", "receptor")
  pair_summary <- stats::aggregate(score ~ sender + receiver, data = comm, FUN = sum)
  list(
    communication_table = comm,
    communication_pair_summary = pair_summary
  )
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
  dataset_summary <- data.frame(
    metric = c("cells", "genes", "annotations", "samples", "value_type"),
    value = c(
      ncol(obj$expr),
      nrow(obj$expr),
      length(unique(obj$cells[[obj$annotation_col]])),
      length(unique(obj$cells[[obj$sample_col]])),
      obj$value_type
    ),
    stringsAsFactors = FALSE
  )
  annotation_count <- as.data.frame(table(obj$cells[[obj$annotation_col]]), stringsAsFactors = FALSE)
  names(annotation_count) <- c("annotation", "n_cells")
  sample_count <- as.data.frame(table(obj$cells[[obj$sample_col]]), stringsAsFactors = FALSE)
  names(sample_count) <- c("sample", "n_cells")
  tables <- list(
    dataset_summary = dataset_summary,
    annotation_count = annotation_count,
    sample_count = sample_count
  )
  .write_tables(tables, outdir)
  list(tables = tables, plots = character())
}

#' Plot cluster or annotation maps
#'
#' @param obj A `scdown_obj`.
#' @param outdir Optional output directory.
#'
#' @return A list with tables and plot paths.
#' @export
plot_map <- function(obj, outdir = NULL) {
  stopifnot(inherits(obj, "scdown_obj"))
  embedding_cells <- merge(obj$embedding, obj$cells, by = "cell", all.x = TRUE, sort = FALSE)
  tables <- list(embedding_cells = embedding_cells)
  .write_tables(tables, outdir)
  plotters <- list(
    map_by_annotation = function() {
      info <- .factor_info(embedding_cells[[obj$annotation_col]])
      cols <- grDevices::hcl.colors(length(info$levels), "Set 2")
      graphics::plot(
        embedding_cells$UMAP1, embedding_cells$UMAP2,
        col = cols[info$index], pch = 19,
        xlab = "UMAP1", ylab = "UMAP2", main = "Cells by annotation"
      )
      graphics::legend("topright", legend = info$levels, col = cols, pch = 19, cex = 0.75)
    }
  )
  list(tables = tables, plots = .save_plot_paths(plotters, outdir))
}

#' Find marker genes for each annotation
#'
#' @param obj A `scdown_obj`.
#' @param top_n Number of top marker genes per annotation.
#' @param exclude_pattern Pattern for low-information genes to exclude.
#' @param min_pct_in Minimum detection rate inside the target annotation.
#' @param min_pct_diff Minimum detection improvement over the rest.
#' @param min_marker_score Minimum expression improvement over the rest.
#' @param outdir Optional output directory.
#'
#' @return A list with marker tables and plot paths.
#' @export
find_markers <- function(obj,
                         top_n = 10L,
                         exclude_pattern = "^(RPL|RPS|MT-|MALAT1$|MTRNR|HSP)",
                         min_pct_in = 0.1,
                         min_pct_diff = 0.05,
                         min_marker_score = 0.25,
                         outdir = NULL) {
  stopifnot(inherits(obj, "scdown_obj"))
  marker_table <- .marker_statistics(obj, exclude_pattern = exclude_pattern)
  keep <- marker_table$pct_in >= min_pct_in &
    (marker_table$pct_in - marker_table$pct_rest) >= min_pct_diff &
    marker_table$marker_score >= min_marker_score
  marker_hits <- marker_table[keep, , drop = FALSE]
  if (!nrow(marker_hits)) {
    marker_hits <- marker_table[marker_table$marker_score > 0, , drop = FALSE]
  }
  if (!nrow(marker_hits)) {
    marker_hits <- marker_table
  }
  top_markers <- .top_n_per_group(marker_hits, "annotation", "marker_rank", n = top_n)
  tables <- list(
    marker_table = marker_table,
    top_markers = top_markers
  )
  .write_tables(tables, outdir)
  plotters <- list(
    marker_dotplot = function() {
      top <- utils::head(top_markers, 40)
      ann_info <- .factor_info(top$annotation)
      gene_info <- .factor_info(top$gene)
      graphics::plot(
        ann_info$index, gene_info$index,
        cex = pmax(top$marker_score, 0.2), pch = 19, col = "steelblue",
        xlab = "Annotation", ylab = "Gene", xaxt = "n", yaxt = "n",
        main = "Top markers"
      )
      graphics::axis(1, at = seq_along(ann_info$levels), labels = ann_info$levels, las = 2, cex.axis = 0.7)
      graphics::axis(2, at = seq_along(gene_info$levels), labels = gene_info$levels, las = 2, cex.axis = 0.7)
    }
  )
  list(tables = tables, plots = .save_plot_paths(plotters, outdir))
}

#' Check annotation with known marker panels
#'
#' @param obj A `scdown_obj`.
#' @param panels Character vector of built-in panel names or a named list of
#'   custom panels.
#' @param marker_panels Optional named list of custom panels.
#' @param outdir Optional output directory.
#'
#' @return A list with annotation-check tables and plot paths.
#' @export
check_annotation <- function(obj, panels = NULL, marker_panels = NULL, outdir = NULL) {
  stopifnot(inherits(obj, "scdown_obj"))
  panels <- .resolve_gene_sets(obj, panels, resource = "marker_panels", override = marker_panels)
  panel_names <- names(panels)

  cell_scores <- do.call(rbind, lapply(panel_names, function(panel_name) {
    .score_gene_set(obj, panels[[panel_name]], panel_name)
  }))
  cell_scores <- merge(
    cell_scores,
    obj$cells[, c("cell", obj$annotation_col), drop = FALSE],
    by = "cell",
    all.x = TRUE
  )
  names(cell_scores)[names(cell_scores) == obj$annotation_col] <- "annotation"
  annotation_scores <- stats::aggregate(score ~ annotation + score_name, data = cell_scores, FUN = mean)
  panel_table <- data.frame(
    panel = rep(panel_names, lengths(panels)),
    gene = unlist(panels, use.names = FALSE),
    stringsAsFactors = FALSE
  )
  tables <- list(
    marker_panel_table = panel_table,
    annotation_panel_scores = annotation_scores
  )
  .write_tables(tables, outdir)
  plotters <- list(
    annotation_check_heatmap = function() {
      mat <- stats::xtabs(score ~ annotation + score_name, data = annotation_scores)
      .image_heatmap(mat, main = "Marker panel scores")
    },
    annotation_check_dotplot = function() {
      ann_info <- .factor_info(annotation_scores$annotation)
      panel_info <- .factor_info(annotation_scores$score_name)
      graphics::plot(
        ann_info$index, panel_info$index,
        cex = pmax(annotation_scores$score, 0.3),
        pch = 19, col = "firebrick",
        xlab = "Annotation", ylab = "Marker panel", xaxt = "n", yaxt = "n",
        main = "Annotation check"
      )
      graphics::axis(1, at = seq_along(ann_info$levels), labels = ann_info$levels, las = 2, cex.axis = 0.7)
      graphics::axis(2, at = seq_along(panel_info$levels), labels = panel_info$levels, las = 2, cex.axis = 0.8)
    }
  )
  list(tables = tables, plots = .save_plot_paths(plotters, outdir))
}

#' Plot average expression by annotation
#'
#' @param obj A `scdown_obj`.
#' @param features Optional character vector of genes to show.
#' @param top_n_markers Number of markers per annotation used when `features` is
#'   `NULL`.
#' @param outdir Optional output directory.
#'
#' @return A list with average-expression tables and plot paths.
#' @export
plot_average_expression <- function(obj, features = NULL, top_n_markers = 5L, outdir = NULL) {
  stopifnot(inherits(obj, "scdown_obj"))
  if (is.null(features)) {
    marker_res <- find_markers(obj, top_n = top_n_markers)
    features <- unique(marker_res$tables$top_markers$gene)
  }
  average_expression <- .average_expression_table(obj, features = features)
  average_expression_matrix <- stats::xtabs(mean_expression ~ annotation + gene, data = average_expression)
  tables <- list(
    average_expression = average_expression,
    average_expression_matrix = as.data.frame.matrix(average_expression_matrix)
  )
  .write_tables(tables, outdir)
  plotters <- list(
    average_expression_heatmap = function() {
      .image_heatmap(average_expression_matrix, main = "Average expression")
    }
  )
  list(tables = tables, plots = .save_plot_paths(plotters, outdir))
}

#' Plot similarity between annotations
#'
#' @param obj A `scdown_obj`.
#' @param top_n_features Number of variable genes used for similarity.
#' @param outdir Optional output directory.
#'
#' @return A list with similarity tables and plot paths.
#' @export
plot_cluster_similarity <- function(obj, top_n_features = 200L, outdir = NULL) {
  stopifnot(inherits(obj, "scdown_obj"))
  mat <- .analysis_matrix(obj)
  vars <- .matrix_row_vars(mat)
  vars[is.na(vars)] <- 0
  keep <- order(vars, decreasing = TRUE)[seq_len(min(length(vars), top_n_features))]
  avg <- .average_expression_table(obj, features = rownames(mat)[keep])
  sim_mat <- stats::xtabs(mean_expression ~ annotation + gene, data = avg)
  sim <- stats::cor(t(sim_mat), method = "pearson")
  similarity_long <- as.data.frame(as.table(sim), stringsAsFactors = FALSE)
  names(similarity_long) <- c("annotation_1", "annotation_2", "correlation")
  tables <- list(
    cluster_similarity = similarity_long,
    cluster_similarity_matrix = as.data.frame.matrix(sim)
  )
  .write_tables(tables, outdir)
  plotters <- list(
    cluster_similarity_heatmap = function() {
      .image_heatmap(sim, main = "Cluster similarity")
    }
  )
  list(tables = tables, plots = .save_plot_paths(plotters, outdir))
}

#' Plot signature or pathway scores
#'
#' @param obj A `scdown_obj`.
#' @param signatures Character vector of built-in signature names or a named
#'   list of signatures.
#' @param signature_sets Optional named list of custom signatures.
#' @param outdir Optional output directory.
#'
#' @return A list with signature tables and plot paths.
#' @export
plot_signature <- function(obj, signatures = NULL, signature_sets = NULL, outdir = NULL) {
  stopifnot(inherits(obj, "scdown_obj"))
  signatures <- .resolve_gene_sets(obj, signatures, resource = "signatures", override = signature_sets)
  sig_names <- names(signatures)
  cell_scores <- do.call(rbind, lapply(sig_names, function(sig_name) {
    .score_gene_set(obj, signatures[[sig_name]], sig_name)
  }))
  cell_scores <- merge(
    cell_scores,
    obj$cells[, c("cell", obj$annotation_col), drop = FALSE],
    by = "cell",
    all.x = TRUE
  )
  names(cell_scores)[names(cell_scores) == obj$annotation_col] <- "annotation"
  score_by_annotation <- stats::aggregate(score ~ annotation + score_name, data = cell_scores, FUN = mean)
  tables <- list(
    signature_score_by_cell = cell_scores,
    signature_score_by_annotation = score_by_annotation
  )
  .write_tables(tables, outdir)
  plotters <- list(
    signature_dotplot = function() {
      ann_info <- .factor_info(score_by_annotation$annotation)
      sig_info <- .factor_info(score_by_annotation$score_name)
      graphics::plot(
        ann_info$index, sig_info$index,
        cex = pmax(score_by_annotation$score, 0.3),
        pch = 19, col = "navy",
        xlab = "Annotation", ylab = "Signature", xaxt = "n", yaxt = "n",
        main = "Signature scores"
      )
      graphics::axis(1, at = seq_along(ann_info$levels), labels = ann_info$levels, las = 2, cex.axis = 0.7)
      graphics::axis(2, at = seq_along(sig_info$levels), labels = sig_info$levels, las = 2, cex.axis = 0.8)
    },
    signature_heatmap = function() {
      mat <- stats::xtabs(score ~ annotation + score_name, data = score_by_annotation)
      .image_heatmap(mat, main = "Signature heatmap")
    }
  )
  list(tables = tables, plots = .save_plot_paths(plotters, outdir))
}

#' Infer exploratory cell-cell communication
#'
#' @param obj A `scdown_obj`.
#' @param lr_pairs Optional ligand-receptor table.
#' @param min_detected Minimum detection rate for ligand and receptor.
#' @param min_specificity Minimum specificity score.
#' @param outdir Optional output directory.
#'
#' @return A list with communication tables and plot paths.
#' @export
infer_communication <- function(obj,
                                lr_pairs = NULL,
                                min_detected = 0.1,
                                min_specificity = 0.15,
                                outdir = NULL) {
  stopifnot(inherits(obj, "scdown_obj"))
  tables <- .communication_tables(
    obj = obj,
    lr_pairs = lr_pairs,
    min_detected = min_detected,
    min_specificity = min_specificity
  )
  .write_tables(tables, outdir)
  plotters <- list(
    communication_network = function() {
      df <- tables$communication_pair_summary
      if (!nrow(df)) {
        return(.plot_message("No communication pairs"))
      }
      graphics::barplot(
        df$score,
        names.arg = paste(df$sender, df$receiver, sep = " -> "),
        las = 2, col = "goldenrod",
        main = "Communication pairs", ylab = "Score"
      )
    }
  )
  list(tables = tables, plots = .save_plot_paths(plotters, outdir))
}

#' Run the minimal downstream workflow
#'
#' @param obj A `scdown_obj`.
#' @param signatures Character vector of signatures to score.
#' @param panels Character vector of marker panels to use for annotation checks.
#' @param signature_sets Optional named list of custom signatures.
#' @param marker_panels Optional named list of custom marker panels.
#' @param outdir Optional output directory.
#'
#' @return A named list of downstream results.
#' @export
run_core <- function(obj,
                     signatures = c("cytotoxic", "antigen_presentation"),
                     panels = c("t_cell", "nk", "b_cell", "monocyte", "macrophage", "dendritic"),
                     signature_sets = NULL,
                     marker_panels = NULL,
                     outdir = NULL) {
  stopifnot(inherits(obj, "scdown_obj"))
  if (!is.null(outdir)) {
    .dir_create(outdir)
  }
  list(
    summary = summarize_dataset(obj, outdir = .subdir_or_null(outdir, "summary")),
    map = plot_map(obj, outdir = .subdir_or_null(outdir, "map")),
    markers = find_markers(obj, outdir = .subdir_or_null(outdir, "markers")),
    annotation_check = check_annotation(
      obj,
      panels = panels,
      marker_panels = marker_panels,
      outdir = .subdir_or_null(outdir, "annotation_check")
    ),
    average_expression = plot_average_expression(obj, outdir = .subdir_or_null(outdir, "average_expression")),
    cluster_similarity = plot_cluster_similarity(obj, outdir = .subdir_or_null(outdir, "cluster_similarity")),
    signatures = plot_signature(
      obj,
      signatures = signatures,
      signature_sets = signature_sets,
      outdir = .subdir_or_null(outdir, "signatures")
    ),
    communication = infer_communication(obj, outdir = .subdir_or_null(outdir, "communication"))
  )
}
