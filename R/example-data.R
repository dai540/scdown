#' Example annotated scRNA-seq dataset
#'
#' @return A list with `expr`, `cells`, and `embedding`.
#' @export
scdown_example <- function() {
  cells <- data.frame(
    cell = paste0("cell_", seq_len(24)),
    sample = rep(c("S1", "S2", "S3"), each = 8),
    group = rep(c("control", "case"), each = 12),
    cell_type = rep(c("CD4_T", "CD8_T", "NK", "B_cell", "Mono", "Macro"), each = 4),
    stringsAsFactors = FALSE
  )
  embedding <- data.frame(
    cell = cells$cell,
    UMAP1 = rep(c(-2.2, -1.8, -0.1, 0.2, 1.8, 2.2), each = 4),
    UMAP2 = c(2.2, 1.8, 2.4, 1.6, 1.9, 1.5, 2.3, 1.7, -1.1, -0.8, -1.3, -0.7, -1.0, -0.6, -1.2, -0.5, 1.0, 1.3, 0.8, 1.5, 1.2, 1.4, 0.9, 1.6),
    stringsAsFactors = FALSE
  )
  genes <- c("CD3D", "IL7R", "NKG7", "GNLY", "MS4A1", "CD79A", "LYZ", "FCN1", "C1QC", "APOE", "HLA-DRA", "CCL5", "CCR5", "CXCL8", "CXCR2")
  marker_map <- list(
    CD4_T = c("CD3D", "IL7R"),
    CD8_T = c("CD3D", "CCL5"),
    NK = c("NKG7", "GNLY"),
    B_cell = c("MS4A1", "CD79A"),
    Mono = c("LYZ", "FCN1", "CXCL8"),
    Macro = c("C1QC", "APOE", "HLA-DRA")
  )
  expr <- expand.grid(cell = cells$cell, gene = genes, stringsAsFactors = FALSE)
  expr <- merge(expr, cells, by = "cell", all.x = TRUE, sort = FALSE)
  expr$value <- 0.2
  for (ct in names(marker_map)) {
    hit <- expr$cell_type == ct & expr$gene %in% marker_map[[ct]]
    expr$value[hit] <- expr$value[hit] + 3
  }
  expr$value[expr$gene == "CCR5" & expr$cell_type %in% c("CD4_T", "CD8_T", "NK")] <- 2.4
  expr$value[expr$gene == "CXCR2" & expr$cell_type == "Mono"] <- 2.4
  expr$value <- expr$value + ((seq_len(nrow(expr)) - 1L) %% 5L) / 20
  expr <- expr[, c("sample", "cell", "cell_type", "group", "gene", "value")]
  list(expr = expr, cells = cells, embedding = embedding)
}
