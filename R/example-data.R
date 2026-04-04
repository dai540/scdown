#' Example annotated scRNA-seq dataset
#'
#' @return A list with `expr`, `cells`, and `embedding`.
#' @export
scdown_example <- function() {
  design <- expand.grid(
    sample = c("S1", "S2", "S3", "S4"),
    cell_type = c("CD4_T", "CD8_T", "NK", "B_cell", "Mono", "Macro"),
    replicate = seq_len(2),
    stringsAsFactors = FALSE
  )
  design$group <- ifelse(design$sample %in% c("S1", "S2"), "control", "case")
  design$cell <- paste0("cell_", seq_len(nrow(design)))
  cells <- design[, c("cell", "sample", "group", "cell_type")]
  centers <- data.frame(
    cell_type = c("CD4_T", "CD8_T", "NK", "B_cell", "Mono", "Macro"),
    UMAP1 = c(-2.2, -1.5, -0.3, 0.6, 1.7, 2.3),
    UMAP2 = c(2.0, 1.4, 1.7, -0.6, -1.1, -0.3),
    stringsAsFactors = FALSE
  )
  embedding <- merge(cells, centers, by = "cell_type", all.x = TRUE, sort = FALSE)
  sample_offset <- c(S1 = -0.15, S2 = 0.05, S3 = 0.12, S4 = 0.22)
  embedding$UMAP1 <- embedding$UMAP1 + unname(sample_offset[embedding$sample]) + c(-0.08, 0.08)[((seq_len(nrow(embedding)) - 1L) %% 2L) + 1L]
  embedding$UMAP2 <- embedding$UMAP2 + c(0.1, -0.1)[((seq_len(nrow(embedding)) - 1L) %% 2L) + 1L]
  embedding <- embedding[, c("cell", "UMAP1", "UMAP2")]
  genes <- c("CD3D", "IL7R", "NKG7", "GNLY", "MS4A1", "CD79A", "LYZ", "FCN1", "C1QC", "APOE", "HLA-DRA", "CCL5", "CCR5", "CXCL8", "CXCR2")
  marker_map <- list(
    CD4_T = c("CD3D", "IL7R"),
    CD8_T = c("CD3D", "CCL5"),
    NK = c("NKG7", "GNLY"),
    B_cell = c("MS4A1", "CD79A"),
    Mono = c("LYZ", "FCN1", "CXCL8"),
    Macro = c("C1QC", "APOE", "HLA-DRA")
  )
  expr <- matrix(1L, nrow = length(genes), ncol = nrow(cells), dimnames = list(genes, cells$cell))
  for (cell_type in names(marker_map)) {
    cell_ids <- cells$cell[cells$cell_type == cell_type]
    genes_hit <- marker_map[[cell_type]]
    expr[genes_hit, cell_ids] <- expr[genes_hit, cell_ids] + 12L
  }
  expr["CCR5", cells$cell[cells$cell_type %in% c("CD4_T", "CD8_T", "NK")]] <- expr["CCR5", cells$cell[cells$cell_type %in% c("CD4_T", "CD8_T", "NK")]] + 8L
  expr["CXCR2", cells$cell[cells$cell_type == "Mono"]] <- expr["CXCR2", cells$cell[cells$cell_type == "Mono"]] + 7L
  expr["CCL5", cells$cell[cells$cell_type %in% c("CD8_T", "NK") & cells$group == "case"]] <- expr["CCL5", cells$cell[cells$cell_type %in% c("CD8_T", "NK") & cells$group == "case"]] + 4L
  expr["HLA-DRA", cells$cell[cells$cell_type %in% c("Mono", "Macro") & cells$group == "case"]] <- expr["HLA-DRA", cells$cell[cells$cell_type %in% c("Mono", "Macro") & cells$group == "case"]] + 3L
  expr["MS4A1", cells$cell[cells$cell_type == "B_cell" & cells$group == "control"]] <- expr["MS4A1", cells$cell[cells$cell_type == "B_cell" & cells$group == "control"]] + 2L
  expr <- expr + matrix(rep(c(0L, 1L, 2L, 0L), length.out = length(expr)), nrow = nrow(expr))
  if (requireNamespace("Matrix", quietly = TRUE)) {
    expr <- Matrix::Matrix(expr, sparse = TRUE)
  }
  list(expr = expr, cells = cells, embedding = embedding)
}
