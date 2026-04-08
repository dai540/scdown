#' Example annotated scRNA-seq dataset
#'
#' @return A list with `expr`, `cells`, and `embedding`.
#' @export
scdown_example <- function() {
  design <- expand.grid(
    sample = c("S1", "S2", "S3"),
    annotation = c("CD4_T", "CD8_T", "NK", "B_cell", "Monocyte", "Macrophage"),
    replicate = seq_len(6),
    stringsAsFactors = FALSE
  )
  design$cell <- paste0("cell_", seq_len(nrow(design)))
  cells <- design[, c("cell", "sample", "annotation")]

  centers <- data.frame(
    annotation = c("CD4_T", "CD8_T", "NK", "B_cell", "Monocyte", "Macrophage"),
    UMAP1 = c(-2.4, -1.6, -0.5, 0.7, 1.8, 2.6),
    UMAP2 = c(2.0, 1.2, 1.7, -0.4, -1.2, -0.3),
    stringsAsFactors = FALSE
  )
  embedding <- merge(cells, centers, by = "annotation", all.x = TRUE, sort = FALSE)
  sample_shift <- c(S1 = -0.15, S2 = 0.00, S3 = 0.15)
  embedding$UMAP1 <- embedding$UMAP1 + unname(sample_shift[embedding$sample]) +
    rep(c(-0.06, 0.00, 0.06), length.out = nrow(embedding))
  embedding$UMAP2 <- embedding$UMAP2 + rep(c(0.08, -0.08), length.out = nrow(embedding))
  embedding <- embedding[, c("cell", "UMAP1", "UMAP2")]

  genes <- c(
    "CD3D", "IL7R", "LTB",
    "CCL5", "NKG7", "GNLY", "PRF1",
    "MS4A1", "CD79A", "CD74",
    "LYZ", "FCN1", "S100A8", "S100A9",
    "C1QC", "APOE", "HLA-DRA", "FCER1A",
    "CCR5", "CXCL8", "CXCR2"
  )
  marker_map <- list(
    CD4_T = c("CD3D", "IL7R", "LTB"),
    CD8_T = c("CD3D", "CCL5", "PRF1"),
    NK = c("NKG7", "GNLY", "CCL5"),
    B_cell = c("MS4A1", "CD79A", "CD74"),
    Monocyte = c("LYZ", "FCN1", "S100A8", "S100A9", "CXCL8", "CXCR2"),
    Macrophage = c("C1QC", "APOE", "HLA-DRA", "FCER1A")
  )

  expr <- matrix(0, nrow = length(genes), ncol = nrow(cells), dimnames = list(genes, cells$cell))
  for (annotation in names(marker_map)) {
    cell_ids <- cells$cell[cells$annotation == annotation]
    expr[marker_map[[annotation]], cell_ids] <- expr[marker_map[[annotation]], cell_ids] + 12
  }
  expr["CCR5", cells$cell[cells$annotation %in% c("CD4_T", "CD8_T", "NK")]] <-
    expr["CCR5", cells$cell[cells$annotation %in% c("CD4_T", "CD8_T", "NK")]] + 6
  expr["HLA-DRA", cells$cell[cells$annotation %in% c("Monocyte", "Macrophage", "B_cell")]] <-
    expr["HLA-DRA", cells$cell[cells$annotation %in% c("Monocyte", "Macrophage", "B_cell")]] + 3
  expr <- expr + matrix(rep(c(0, 1, 2), length.out = length(expr)), nrow = nrow(expr))
  if (requireNamespace("Matrix", quietly = TRUE)) {
    expr <- Matrix::Matrix(expr, sparse = TRUE)
  }

  list(expr = expr, cells = cells, embedding = embedding)
}
