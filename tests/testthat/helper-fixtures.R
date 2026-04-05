example_without_group <- function() {
  demo <- scdown_example()
  demo$cells$group <- NULL
  demo
}

example_sparse_large <- function() {
  skip_if_not_installed("Matrix")
  set.seed(1)
  genes <- paste0("gene_", seq_len(300))
  cells <- paste0("cell_", seq_len(80))
  mat <- Matrix::rsparsematrix(
    nrow = length(genes),
    ncol = length(cells),
    density = 0.04
  )
  mat@x <- abs(round(mat@x * 10))
  rownames(mat) <- genes
  colnames(mat) <- cells
  meta <- data.frame(
    cell = cells,
    sample = rep(paste0("S", 1:4), each = 20),
    group = rep(c("control", "case"), each = 40),
    cell_type = rep(c("T_cell", "B_cell", "Myeloid", "NK"), each = 20),
    stringsAsFactors = FALSE
  )
  list(expr = mat, cells = meta)
}
