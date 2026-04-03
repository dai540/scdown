#' Built-in immune signatures
#'
#' @return A named list of character vectors.
#' @export
immune_signatures <- function() {
  list(
    naive_t = c("CCR7", "IL7R", "LTB"),
    cytotoxic = c("NKG7", "GNLY", "PRF1", "GZMB"),
    exhaustion = c("PDCD1", "LAG3", "HAVCR2", "CTLA4"),
    b_cell = c("MS4A1", "CD79A", "CD79B"),
    plasma = c("JCHAIN", "MZB1", "SDC1"),
    inflammatory_monocyte = c("LYZ", "FCN1", "S100A8", "S100A9"),
    macrophage = c("C1QC", "APOE", "CTSB", "HLA-DRA"),
    antigen_presentation = c("HLA-DRA", "HLA-DPA1", "CD74")
  )
}

#' Built-in immune marker panels
#'
#' @return A named list of character vectors.
#' @export
immune_marker_panels <- function() {
  list(
    t_cell = c("CD3D", "TRAC", "LTB"),
    nk = c("NKG7", "GNLY", "PRF1"),
    b_cell = c("MS4A1", "CD79A", "CD79B"),
    plasma = c("JCHAIN", "MZB1", "SDC1"),
    monocyte = c("LYZ", "FCN1", "S100A8"),
    macrophage = c("C1QC", "APOE", "CTSB"),
    dendritic = c("FCER1A", "CLEC10A", "CD1C")
  )
}

#' Default lineage mapping
#'
#' @return A data frame with pattern and lineage.
#' @export
default_lineage_mapping <- function() {
  data.frame(
    pattern = c("treg", "cd4", "cd8", "t_cell", "t cell", "nk", "b_cell", "b cell", "plasma", "mono", "monocyte", "macro", "macrophage", "dc", "dendritic"),
    lineage = c("T/NK", "T/NK", "T/NK", "T/NK", "T/NK", "T/NK", "B/Plasma", "B/Plasma", "B/Plasma", "Myeloid", "Myeloid", "Myeloid", "Myeloid", "Myeloid", "Myeloid"),
    stringsAsFactors = FALSE
  )
}

#' List available signatures
#'
#' @return A character vector.
#' @export
list_signatures <- function() {
  names(immune_signatures())
}

#' Get one built-in signature
#'
#' @param name Signature name.
#'
#' @return Character vector of genes.
#' @export
get_signature <- function(name) {
  sigs <- immune_signatures()
  name <- as.character(name)[1]
  if (!name %in% names(sigs)) {
    stop("Unknown signature.", call. = FALSE)
  }
  sigs[[name]]
}

.lr_pairs <- function() {
  data.frame(
    ligand = c("CCL5", "CCL4", "CCL19", "CCL2", "CXCL8", "CXCL9", "CXCL10", "CXCL12", "IFNG", "TNF", "TGFB1", "IL7", "XCL1", "MIF"),
    receptor = c("CCR5", "CCR5", "CCR7", "CCR2", "CXCR2", "CXCR3", "CXCR3", "CXCR4", "IFNGR1", "TNFRSF1A", "TGFBR1", "IL7R", "XCR1", "CD74"),
    pathway = c("chemokine", "chemokine", "chemokine", "chemokine", "chemokine", "chemokine", "chemokine", "chemokine", "ifn", "tnf", "tgfb", "il7", "xcl", "mif"),
    stringsAsFactors = FALSE
  )
}
