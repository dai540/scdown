#' Built-in immune signatures
#'
#' @return A named list of character vectors.
#' @export
immune_signatures <- function() {
  list(
    cytotoxic = c("NKG7", "GNLY", "PRF1", "GZMB"),
    exhaustion = c("PDCD1", "LAG3", "HAVCR2", "CTLA4"),
    b_cell = c("MS4A1", "CD79A", "CD79B"),
    inflammatory_monocyte = c("LYZ", "FCN1", "S100A8", "S100A9"),
    antigen_presentation = c("HLA-DRA", "HLA-DPA1", "CD74")
  )
}

#' Built-in immune marker panels
#'
#' @return A named list of character vectors.
#' @export
immune_marker_panels <- function() {
  list(
    t_cell = c("CD3D", "TRAC", "IL7R"),
    nk = c("NKG7", "GNLY", "PRF1"),
    b_cell = c("MS4A1", "CD79A", "CD79B"),
    monocyte = c("LYZ", "FCN1", "S100A8"),
    macrophage = c("C1QC", "APOE", "CTSB"),
    dendritic = c("FCER1A", "CLEC10A", "CD1C")
  )
}

.default_lr_pairs <- function() {
  data.frame(
    ligand = c("CCL5", "CCL19", "CCL2", "CXCL8", "CXCL9", "CXCL10", "IFNG", "TNF", "MIF"),
    receptor = c("CCR5", "CCR7", "CCR2", "CXCR2", "CXCR3", "CXCR3", "IFNGR1", "TNFRSF1A", "CD74"),
    pathway = c("chemokine", "chemokine", "chemokine", "chemokine", "chemokine", "chemokine", "ifn", "tnf", "mif"),
    stringsAsFactors = FALSE
  )
}

.validate_named_gene_sets <- function(x, what) {
  if (is.null(x)) {
    return(NULL)
  }
  if (!is.list(x) || is.null(names(x)) || any(!nzchar(names(x)))) {
    stop(sprintf("`%s` must be a named list of character vectors.", what), call. = FALSE)
  }
  lapply(x, function(genes) unique(as.character(genes[!is.na(genes) & nzchar(genes)])))
}

.validate_lr_pairs <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  .assert_columns(x, c("ligand", "receptor", "pathway"), "lr_pairs")
  x[, c("ligand", "receptor", "pathway"), drop = FALSE]
}

.merge_named_lists <- function(base, override = NULL) {
  if (is.null(override)) {
    return(base)
  }
  base[names(override)] <- override
  base
}

.resource_signatures <- function(obj = NULL, override = NULL) {
  base <- immune_signatures()
  if (!is.null(obj) && !is.null(obj$resources$signatures)) {
    base <- .merge_named_lists(base, obj$resources$signatures)
  }
  .merge_named_lists(base, .validate_named_gene_sets(override, "signatures"))
}

.resource_marker_panels <- function(obj = NULL, override = NULL) {
  base <- immune_marker_panels()
  if (!is.null(obj) && !is.null(obj$resources$marker_panels)) {
    base <- .merge_named_lists(base, obj$resources$marker_panels)
  }
  .merge_named_lists(base, .validate_named_gene_sets(override, "marker_panels"))
}

.resource_lr_pairs <- function(obj = NULL, override = NULL) {
  override <- .validate_lr_pairs(override)
  if (!is.null(override)) {
    return(override)
  }
  if (!is.null(obj) && !is.null(obj$resources$lr_pairs)) {
    return(obj$resources$lr_pairs)
  }
  .default_lr_pairs()
}

.resolve_gene_sets <- function(obj, gene_sets, resource, override = NULL) {
  lookup <- switch(
    resource,
    signatures = .resource_signatures(obj, override = override),
    marker_panels = .resource_marker_panels(obj, override = override)
  )
  if (missing(gene_sets) || is.null(gene_sets)) {
    return(lookup)
  }
  if (is.list(gene_sets)) {
    return(.validate_named_gene_sets(gene_sets, resource))
  }
  gene_sets <- unique(as.character(gene_sets))
  missing_names <- setdiff(gene_sets, names(lookup))
  if (length(missing_names)) {
    stop(sprintf("Unknown %s: %s.", resource, paste(missing_names, collapse = ", ")), call. = FALSE)
  }
  lookup[gene_sets]
}

#' List available signatures
#'
#' @param obj Optional `scdown_obj` with custom signatures.
#'
#' @return A character vector.
#' @export
list_signatures <- function(obj = NULL) {
  names(.resource_signatures(obj))
}

#' Get one signature
#'
#' @param name Signature name.
#' @param obj Optional `scdown_obj` with custom signatures.
#'
#' @return Character vector of genes.
#' @export
get_signature <- function(name, obj = NULL) {
  sigs <- .resource_signatures(obj)
  name <- as.character(name)[1]
  if (!name %in% names(sigs)) {
    stop("Unknown signature.", call. = FALSE)
  }
  sigs[[name]]
}
