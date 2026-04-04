#' Scope table for `scdown`
#'
#' `scdown` is strongest after annotation, when the goal is fast review,
#' lightweight testing, and compact reporting. This helper makes the intended
#' boundary explicit so users can see where `scdown` is enough and where a
#' specialist method should take over.
#'
#' @param scenario One of `"generic"`, `"pbmc"`, `"tumor_immune"`, or
#'   `"patient_comparison"`.
#'
#' @return A data frame with workflow stages, whether `scdown` is a good fit,
#'   and which specialist tools are better suited when the analysis needs to go
#'   further.
#' @export
scope_table <- function(scenario = c("generic", "pbmc", "tumor_immune", "patient_comparison")) {
  scenario <- match.arg(scenario)

  base <- data.frame(
    stage = c(
      "qc_and_cleanup",
      "normalize_cluster_integrate",
      "label_assignment",
      "exploratory_review",
      "lightweight_group_tests",
      "publication_grade_state_models",
      "publication_grade_abundance",
      "publication_grade_communication"
    ),
    use_scdown = c("no", "no", "no", "yes", "yes", "handoff", "handoff", "handoff"),
    scdown_role = c(
      "not designed for this stage",
      "not designed for this stage",
      "annotation sanity check only",
      "core use case",
      "good for quick sample-aware screening",
      "handoff point",
      "handoff point",
      "handoff point"
    ),
    specialist_tools = c(
      "SoupX, scDblFinder, scater, scuttle",
      "Seurat, scran, batchelor, scater",
      "SingleR, Azimuth, manual curation",
      "scdown",
      "scdown",
      "muscat, dreamlet",
      "miloR",
      "CellChat"
    ),
    note = c(
      "Handle droplets, ambient RNA, doublets, and low-quality cells upstream.",
      "Settle normalization, embedding, clustering, and integration before `scdown`.",
      "Use `scdown` to review labels, not to create them from scratch.",
      "Maps, markers, genes, composition, signatures, communication, and reports are the main job of `scdown`.",
      "Useful for prioritizing endpoints, but not a substitute for formal specialist models.",
      "Use specialist models when differential state becomes a primary claim.",
      "Use neighborhood- or model-based abundance tools when abundance is central.",
      "Use a dedicated communication framework when communication is a paper-level result."
    ),
    stringsAsFactors = FALSE
  )

  scenario_note <- switch(
    scenario,
    generic = c(
      "Generic single-cell workflow boundary.",
      "Generic single-cell workflow boundary.",
      "Generic single-cell workflow boundary.",
      "Review any annotated object quickly.",
      "Good for initial prioritization across groups.",
      "Promote key marker or state claims into specialist DE/DS models.",
      "Promote abundance claims into specialist DA methods.",
      "Promote communication claims into a specialist communication framework."
    ),
    pbmc = c(
      "Standard PBMC preprocessing stays upstream.",
      "PBMC clustering and embedding are typically handled upstream.",
      "PBMC labels are usually stabilized with reference mapping plus manual review.",
      "PBMC is a strong fit for `scdown` review.",
      "Often enough for internal review or targeted follow-up.",
      "Escalate only when PBMC state differences become central.",
      "Escalate when abundance changes are part of the main claim.",
      "Escalate when communication is more than exploratory."
    ),
    tumor_immune = c(
      "Tumor immune preprocessing and integration stay upstream.",
      "Tumor immune batch handling belongs upstream.",
      "Tumor immune labeling needs specialist annotation plus manual curation.",
      "A strong fit for marker, signature, composition, and exploratory communication review.",
      "Good for deciding which tumor-immune endpoints deserve formal modeling.",
      "Escalate tumor state claims into donor-aware specialist models.",
      "Escalate abundance shifts into neighborhood or compositional specialist methods.",
      "Escalate communication claims into a richer ligand-receptor framework."
    ),
    patient_comparison = c(
      "Patient-aware QC must be solved upstream.",
      "Patient-aware normalization, integration, and embedding stay upstream.",
      "Patient labels should be stabilized before comparison.",
      "Use `scdown` to review and prioritize patient comparison endpoints.",
      "Good for a first pass, not for the final claim set.",
      "Use specialist mixed or pseudobulk models for patient comparisons.",
      "Use specialist DA models for patient abundance claims.",
      "Use specialist communication models when patient-group signaling is central."
    )
  )
  base$scenario_note <- scenario_note
  base
}

.sample_metadata_table <- function(obj) {
  sample_cols <- unique(c(obj$sample_col, obj$group_col, obj$batch_col))
  sample_cols <- sample_cols[!is.na(sample_cols) & nzchar(sample_cols) & sample_cols %in% names(obj$cells)]
  if (!length(sample_cols)) {
    return(data.frame())
  }
  sample_meta <- unique(obj$cells[, sample_cols, drop = FALSE])
  rownames(sample_meta) <- NULL
  sample_meta
}

.as_sce_for_handoff <- function(obj, assay = c("counts", "analysis")) {
  assay <- match.arg(assay)
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE) ||
      !requireNamespace("SummarizedExperiment", quietly = TRUE) ||
      !requireNamespace("S4Vectors", quietly = TRUE)) {
    return(NULL)
  }
  assays <- list()
  assays$counts <- obj$expr
  assays$logcounts <- .analysis_matrix(obj)
  if (identical(assay, "analysis")) {
    assays <- assays[c("logcounts", "counts")]
  }
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = assays,
    colData = S4Vectors::DataFrame(obj$cells)
  )
  rownames(sce) <- rownames(obj$expr)
  colnames(sce) <- colnames(obj$expr)
  if (nrow(obj$embedding)) {
    emb <- as.matrix(obj$embedding[, c("UMAP1", "UMAP2"), drop = FALSE])
    rownames(emb) <- obj$embedding$cell
    SingleCellExperiment::reducedDim(sce, "UMAP") <- emb
  }
  sce
}

#' Export `scdown` data for specialist handoff
#'
#' This helper does not run `muscat`, `dreamlet`, `miloR`, or CellChat for you.
#' Instead, it prepares a compact handoff object with the matrices, metadata,
#' and embeddings that those downstream specialist tools usually need.
#'
#' @param obj A `scdown_obj`.
#' @param target One of `"muscat"`, `"dreamlet"`, `"miloR"`, or `"CellChat"`.
#' @param assay Which matrix to emphasize in the handoff: `"counts"` or
#'   `"analysis"`.
#' @param include_lineage Whether to ensure a lineage column is included in the
#'   exported cell metadata.
#'
#' @return A named list with matrices, metadata, guidance, and, when the needed
#'   packages are installed, a `SingleCellExperiment` object for Bioconductor
#'   handoff targets.
#' @export
export_for_handoff <- function(obj,
                               target = c("muscat", "dreamlet", "miloR", "CellChat"),
                               assay = c("counts", "analysis"),
                               include_lineage = TRUE) {
  stopifnot(inherits(obj, "scdown_obj"))
  target <- match.arg(target)
  assay <- match.arg(assay)
  if (isTRUE(include_lineage) && !obj$lineage_col %in% names(obj$cells)) {
    obj <- collapse_lineage(obj)
  }

  assay_matrix <- if (identical(assay, "counts")) obj$expr else .analysis_matrix(obj)
  sample_meta <- .sample_metadata_table(obj)
  scope <- scope_table(
    if (obj$group_col %in% names(obj$cells)) "patient_comparison" else "generic"
  )

  handoff <- list(
    target = target,
    assay = assay,
    expr = assay_matrix,
    cells = obj$cells,
    sample_metadata = sample_meta,
    embedding = obj$embedding,
    cluster_col = obj$annotation_col,
    sample_col = obj$sample_col,
    group_col = if (obj$group_col %in% names(obj$cells)) obj$group_col else NULL,
    batch_col = if (!is.null(obj$batch_col) && obj$batch_col %in% names(obj$cells)) obj$batch_col else NULL,
    lineage_col = if (obj$lineage_col %in% names(obj$cells)) obj$lineage_col else NULL,
    recommended_scope = scope,
    guidance = switch(
      target,
      muscat = "Use this export as the starting point for pseudobulk or cell-level DS analyses in muscat.",
      dreamlet = "Use this export as the starting point for donor-aware pseudobulk mixed models in dreamlet.",
      miloR = "Use this export as the starting point for neighborhood-based abundance testing in miloR.",
      CellChat = "Use this export as the starting point for specialist communication modeling in CellChat."
    )
  )

  if (target %in% c("muscat", "dreamlet", "miloR")) {
    handoff$sce <- .as_sce_for_handoff(obj, assay = assay)
  }
  if (identical(target, "CellChat")) {
    handoff$cellchat_input <- list(
      data = assay_matrix,
      meta = obj$cells,
      group.by = obj$annotation_col
    )
  }

  class(handoff) <- c(paste0("scdown_", tolower(target), "_export"), "scdown_handoff_export")
  handoff
}
