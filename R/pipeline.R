#' Recommended Pipeline Around `scdown`
#'
#' @param scenario One of `"generic"`, `"pbmc"`, `"tumor_immune"`, or
#'   `"patient_comparison"`.
#'
#' @return A data frame describing suggested stages, goals, tools, and where
#'   `scdown` fits.
#' @export
recommended_pipeline <- function(scenario = c("generic", "pbmc", "tumor_immune", "patient_comparison")) {
  scenario <- match.arg(scenario)

  generic <- data.frame(
    stage = c(
      "raw_qc",
      "normalize_cluster",
      "annotate",
      "review_with_scdown",
      "rigorous_testing"
    ),
    goal = c(
      "remove technical artifacts and low-quality droplets",
      "normalize, reduce dimensions, cluster, and integrate batches",
      "assign biologically meaningful cell labels",
      "inspect maps, markers, genes, composition, signatures, and communication",
      "run publication-grade differential state, abundance, or communication analyses"
    ),
    recommended_tools = c(
      "SoupX, scDblFinder, scater, scuttle",
      "Seurat, scater, scran, batchelor",
      "SingleR, Azimuth, manual curation",
      "scdown",
      "muscat, dreamlet, miloR, CellChat"
    ),
    scdown_role = c("no", "no", "no", "core", "handoff"),
    stringsAsFactors = FALSE
  )

  pbmc <- data.frame(
    stage = c(
      "preprocess",
      "annotate",
      "review_with_scdown",
      "targeted_testing"
    ),
    goal = c(
      "standard PBMC QC, normalization, and clustering",
      "assign canonical immune populations",
      "rapidly review markers, cytotoxic programs, composition, and report",
      "test a few predefined signatures or compositions across groups"
    ),
    recommended_tools = c(
      "Seurat or scater/scran",
      "SingleR, Azimuth, manual marker review",
      "scdown",
      "scdown, then muscat or dreamlet if the question becomes formal"
    ),
    scdown_role = c("no", "no", "core", "optional_test"),
    stringsAsFactors = FALSE
  )

  tumor_immune <- data.frame(
    stage = c(
      "preprocess",
      "annotate",
      "review_with_scdown",
      "rigorous_followup"
    ),
    goal = c(
      "clean tumor-infiltrating immune profiles and integrate samples",
      "label T/NK, B/plasma, myeloid, and dendritic compartments",
      "review immune state markers, composition, signatures, and exploratory communication",
      "promote key abundance, state, or communication claims into specialist methods"
    ),
    recommended_tools = c(
      "Seurat, scater/scran, batchelor",
      "SingleR plus manual tumor-immune curation",
      "scdown",
      "muscat, dreamlet, miloR, CellChat"
    ),
    scdown_role = c("no", "no", "core", "handoff"),
    stringsAsFactors = FALSE
  )

  patient_comparison <- data.frame(
    stage = c(
      "preprocess",
      "annotate",
      "review_with_scdown",
      "formal_models"
    ),
    goal = c(
      "prepare a donor-aware object with clean sample metadata",
      "stabilize labels before comparison",
      "use `scdown` to decide which endpoints are worth formal modeling",
      "fit sample-aware models for state, abundance, and pathway differences"
    ),
    recommended_tools = c(
      "Seurat or Bioconductor upstream stack",
      "SingleR, Azimuth, manual review",
      "scdown",
      "dreamlet, muscat, miloR"
    ),
    scdown_role = c("no", "no", "core", "handoff"),
    stringsAsFactors = FALSE
  )

  switch(
    scenario,
    generic = generic,
    pbmc = pbmc,
    tumor_immune = tumor_immune,
    patient_comparison = patient_comparison
  )
}
