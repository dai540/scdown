#' Build a compact HTML report
#'
#' @param obj A `scdown_obj`.
#' @param signatures Character vector of built-in signature names or a named
#'   list of signatures.
#' @param panels Character vector of built-in marker panel names or a named list
#'   of marker panels.
#' @param signature_sets Optional named list of custom signatures.
#' @param marker_panels Optional named list of custom marker panels.
#' @param outdir Output directory.
#'
#' @return Path to the generated HTML report.
#' @export
build_scdown_report <- function(obj,
                                signatures = c("cytotoxic", "antigen_presentation"),
                                panels = c("t_cell", "nk", "b_cell", "monocyte", "macrophage", "dendritic"),
                                signature_sets = NULL,
                                marker_panels = NULL,
                                outdir = tempfile("scdown-report-")) {
  stopifnot(inherits(obj, "scdown_obj"))
  .dir_create(outdir)
  results <- run_core(
    obj = obj,
    signatures = signatures,
    panels = panels,
    signature_sets = signature_sets,
    marker_panels = marker_panels,
    outdir = file.path(outdir, "artifacts")
  )

  html <- c(
    "<html><head><meta charset='utf-8'><title>scdown report</title>",
    "<style>body{font-family:Arial,sans-serif;max-width:1100px;margin:40px auto;padding:0 20px;line-height:1.6} h1,h2,h3{color:#2c3e50} img{max-width:100%;border:1px solid #ddd;margin:16px 0} pre{background:#f8f9fa;padding:16px;border-radius:8px;overflow:auto} table{border-collapse:collapse} td,th{border:1px solid #ddd;padding:6px 8px}</style>",
    "</head><body>",
    "<h1>scdown report</h1>",
    "<p>Minimal downstream report for one annotated single-cell dataset.</p>",
    "<p><strong>Included analyses:</strong> map, markers, marker-panel annotation check, average expression, cluster similarity, signature scoring, and exploratory communication.</p>",
    "<p><strong>Scope:</strong> this report is for within-dataset structure and interpretation. It is not a substitute for group-aware statistical modeling.</p>",
    "<h2>Dataset summary</h2>",
    paste0("<pre>", paste(utils::capture.output(print(results$summary$tables$dataset_summary)), collapse = "\n"), "</pre>")
  )

  section_order <- c(
    map = "Cluster map",
    markers = "Marker genes",
    annotation_check = "Annotation check",
    average_expression = "Average expression",
    cluster_similarity = "Cluster similarity",
    signatures = "Signature scores",
    communication = "Cell-cell communication"
  )

  for (nm in names(section_order)) {
    res <- results[[nm]]
    html <- c(html, sprintf("<h2>%s</h2>", section_order[[nm]]))
    if (length(res$tables)) {
      first_table <- res$tables[[1]]
      html <- c(
        html,
        "<h3>Table preview</h3>",
        paste0("<pre>", paste(utils::capture.output(print(utils::head(first_table, 12))), collapse = "\n"), "</pre>")
      )
    }
    if (!length(res$plots)) {
      html <- c(html, "<p>No plots generated.</p>")
    } else {
      html <- c(html, sprintf("<img src='%s' alt='%s'>", file.path("artifacts", nm, basename(res$plots)), names(res$plots)))
    }
  }

  html <- c(html, "</body></html>")
  path <- file.path(outdir, "scdown_report.html")
  writeLines(html, path, useBytes = TRUE)
  path
}
