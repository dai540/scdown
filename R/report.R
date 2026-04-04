#' Build a compact HTML report
#'
#' @param obj A `scdown_obj`.
#' @param genes User-selected genes.
#' @param signatures User-selected signatures.
#' @param signature_sets Optional custom signatures.
#' @param run_tests Whether to include sample-aware tests when possible.
#' @param outdir Output directory.
#'
#' @return Path to the generated HTML file.
#' @export
build_scdown_report <- function(obj,
                                genes = c("CD3D", "NKG7", "LYZ", "MS4A1"),
                                signatures = c("cytotoxic", "exhaustion", "antigen_presentation"),
                                signature_sets = NULL,
                                run_tests = TRUE,
                                outdir = tempfile("scdown-report-")) {
  stopifnot(inherits(obj, "scdown_obj"))
  .dir_create(outdir)
  results <- run_core(
    obj,
    genes = genes,
    signatures = signatures,
    signature_sets = signature_sets,
    outdir = file.path(outdir, "artifacts")
  )
  test_results <- list()
  if (isTRUE(run_tests) && obj$group_col %in% names(obj$cells)) {
    test_results$composition <- tryCatch(test_composition(obj), error = function(e) NULL)
    test_results$signatures <- tryCatch(test_signature(obj, signatures = signatures, signature_sets = signature_sets), error = function(e) NULL)
  }
  html <- c(
    "<html><head><meta charset='utf-8'><title>scdown report</title>",
    "<style>body{font-family:Arial,sans-serif;max-width:1100px;margin:40px auto;padding:0 20px;line-height:1.6} h1,h2{color:#2c3e50} img{max-width:100%;border:1px solid #ddd;margin:16px 0} pre{background:#f8f9fa;padding:16px;border-radius:8px;overflow:auto}</style>",
    "</head><body>",
    "<h1>scdown report</h1>",
    "<p>Simple downstream summary for annotated single-cell RNA-seq.</p>",
    sprintf("<p><strong>Selected genes:</strong> %s</p>", paste(genes, collapse = ", ")),
    sprintf("<p><strong>Selected signatures:</strong> %s</p>", paste(names(.resolve_gene_sets(obj, signatures, resource = "signatures", override = signature_sets)), collapse = ", ")),
    "<h2>Dataset summary</h2>",
    paste0("<pre>", paste(utils::capture.output(print(results$summary$tables$dataset_summary)), collapse = "\n"), "</pre>")
  )
  sections <- list(
    map = results$map$plots,
    markers = results$markers$plots,
    genes = results$genes$plots,
    composition = results$composition$plots,
    annotation = results$annotation$plots,
    signatures = results$signatures$plots,
    heatmap = results$heatmap$plots,
    communication = results$communication$plots
  )
  for (nm in names(sections)) {
    html <- c(html, sprintf("<h2>%s</h2>", nm))
    if (!length(sections[[nm]])) {
      html <- c(html, "<p>No plots generated.</p>")
    } else {
      html <- c(html, sprintf("<img src='%s' alt='%s'>", normalizePath(sections[[nm]], winslash = "/"), names(sections[[nm]])))
    }
  }
  if (length(test_results)) {
    html <- c(html, "<h2>Sample-aware tests</h2>")
    if (!is.null(test_results$composition)) {
      html <- c(html, "<h3>Composition</h3>", paste0("<pre>", paste(utils::capture.output(print(utils::head(test_results$composition$tables$composition_test_table, 10))), collapse = "\n"), "</pre>"))
    }
    if (!is.null(test_results$signatures)) {
      html <- c(html, "<h3>Signatures</h3>", paste0("<pre>", paste(utils::capture.output(print(utils::head(test_results$signatures$tables$signature_test_table, 10))), collapse = "\n"), "</pre>"))
    }
  }
  html <- c(html, "</body></html>")
  path <- file.path(outdir, "scdown_report.html")
  writeLines(html, path, useBytes = TRUE)
  path
}
