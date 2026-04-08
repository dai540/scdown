test_that("report writes html and artifact directories", {
  obj <- scdown(scdown_example(), annotation_col = "annotation", sample_col = "sample")
  outdir <- tempfile("scdown-report-")
  path <- build_scdown_report(obj, outdir = outdir)

  expect_true(file.exists(path))
  html <- readLines(path, warn = FALSE)
  expect_true(any(grepl("Cluster map", html, fixed = TRUE)))
  expect_true(any(grepl("Marker genes", html, fixed = TRUE)))
  expect_true(dir.exists(file.path(outdir, "artifacts", "markers")))
  expect_true(dir.exists(file.path(outdir, "artifacts", "communication")))
})
