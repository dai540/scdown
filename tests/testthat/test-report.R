test_that("report builder creates html file with scope and screening disclaimers", {
  obj <- scdown(scdown_example())
  path <- build_scdown_report(obj, outdir = tempfile("scdown-report-"))
  expect_true(file.exists(path))
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_match(text, "Scope and handoff")
  expect_match(text, "Markers")
  expect_match(text, "Communication")
  expect_match(text, "screening-oriented")
})

test_that("report builder supports explicit scenario selection", {
  obj <- scdown(scdown_example())
  path <- build_scdown_report(obj, scenario = "tumor_immune", outdir = tempfile("scdown-report-"))
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_match(text, "Recommended pipeline")
  expect_match(text, "tumor")
})
