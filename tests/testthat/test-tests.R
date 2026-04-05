test_that("sample-aware test functions return expected tables", {
  obj <- collapse_lineage(scdown(scdown_example()))
  expect_true("top_marker_tests" %in% names(test_markers(obj, top_n = 3)$tables))
  expect_true("composition_test_table" %in% names(test_composition(obj)$tables))
  expect_true("signature_test_table" %in% names(test_signature(obj, signatures = c("cytotoxic", "b_cell"))$tables))
  expect_true("communication_test_table" %in% names(test_communication(obj, n_perm = 5)$tables))
})

test_that("group-dependent tests fail clearly when group metadata is missing", {
  obj <- collapse_lineage(scdown(example_without_group(), group_col = "group"))
  expect_error(test_composition(obj), "group_col")
  expect_error(test_signature(obj, signatures = "cytotoxic"), "group_col")
})

test_that("communication test returns empty but valid table when no LR pairs match", {
  obj <- scdown(scdown_example(), lr_pairs = data.frame(
    ligand = "NOT_A_GENE",
    receptor = "ALSO_NOT_A_GENE",
    pathway = "none",
    stringsAsFactors = FALSE
  ))
  res <- test_communication(obj, n_perm = 3)
  expect_true("communication_test_table" %in% names(res$tables))
  expect_equal(nrow(res$tables$communication_test_table), 0)
})
