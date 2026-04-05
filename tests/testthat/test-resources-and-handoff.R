test_that("built-ins and custom resources are available", {
  obj <- scdown(
    scdown_example(),
    signatures = list(my_state = c("CD3D", "IL7R")),
    marker_panels = list(my_panel = c("CD3D", "TRAC"))
  )
  expect_true("cytotoxic" %in% list_signatures())
  expect_true("my_state" %in% list_signatures(obj))
  expect_true("NKG7" %in% get_signature("cytotoxic"))
  expect_true("CD3D" %in% get_signature("my_state", obj))
})

test_that("recommended pipeline and scope cover all public scenarios", {
  scenarios <- c("generic", "pbmc", "tumor_immune", "patient_comparison")
  for (scenario in scenarios) {
    tbl <- recommended_pipeline(scenario)
    expect_true(all(c("stage", "goal", "recommended_tools", "scdown_role") %in% names(tbl)))
    expect_true(any(tbl$scdown_role == "core"))

    scope <- scope_table(scenario)
    expect_true(all(c("stage", "use_scdown", "specialist_tools", "note", "scenario_note") %in% names(scope)))
    expect_true(any(scope$use_scdown == "handoff"))
  }
})

test_that("handoff export returns target-ready structures", {
  obj <- collapse_lineage(scdown(scdown_example()))

  muscat_export <- export_for_handoff(obj, target = "muscat")
  expect_true(all(c("expr", "cells", "sample_metadata", "embedding", "guidance") %in% names(muscat_export)))
  expect_equal(muscat_export$cluster_col, obj$annotation_col)

  dreamlet_export <- export_for_handoff(obj, target = "dreamlet")
  expect_true("sce" %in% names(dreamlet_export))

  milor_export <- export_for_handoff(obj, target = "miloR")
  expect_true("sce" %in% names(milor_export))

  cellchat_export <- export_for_handoff(obj, target = "CellChat")
  expect_true("cellchat_input" %in% names(cellchat_export))
  expect_equal(cellchat_export$cellchat_input$group.by, obj$annotation_col)
})

test_that("resource lookups fail clearly for unknown names", {
  obj <- scdown(scdown_example())
  expect_error(get_signature("does_not_exist", obj), "Unknown signature")
  expect_error(plot_signature(obj, signatures = "does_not_exist"), "Unknown signatures")
  expect_error(check_annotation(obj, panels = "does_not_exist"), "Unknown marker_panels")
})
