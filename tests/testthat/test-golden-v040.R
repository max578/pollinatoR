# Golden values: numbers carried from 0.4.0 -----------------------------------
#
# Oracle `_golden/golden_v0.4.0.rds`, written once from the 0.4.0 code by
# `validation/golden/make_golden_v0.4.0.R`; model-based tolerance 1e-10.

test_that("joint delta effects, draws, value and total match 0.4.0", {
  skip_on_cran()
  gold <- readRDS(golden_path_v040())
  now <- suppressWarnings(golden_capture_v040())
  expect_identical(names(now), names(gold))
  for (k in names(gold)) {
    expect_equal(now[[k]], gold[[k]], tolerance = 1e-10, info = k)
  }
})
