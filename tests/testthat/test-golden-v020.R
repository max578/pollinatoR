# Golden values: numbers carried from 0.2.0 -----------------------------------
#
# Oracle `_golden/golden_v0.2.0.rds`, written once from the 0.2.0 code by
# `validation/golden/make_golden_v0.2.0.R`; model-based tolerance 1e-10.

test_that("draws, population, value and total match 0.2.0", {
  skip_on_cran()
  gold <- readRDS(golden_path_v020())
  now <- suppressWarnings(golden_capture_v020())
  expect_identical(names(now), names(gold))
  for (k in names(gold)) {
    expect_equal(now[[k]], gold[[k]], tolerance = 1e-10, info = k)
  }
})
