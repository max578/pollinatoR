# Golden values: numbers carried from 0.2.0 -----------------------------------
#
# Oracle `_golden/golden_v0.2.0.rds`, written once from the 0.2.0 code by
# `validation/golden/make_golden_v0.2.0.R`; model-based tolerance 1e-10.
# Since 0.4.0 the simulation draws come from the joint covariance of the fixed
# and random effects, so the outputs built on them are guarded by the 0.4.0
# golden values; the bootstrap draws are unchanged and still match 0.2.0.

test_that("bootstrap draws match 0.2.0", {
  skip_on_cran()
  gold <- readRDS(golden_path_v020())
  now <- suppressWarnings(golden_capture_v020())
  expect_identical(names(now), names(gold))
  expect_equal(
    now$draws_bootstrap,
    gold$draws_bootstrap,
    tolerance = 1e-10
  )
})
