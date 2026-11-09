# Golden values: numbers carried from 0.4.0 -----------------------------------
#
# Oracle `_golden/golden_v0.4.0.rds`, written once from the 0.4.0 code by
# `validation/golden/make_golden_v0.4.0.R`; model-based tolerance 1e-10.

test_that("joint delta effects, draws, value and total match 0.4.0", {
  skip_on_cran()
  gold <- readRDS(golden_path_v040())
  now <- golden_total_only(suppressWarnings(golden_capture_v040()))
  expect_identical(names(now), names(gold))
  for (k in names(gold)) {
    expect_equal(now[[k]], gold[[k]], tolerance = 1e-10, info = k)
  }
})

# Oracle `_golden/golden_v0.4.0_moderation.rds`, written once from the 0.4.0
# code by `validation/golden/make_golden_v0.4.0_moderation.R`. Columns added
# after 0.4.0 are allowed; the 0.4.0 columns must not move.

test_that("moderation matches 0.4.0", {
  skip_on_cran()
  gold <- readRDS(golden_path_v040_moderation())
  now <- suppressWarnings(suppressMessages(golden_capture_v040_moderation()))
  expect_identical(names(now), names(gold))
  for (k in names(gold)) {
    expect_identical(now[[k]]$outcome, gold[[k]]$outcome, info = k)
    expect_equal(
      now[[k]]$table[names(gold[[k]]$table)],
      gold[[k]]$table,
      tolerance = 1e-10,
      info = k
    )
  }
})

# Oracle `_golden/golden_v0.4.0_cate.rds`, written once from the 0.4.0 code
# (2a5748b) by `validation/golden/make_golden_v0.4.0_cate.R`. Columns added
# after 0.4.0 are allowed; the 0.4.0 columns must not move.

test_that("CATE matches 0.4.0", {
  skip_on_cran()
  gold <- readRDS(golden_path_v040_cate())
  now <- suppressWarnings(golden_capture_v040_cate())
  expect_identical(names(now), names(gold))
  for (k in names(gold)) {
    expect_equal(now[[k]][names(gold[[k]])], gold[[k]], tolerance = 1e-10,
      info = k)
  }
})
