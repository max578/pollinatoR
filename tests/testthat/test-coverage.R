# Edge branches, errors, and print methods ------------------------------------

test_that("print methods run", {
  expect_output(print(pollin_rii(
    c(8, 9, 10),
    c(6, 7, 5),
    n_boot = 50,
    seed = 1
  )))
  df <- data.frame(
    y = 1:4,
    trt = rep(c("Cage", "Open"), 2),
    Site = rep(c("a", "b"), each = 2)
  )
  sp <- pollin_spec(
    df,
    "y",
    "trt",
    "Cage",
    "Site",
    exposures = "y",
    adjust = "y",
    group = "Site"
  )
  expect_output(print(sp))
  expect_output(print(pollin_validate(sp)))
  expect_output(print(pollin_compensation(data.frame(
    a = c(1, 2, 3),
    b = c(3, 2, 1)
  ))))
  g <- pollin_guard_spurious(
    c(10, 20, 30, 40),
    rep(5, 4),
    n_perm = 20,
    seed = 1
  )
  expect_output(print(g))
})

test_that(".pick_family covers each branch", {
  expect_equal(pollinatoR:::.pick_family(c(0L, 1L, 2L, 3L))$family, "nbinom2")
  expect_equal(pollinatoR:::.pick_family(c(0.1, 0.5, 0.9))$family, "beta")
  expect_equal(pollinatoR:::.pick_family(c(0, 1.5, 2.5))$family, "tweedie")
  expect_equal(pollinatoR:::.pick_family(c(1.2, 3.4, 5.6))$family, "Gamma")
  expect_s3_class(pollinatoR:::.pick_family(numeric(0)), "family")
})

test_that("pollin_index handles a zero exclusion mean (lnRR undefined)", {
  idx <- pollin_index(c(5, 6, 7), c(0, 0, 0), n_boot = 50, seed = 1)
  expect_true(is.na(idx$estimate[idx$index == "lnRR"]))
})

test_that("input-validation errors fire", {
  expect_error(pollin_within_between(1:3, c("a", "b")))
  expect_error(pollin_guard_spurious(1:3, 1:2))
  expect_error(pollin_compensation(data.frame(a = 1:3)))
  expect_error(pollin_moderation(
    pollin_spec(
      data.frame(y = 1:2, trt = c("Cage", "Open"), Site = c("a", "b")),
      "y",
      "trt",
      "Cage",
      "Site"
    ),
    outcome = "y"
  ))
})

test_that("validate flags NA columns and a decoupled guard runs", {
  df <- data.frame(
    y = c(1, NA, 3, 4),
    trt = rep(c("Cage", "Open"), 2),
    Site = rep(c("a", "b"), each = 2)
  )
  sp <- pollin_spec(df, "y", "trt", "Cage", "Site")
  v <- pollin_validate(sp)
  expect_true(any(grepl("NA:y", v$table$check)))
  g <- pollin_guard_spurious(
    c(10, 20, 30, 40),
    c(5, 6, 7, 8),
    decouple = c(1, 2, 3, 4),
    n_perm = 50,
    seed = 1
  )
  expect_false(is.null(g$decoupled))
})
