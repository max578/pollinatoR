test_that("the terms add up to the yield log ratio in every replicate", {
  skip_on_cran()
  fit <- decompose_fit(decompose_fixture(10L, 6L, seed = 2L))
  r <- pollin_decompose(fit, "yield", c("pods", "seeds", "weight"),
    n_boot = 20L, seed = 1L)
  parts <- c("pods", "seeds", "weight", "interaction")
  expect_equal(rowSums(r$draws[, parts]), r$draws[, "yield"],
    tolerance = 1e-12)
  tab <- r$table
  expect_identical(tab$term, c(parts, "yield"))
  expect_equal(sum(tab$log_ratio[1:4]), tab$log_ratio[5], tolerance = 1e-12)
  expect_equal(sum(tab$share[1:4]), 1, tolerance = 1e-12)
  expect_lt(r$product_check, 1e-12)
})

# Analytic oracle: with independent multiplicative components the planted log
# ratio of yield is 0.10 + 0.30 - 0.10 = 0.30 and the interaction is 0; each
# estimate must lie within 3 bootstrap SE of its planted value.
test_that("planted multiplicative effects are recovered on the log scale", {
  skip_on_cran()
  fit <- decompose_fit(decompose_fixture(30L, 15L, seed = 11L))
  r <- pollin_decompose(fit, "yield", c("pods", "seeds", "weight"),
    n_boot = 30L, seed = 1L)
  tab <- r$table
  truth <- c(pods = 0.10, seeds = 0.30, weight = -0.10, interaction = 0,
    yield = 0.30)
  se <- apply(r$draws, 2L, stats::sd)
  expect_true(all(abs(tab$log_ratio - truth) < 3 * se))
  expect_lt(tab$conf.low[tab$term == "interaction"], 0)
  expect_gt(tab$conf.high[tab$term == "interaction"], 0)
  expect_identical(tab$compensating[1:3], c(FALSE, FALSE, TRUE))
})

test_that("the resampling is reproducible and leaves the session's seed", {
  skip_on_cran()
  fit <- decompose_fit(decompose_fixture(8L, 5L, seed = 3L))
  set.seed(99)
  before <- .Random.seed
  a <- pollin_decompose(fit, "yield", c("pods", "seeds", "weight"),
    n_boot = 5L, seed = 7L)
  expect_identical(.Random.seed, before)
  b <- pollin_decompose(fit, "yield", c("pods", "seeds", "weight"),
    n_boot = 5L, seed = 7L)
  expect_identical(a$draws, b$draws)
  expect_output(print(a), "compensating: weight")
})

test_that("a unit drawn twice enters the refit as two units", {
  d <- data.frame(Site = c("a", "a", "b"), y = 1:3)
  rows <- split(seq_len(3L), factor(d$Site, c("a", "b")))
  out <- .decomp_resampled_data(d, rows, c(1L, 1L, 2L), "Site")
  expect_identical(nlevels(out$Site), 3L)
  expect_identical(out$y, c(1L, 2L, 1L, 2L, 3L))
  expect_identical(.decomp_groups(c("trt", "(1 | Site)", "(trt | f/p)")),
    c("Site", "f", "p"))
})

test_that("inputs that cannot be decomposed are refused", {
  skip_on_cran()
  d <- decompose_fixture(8L, 5L, seed = 4L)
  fit <- decompose_fit(d)
  expect_error(pollin_decompose(fit, "yeld", c("pods", "seeds")), "`yield`")
  expect_error(pollin_decompose(fit, "yield", "pods"), "at least two")
  expect_error(pollin_decompose(fit, "yield", c("pods", "yield")),
    "at least two")
  d$yield <- d$yield * stats::runif(nrow(d), 0.9, 1.1)
  expect_warning(
    pollin_decompose(decompose_fit(d), "yield", c("pods", "seeds", "weight"),
      n_boot = 3L, seed = 1L),
    "not a constant multiple"
  )
  pd <- pollin_data(d, c("yield", "pods", "seeds"), "trt", unit = "Site")
  ds <- pollin_design(pd, c(excluded = "Cage", open = "Open"), "rcbd")
  est <- pollin_estimand(ds, "dependence", "ratio", "loss")
  fit_id <- pollin_fit(est, family = list(
    yield = stats::gaussian(), pods = glmmTMB::nbinom2(),
    seeds = stats::Gamma(link = "log")))
  expect_error(
    pollin_decompose(fit_id, "yield", c("pods", "seeds"), n_boot = 3L),
    "log link"
  )
})

# Oracle `_golden/golden_v0.5.0_decompose.rds`, written once from the 0.5.0
# code by `validation/golden/make_golden_v0.5.0_decompose.R`.

test_that("the decomposition matches 0.5.0", {
  skip_on_cran()
  gold <- readRDS(golden_path_v050_decompose())
  now <- golden_capture_v050_decompose()
  expect_identical(names(now), names(gold))
  for (k in names(gold)) {
    expect_equal(now[[k]], gold[[k]], tolerance = 1e-10, info = k,
      ignore_attr = TRUE)
  }
})
