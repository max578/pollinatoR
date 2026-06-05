# Smoke tests for the model-based path (simulated data) -----------------------

test_that("pollin_contribution fits and returns a treatment effect", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("marginaleffects")
  set.seed(7)
  ns <- 12L
  per <- 20L
  Site <- rep(paste0("s", seq_len(ns)), each = per)
  trt <- rep(rep(c("Cage", "Open"), each = per / 2L), ns)
  site_eff <- rnorm(ns, 0, 0.3)[as.integer(factor(Site))]
  mu <- exp(1.7 + 0.35 * (trt == "Open") + site_eff)
  df <- data.frame(y = rpois(length(mu), mu), trt = trt, Site = Site)
  sp <- pollin_spec(
    df,
    outcomes = "y",
    treatment = "trt",
    control = "Cage",
    random = "Site"
  )
  ct <- suppressWarnings(pollin_contribution(sp))
  expect_s3_class(ct, "pollin_contribution")
  expect_true("estimate" %in% names(ct$table))
  expect_true(is.finite(ct$table$estimate[1]))
})

test_that("pollin_moderation returns within and between interaction terms", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")
  set.seed(9)
  ns <- 14L
  per <- 20L
  Site <- rep(paste0("s", seq_len(ns)), each = per)
  trt <- rep(rep(c("Cage", "Open"), each = per / 2L), ns)
  bees <- rep(rpois(ns, 8), each = per) + rpois(length(Site), 1)
  site_eff <- rnorm(ns, 0, 0.3)[as.integer(factor(Site))]
  mu <- exp(
    1.6 + 0.3 * (trt == "Open") + 0.02 * bees * (trt == "Open") + site_eff
  )
  df <- data.frame(
    y = rpois(length(mu), mu),
    trt = trt,
    Site = Site,
    bees = bees
  )
  sp <- pollin_spec(
    df,
    outcomes = "y",
    treatment = "trt",
    control = "Cage",
    random = "Site",
    exposures = "bees"
  )
  md <- suppressWarnings(pollin_moderation(sp, outcome = "y"))
  expect_s3_class(md, "pollin_moderation")
  expect_true(any(grepl("within", md$table$component)))
  expect_true(any(grepl("between", md$table$component)))
})
