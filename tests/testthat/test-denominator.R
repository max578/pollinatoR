# A response counted out of a denominator: binomial fits on both engines -------

# Seed set: blocks x two arms x plants, seeds out of 20 to 40 ovules, logit p.
den_df <- function(seed = 3L, n_block = 8L, n_plant = 6L, b0 = -0.4, b1 = 0.9) {
  set.seed(seed)
  d <- expand.grid(
    plant = seq_len(n_plant),
    trt = c("Cage", "Open"),
    Block = sprintf("b%02d", seq_len(n_block)),
    stringsAsFactors = FALSE
  )
  u <- stats::rnorm(n_block, 0, 0.3)
  d$u <- u[as.integer(factor(d$Block))]
  d$ovules <- sample(20:40, nrow(d), replace = TRUE)
  d$seeds <- stats::rbinom(
    nrow(d),
    d$ovules,
    stats::plogis(b0 + b1 * (d$trt == "Open") + d$u)
  )
  d$plant <- NULL
  d
}

den_estimand <- function(d = den_df(), scale = "difference", basis = "loss") {
  pd <- pollin_data(d, "seeds", "trt", unit = "Block", denominator = "ovules")
  ds <- pollin_design(pd, c(excluded = "Cage", open = "Open"), "rcbd")
  pollin_estimand(ds, "dependence", scale, basis)
}

# Hand G-computation: each fitted row predicted under each arm, averaged.
den_hand <- function(fit) {
  m <- fit$fits$seeds
  d <- fit$data$seeds
  arm <- function(level) {
    nd <- d
    nd$trt <- factor(rep(level, nrow(d)), levels = levels(d$trt))
    mean(stats::predict(m, newdata = nd, type = "response"))
  }
  c(high = arm("Open"), low = arm("Cage"))
}

test_that("a denominator fits the count as binomial on glmmTMB", {
  skip_on_cran()
  fit <- pollin_fit(den_estimand())
  expect_identical(fit$denominator, "ovules")
  expect_identical(fit$families$seeds$family, "binomial")
  expect_identical(fit$families$seeds$link, "logit")
  expect_identical(
    deparse(stats::formula(fit$fits$seeds)[[2L]]),
    "cbind(seeds, ovules - seeds)"
  )
  expect_identical(nrow(fit$data$seeds), 96L)
  expect_output(print(fit), "trials    : ovules")

  hand <- den_hand(fit)
  expect_true(all(hand > 0 & hand < 1))
  m <- .arm_means(fit$fits$seeds, fit$data$seeds, "trt", "Open", "Cage", NULL)
  expect_equal(unname(c(m$high, m$low)), unname(hand), tolerance = 1e-8)

  delta <- suppressWarnings(pollin_effects(fit))$table
  # marginaleffects reads the conditional modes at TMB's last evaluated
  # parameters, which can sit a step away from the optimum.
  expect_equal(delta$estimate, hand[["high"]] - hand[["low"]], tolerance = 1e-4)
  expect_identical(delta$distribution, "binomial")
  expect_true(delta$std.error > 0)

  boot <- suppressWarnings(
    pollin_effects(fit, interval = "bootstrap", n_boot = 20L, seed = 1L)
  )$table
  expect_equal(boot$estimate, hand[["high"]] - hand[["low"]], tolerance = 1e-8)
  expect_true(boot$conf.low < boot$estimate && boot$estimate < boot$conf.high)

  sim <- pollin_draws(fit, n = 400L, seed = 1L)
  expect_identical(dim(sim$estimate), c(400L, 1L))
  expect_true(all(sim$high > 0 & sim$high < 1))
  expect_lt(
    abs(stats::median(sim$estimate) - (hand[["high"]] - hand[["low"]])),
    delta$std.error / 2
  )
  ct <- suppressWarnings(pollin_contribution(den_estimand()))$table
  expect_equal(ct$estimate, delta$estimate, tolerance = 1e-8)
  expect_identical(ct$distribution, "binomial")

  bt <- suppressWarnings(pollin_draws(
    fit,
    n = 10L,
    method = "bootstrap",
    seed = 1L
  ))
  expect_true(all(bt$low > 0 & bt$low < 1))
})

test_that("ratio scales and value run on draws of a binomial fit", {
  skip_on_cran()
  fit <- pollin_fit(den_estimand(scale = "ratio"))
  hand <- den_hand(fit)
  m <- .arm_means(fit$fits$seeds, fit$data$seeds, "trt", "Open", "Cage", NULL)
  expect_equal(
    .effect_from_means(m, "ratio", "loss")[[1L]],
    1 - hand[["low"]] / hand[["high"]],
    tolerance = 1e-8
  )
  dr <- pollin_draws(fit, n = 300L, seed = 2L)
  mk <- pollin_market(1000, 500, "t", "AUD", 2024, "Synthetic", crop = "canola")
  v <- pollin_value(dr, mk, channel = "dependence")
  expect_s3_class(v, "pollin_value")
  expect_equal(
    v$draws,
    (1 - dr$low[, "overall"] / dr$high[, "overall"]) * 500 * 1000
  )
})

test_that("the estimate recovers the planted row-averaged proportion contrast", {
  skip_on_cran()
  d <- den_df(seed = 11L, n_block = 40L, n_plant = 20L)
  truth <- mean(
    stats::plogis(-0.4 + 0.9 + d$u) - stats::plogis(-0.4 + d$u)
  )
  tab <- suppressWarnings(pollin_effects(pollin_fit(den_estimand(d))))$table
  expect_lt(abs(tab$estimate - truth), 3 * tab$std.error)
  expect_lt(abs(tab$estimate - truth), 0.02)
})

test_that("a denominator refuses other families and impossible counts", {
  skip_on_cran()
  expect_error(
    pollin_fit(den_estimand(), family = glmmTMB::nbinom2()),
    "gives 'nbinom2' for 'seeds'.*denominator 'ovules'.*drop `denominator`"
  )
  expect_error(
    pollin_fit(den_estimand(), family = list(seeds = stats::poisson())),
    "gives 'poisson'"
  )
  d <- den_df()
  d$seeds[c(1L, 5L)] <- d$ovules[c(1L, 5L)] + 1L
  expect_error(
    pollin_fit(den_estimand(d)),
    "'seeds' has 2 rows outside 0 to the denominator 'ovules'"
  )
  d <- den_df()
  d$seeds[3L] <- 2.5
  expect_error(pollin_fit(den_estimand(d)), "must be whole numbers")
  d <- den_df()
  d$ovules[c(2L, 4L, 6L)] <- 0
  d$seeds[c(2L, 4L, 6L)] <- 0
  expect_error(
    pollin_fit(den_estimand(d)),
    "'seeds' has 3 rows with the denominator 'ovules' equal to 0"
  )
  d$ovules[c(2L, 4L, 6L)] <- NA
  fit <- pollin_fit(den_estimand(d))
  expect_identical(nrow(fit$data$seeds), 93L)
})

test_that("a binomial fit with a denominator is the binomial family given", {
  skip_on_cran()
  fit <- pollin_fit(den_estimand(), family = stats::binomial(link = "probit"))
  expect_identical(fit$families$seeds$link, "probit")
})

test_that("trials expand to one 0/1 row each, under a row limit", {
  d <- data.frame(y = c(2, 0, 3), n = c(3, 2, 3), g = c("a", "b", "c"))
  ex <- .flexy_expand_trials(d, "y", "n")
  expect_identical(ex$rows, c(1L, 1L, 1L, 2L, 2L, 3L, 3L, 3L))
  expect_identical(ex$data$y, c(1L, 1L, 0L, 0L, 0L, 1L, 1L, 1L))
  expect_identical(ex$data$g, d$g[ex$rows])
  expect_identical(as.vector(tapply(ex$data$y, ex$rows, sum)), c(2L, 0L, 3L))
  op <- options(pollinatoR.flexy_max_trial_rows = 7)
  on.exit(options(op), add = TRUE)
  expect_error(
    .flexy_expand_trials(d, "y", "n"),
    "denominator 'n' gives 8 rows, above the limit of 7.*glmmTMB"
  )
})

test_that("flexyBayes fits a denominator on one row per trial", {
  skip_flexy()
  op <- flexy_options()
  on.exit(options(op), add = TRUE)
  est <- den_estimand()
  op_cap <- options(pollinatoR.flexy_max_trial_rows = 100)
  expect_error(pollin_fit(est, engine = "flexyBayes"), "above the limit")
  options(op_cap)
  fit <- pollin_fit(est, engine = "flexyBayes")
  d <- fit$data$seeds
  rows <- fit$flexy$trial_rows$seeds
  expect_identical(nrow(d), 96L)
  expect_identical(length(rows), as.integer(sum(d$ovules)))
  expect_identical(fit$families$seeds$family, "binomial")

  # Original-row predictor by hand, checked against the expanded rows -------
  latent <- .flexy_latent(fit$fits$seeds, 2L, 5L)
  lv <- flexyBayes::fb_model_terms(fit$fits$seeds)$random_terms[[1L]]$var_levels
  eta <- function(arm, k) {
    latent["(Intercept):1", k] +
      (arm == "Open") * latent["trtOpen:1", k] +
      latent[paste0("Block:", match(d$Block, lv)), k]
  }
  m <- .flexy_arm_draws(fit, "seeds", 2L, 5L, NULL)
  for (k in 1:2) {
    own <- eta(as.character(d$trt), k)
    expect_equal(
      unname(own[rows]),
      unname(latent[paste0("Predictor:", seq_along(rows)), k]),
      tolerance = 1e-10
    )
    expect_equal(
      unname(m$high[k, 1L]),
      mean(stats::plogis(eta("Open", k))),
      tolerance = 1e-10
    )
    expect_equal(
      unname(m$low[k, 1L]),
      mean(stats::plogis(eta("Cage", k))),
      tolerance = 1e-10
    )
  }

  b <- pollin_effects(fit, n_draws = 1000L, seed = 1L)$table
  g <- suppressWarnings(pollin_effects(pollin_fit(est)))$table
  expect_identical(b$distribution, "binomial")
  expect_lt(abs(b$estimate - g$estimate), 0.25 * g$std.error)
})
