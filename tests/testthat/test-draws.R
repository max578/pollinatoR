# Layer L4b: pollin_draws() ----------------------------------------------------
#
# Gate (fixed 2026-11-09, before the first run). With n = 20000 seeded
# simulation draws, on the lnRR scale and for each column of the draws:
#   |2.5% draw quantile  - delta conf.low|  <= 0.03 x delta interval width
#   |97.5% draw quantile - delta conf.high| <= 0.03 x delta interval width
#   |draw median         - delta estimate|  <= 0.03 x delta interval width
#   draw SD / delta standard error in [0.97, 1.03].
# The Monte-Carlo standard error of a 2.5% quantile at n = 20000 is about
# 0.005 interval widths, so 0.03 leaves room for the curvature of log(rr) in
# the fixed effects while catching a misplaced or mis-scaled draw.
# Bootstrap draws must reproduce the bootstrap interval and standard error of
# pollin_effects() exactly (identical numbers) for the same seed and n.

dr_df <- function(seed = 5L) {
  set.seed(seed)
  d <- expand.grid(
    plant = 1:4,
    trt = c("Cage", "Open"),
    Block = sprintf("b%d", 1:8),
    stringsAsFactors = FALSE
  )
  blk <- as.integer(factor(d$Block))
  d$Region <- ifelse(blk <= 4L, "north", "south")
  d$pods <- stats::rpois(
    nrow(d),
    exp(3 + 0.4 * (d$trt == "Open") + stats::rnorm(8L, 0, 0.2)[blk])
  )
  d$mass <- stats::rgamma(
    nrow(d),
    shape = 20,
    rate = 20 / exp(3 + 0.2 * (d$trt == "Open"))
  )
  d$plant <- NULL
  d
}
dr_arms <- c(excluded = "Cage", open = "Open")

dr_fit <- function(
  scale = "lnrr",
  basis = "loss",
  response = "pods",
  by = NULL
) {
  pd <- pollin_data(dr_df(), response, "trt", unit = "Block")
  ds <- pollin_design(pd, dr_arms, "rcbd")
  suppressWarnings(pollin_fit(
    pollin_estimand(ds, "dependence", scale, basis),
    by = by
  ))
}

test_that("simulation draws match the delta interval on the lnRR scale", {
  skip_on_cran()
  for (by in list(NULL, "Region")) {
    fit <- dr_fit(by = by)
    tab <- suppressWarnings(pollin_effects(fit))$table
    dr <- pollin_draws(fit, n = 20000L, seed = 11L)
    expect_identical(dim(dr$estimate), c(20000L, nrow(tab)))
    for (j in seq_len(nrow(tab))) {
      z <- dr$estimate[, j]
      width <- tab$conf.high[j] - tab$conf.low[j]
      q <- stats::quantile(z, c(0.025, 0.5, 0.975), names = FALSE)
      expect_lte(abs(q[1L] - tab$conf.low[j]), 0.03 * width)
      expect_lte(abs(q[3L] - tab$conf.high[j]), 0.03 * width)
      expect_lte(abs(q[2L] - tab$estimate[j]), 0.03 * width)
      ratio <- stats::sd(z) / tab$std.error[j]
      expect_gte(ratio, 0.97)
      expect_lte(ratio, 1.03)
    }
  }
})

test_that("bootstrap draws reproduce the bootstrap interval exactly", {
  skip_on_cran()
  for (by in list(NULL, "Region")) {
    fit <- dr_fit(scale = "ratio", by = by)
    ef <- suppressWarnings(pollin_effects(
      fit,
      interval = "bootstrap",
      n_boot = 19L,
      seed = 7L
    ))$table
    dr <- suppressWarnings(pollin_draws(
      fit,
      n = 19L,
      method = "bootstrap",
      seed = 7L
    ))
    expect_identical(dr$method, "bootstrap")
    a <- (1 - 0.95) / 2
    ci <- apply(
      dr$estimate,
      2L,
      stats::quantile,
      probs = c(a, 1 - a),
      names = FALSE
    )
    expect_identical(unname(ci[1L, ]), ef$conf.low)
    expect_identical(unname(ci[2L, ]), ef$conf.high)
    expect_identical(unname(apply(dr$estimate, 2L, stats::sd)), ef$std.error)
  }
})

test_that("at the estimated fixed effects the simulation rule is .arm_means", {
  skip_on_cran()
  for (by in list(NULL, "Region")) {
    fit <- dr_fit(by = by)
    model <- fit$fits$pods
    d <- fit$data$pods
    beta <- glmmTMB::fixef(model)$cond
    am <- .arm_means(model, d, "trt", "Open", "Cage", by)
    for (arm in c("Open", "Cage")) {
      p <- .draws_arm_linear(model, d, "trt", arm, beta)
      mu <- stats::family(model)$linkinv(as.vector(p$x %*% beta) + p$re)
      g <- if (is.null(by)) rep("all", nrow(d)) else as.character(d[[by]])
      got <- tapply(mu, g, mean)
      want <- if (arm == "Open") am$high else am$low
      expect_equal(unname(got), unname(want), tolerance = 1e-10)
    }
  }
})

test_that("the draws carry the contract fields and form the estimand", {
  skip_on_cran()
  fit <- dr_fit(scale = "ratio", basis = "loss")
  dr <- pollin_draws(fit, n = 50L, seed = 1L)
  expect_s3_class(dr, "pollin_draws")
  expect_identical(
    names(dr),
    c(
      "high",
      "low",
      "estimate",
      "estimand",
      "response",
      "by",
      "method",
      "engine",
      "n",
      "seed",
      "population"
    )
  )
  expect_identical(colnames(dr$high), "overall")
  expect_identical(dim(dr$high), c(50L, 1L))
  expect_equal(dr$estimate, 1 - dr$low / dr$high)
  expect_identical(dr$method, "simulation")
  expect_identical(dr$engine, "glmmTMB")
  expect_identical(dr$response, "pods")
  expect_null(dr$by)
  expect_null(dr$population)
  expect_identical(dr$n, 50L)
  expect_identical(dr$seed, 1L)
  expect_identical(as.matrix(dr), dr$estimate)

  by_dr <- pollin_draws(dr_fit(scale = "difference", by = "Region"), n = 30L)
  expect_identical(colnames(by_dr$estimate), c("north", "south"))
  expect_equal(by_dr$estimate, by_dr$high - by_dr$low)
  expect_identical(by_dr$by, "Region")
})

test_that("a seed reproduces the draws and restores the random state", {
  skip_on_cran()
  fit <- dr_fit()
  set.seed(99L)
  before <- .Random.seed
  a <- pollin_draws(fit, n = 20L, seed = 3L)
  expect_identical(.Random.seed, before)
  b <- pollin_draws(fit, n = 20L, seed = 3L)
  expect_identical(a$estimate, b$estimate)
  c <- pollin_draws(fit, n = 20L, seed = 4L)
  expect_false(identical(a$estimate, c$estimate))
  set.seed(5L)
  d1 <- pollin_draws(fit, n = 20L)
  set.seed(5L)
  d2 <- pollin_draws(fit, n = 20L)
  expect_identical(d1$estimate, d2$estimate)
})

test_that("summary and print name the estimand, scale, basis and method", {
  skip_on_cran()
  dr <- pollin_draws(
    dr_fit(scale = "ratio", basis = "gain"),
    n = 200L,
    seed = 2L
  )
  s <- summary(dr)
  expect_s3_class(s, "data.frame")
  expect_identical(
    names(s),
    c(
      "response",
      "level",
      "median",
      "conf.low",
      "conf.high",
      "family",
      "scale",
      "basis",
      "method",
      "n"
    )
  )
  expect_identical(s$level, "overall")
  expect_identical(s$basis, "gain")
  expect_identical(
    s$median,
    stats::quantile(dr$estimate[, 1L], 0.5, names = FALSE)
  )
  expect_lt(s$conf.low, s$median)
  expect_lt(summary(dr, conf_level = 0.5)$conf.high, s$conf.high)
  expect_error(summary(dr, conf_level = 2), "conf_level")
  out <- capture.output(print(dr))
  expect_true(any(grepl("dependence \\(E2\\)", out)))
  expect_true(any(grepl("basis +: gain", out)))
  expect_true(any(grepl("simulation of the fixed effects", out)))
  expect_true(any(grepl("fixed-effect uncertainty only", out)))
})

test_that("pollin_draws checks its arguments", {
  skip_on_cran()
  fit <- dr_fit()
  expect_error(pollin_draws(list()), "pollin_fit")
  expect_error(pollin_draws(fit, method = "posterior"), "method")
  expect_error(pollin_draws(fit, n = 1L), "`n`")
  expect_error(pollin_draws(fit, n = 10.5), "`n`")
  expect_error(pollin_draws(fit, seed = "a"), "seed")
  expect_error(pollin_draws(fit, response = "mass"), "pods")

  two <- dr_fit(response = c("pods", "mass"))
  expect_error(pollin_draws(two, n = 10L), "several responses")
  dm <- pollin_draws(two, n = 10L, response = "mass", seed = 1L)
  expect_identical(dm$response, "mass")
})

test_that("the constructor validates the shapes of the draws", {
  est <- pollin_estimand(
    pollin_design(
      pollin_data(dr_df(), "pods", "trt", unit = "Block"),
      dr_arms,
      "rcbd"
    ),
    "dependence",
    "lnrr",
    "loss"
  )
  h <- matrix(c(2, 3, 4), ncol = 1L, dimnames = list(NULL, "overall"))
  l <- matrix(c(1, 1, 2), ncol = 1L, dimnames = list(NULL, "overall"))
  ok <- .new_pollin_draws(
    h,
    l,
    est,
    "pods",
    NULL,
    "posterior",
    "other",
    3L,
    NULL,
    NULL
  )
  expect_equal(ok$estimate, log(h / l))
  expect_identical(ok$engine, "other")
  build <- function(...) {
    args <- list(
      high = h,
      low = l,
      estimand = est,
      response = "pods",
      by = NULL,
      method = "posterior",
      engine = "other",
      n = 3L,
      seed = NULL,
      population = NULL
    )
    args[names(list(...))] <- list(...)
    do.call(.new_pollin_draws, args)
  }
  expect_error(build(high = c(2, 3, 4)), "numeric matrix")
  expect_error(build(high = unname(h)), "named")
  expect_error(build(low = l[1:2, , drop = FALSE]), "same dimensions")
  expect_error(build(n = 4L), "`n`")
  expect_error(build(estimand = list()), "pollin_estimand")
  expect_error(build(response = 1), "response")
  expect_error(build(method = NULL), "method")
  expect_error(build(engine = c("a", "b")), "engine")
})
