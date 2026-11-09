# Layer L2: pollin_population() -------------------------------------------------

pop_df <- function(seed = 21L) {
  set.seed(seed)
  d <- expand.grid(
    plant = 1:4,
    trt = c("Cage", "Open"),
    Block = sprintf("b%d", 1:8),
    stringsAsFactors = FALSE
  )
  blk <- as.integer(factor(d$Block))
  d$Band <- ifelse(blk <= 4L, "near", "far")
  # The far band yields less and depends less on pollinators.
  base <- ifelse(d$Band == "near", 3.4, 2.6)
  lift <- ifelse(d$Band == "near", 0.5, 0.15)
  d$pods <- stats::rpois(
    nrow(d),
    exp(base + lift * (d$trt == "Open") + stats::rnorm(8L, 0, 0.1)[blk])
  )
  d$plant <- NULL
  d
}
pop_arms <- c(excluded = "Cage", open = "Open")

pop_estimand <- function(population = "trial", scale = "ratio") {
  pd <- pollin_data(pop_df(), "pods", "trt", unit = "Block")
  ds <- pollin_design(pd, pop_arms, "rcbd")
  pollin_estimand(ds, "dependence", scale, "loss", population = population)
}

test_that("pollin_population records the strata, raw areas and weights", {
  p <- pollin_population("Band", c(near = 300, far = 100), label = "two bands")
  expect_s3_class(p, "pollin_population")
  expect_identical(p$by, "Band")
  expect_identical(p$area, c(near = 300, far = 100))
  expect_identical(p$weights, c(near = 0.75, far = 0.25))
  expect_identical(p$label, "two bands")
  out <- capture.output(print(p))
  expect_true(any(grepl("2 levels of Band", out)))
  expect_true(any(grepl("near 0.75", out)))
  expect_true(any(grepl("two bands", out)))
  expect_null(pollin_population("Band", c(near = 1, far = 1))$label)
})

test_that("pollin_population checks its arguments", {
  w <- c(near = 1, far = 2)
  expect_error(pollin_population(c("a", "b"), w), "`by`")
  expect_error(pollin_population(1, w), "`by`")
  expect_error(pollin_population("Band", c(1, 2)), "named by stratum")
  expect_error(pollin_population("Band", c(a = 1, a = 2)), "each name once")
  expect_error(pollin_population("Band", c(a = "1")), "numeric")
  expect_error(pollin_population("Band", c(a = 1, b = 0)), "positive")
  expect_error(pollin_population("Band", c(a = 1, b = -1)), "positive")
  expect_error(pollin_population("Band", c(a = 1, b = Inf)), "positive")
  expect_error(pollin_population("Band", w, label = 1), "label")
})

test_that("the estimand describes the population and its assumption", {
  p <- pollin_population("Band", c(near = 3, far = 1), label = "two bands")
  est <- pop_estimand(p)
  expect_identical(est$population, p)
  a <- pollin_assumptions(est)
  stmt <- a$statement[a$assumption == "population"]
  expect_match(stmt, "strata of `Band`")
  expect_match(stmt, "weighted by area")
  expect_match(stmt, "comes from the fitted model")
  expect_false(grepl("units in the trial", stmt))
  out <- capture.output(print(est))
  expect_true(any(grepl("two bands: 2 strata of Band, weighted by area", out)))

  trial <- pop_estimand()
  expect_identical(trial$population, "trial")
  expect_match(
    pollin_assumptions(trial)$statement[
      pollin_assumptions(trial)$assumption == "population"
    ],
    "units in the trial"
  )
})

test_that("the estimand refuses a population that does not fit the data", {
  expect_error(pop_estimand("field"), "\"trial\" or a `pollin_population`")
  expect_error(
    pop_estimand(pollin_population("Region", c(a = 1, b = 1))),
    "'Region' is not in the data"
  )
  expect_error(
    pop_estimand(pollin_population("Band", c(near = 1, middle = 1))),
    "levels of 'Band': far, near"
  )
  expect_error(
    pop_estimand(pollin_population("Band", c(near = 1))),
    "levels of 'Band'"
  )
})

test_that("a population estimand must be fitted by its strata", {
  p <- pollin_population("Band", c(near = 3, far = 1))
  est <- pop_estimand(p)
  expect_error(pollin_fit(est), "by = \"Band\"")
  expect_error(pollin_fit(est, by = "Block"), "by = \"Band\"")

  # A fit made without the population's strata is refused at the draws too.
  fit <- suppressWarnings(pollin_fit(pop_estimand(), by = "Block"))
  fit$estimand <- est
  expect_error(pollin_draws(fit, n = 10L), "by = \"Band\"")
})

test_that("C5: a balanced two-stratum population gives the analytic D", {
  # Both strata have the same arm means, so the ratio of population means and
  # the weighted mean of the stratum shares coincide, and both equal
  # 1 - sum(w * low) / sum(w * high).
  p <- pollin_population("Band", c(near = 1, far = 1))
  est <- pop_estimand(p)
  h <- matrix(
    c(40, 40, 36, 36),
    ncol = 2L,
    byrow = TRUE,
    dimnames = list(NULL, c("far", "near"))
  )
  l <- matrix(
    c(30, 30, 27, 27),
    ncol = 2L,
    byrow = TRUE,
    dimnames = list(NULL, c("far", "near"))
  )
  dr <- .new_pollin_draws(h, l, est, "pods", "Band", "x", "x", 2L, NULL, p)
  expect_identical(colnames(dr$estimate), c("far", "near", "population"))
  expect_equal(dr$high[, "population"], c(40, 36))
  expect_equal(dr$low[, "population"], c(30, 27))
  expect_equal(dr$estimate[, "population"], c(0.25, 0.25))
  expect_equal(
    dr$estimate[, "population"],
    1 -
      (0.5 * l[, "far"] + 0.5 * l[, "near"]) /
        (0.5 * h[, "far"] + 0.5 * h[, "near"])
  )
  expect_equal(
    dr$estimate[, "population"],
    0.5 * dr$estimate[, "far"] + 0.5 * dr$estimate[, "near"]
  )
})

test_that("C5: on unbalanced strata the population D is a ratio of means", {
  # near: high 100, low 50 (D = 0.5); far: high 20, low 18 (D = 0.1);
  # weights near 0.25, far 0.75. Population D = 1 - 26 / 40 = 0.35, while the
  # weighted mean of the stratum shares is 0.25 * 0.5 + 0.75 * 0.1 = 0.2.
  p <- pollin_population("Band", c(near = 1, far = 3))
  est <- pop_estimand(p)
  h <- matrix(c(20, 100), ncol = 2L, dimnames = list(NULL, c("far", "near")))
  l <- matrix(c(18, 50), ncol = 2L, dimnames = list(NULL, c("far", "near")))
  dr <- .new_pollin_draws(h, l, est, "pods", "Band", "x", "x", 1L, NULL, p)
  expect_equal(unname(dr$estimate[, "population"]), 0.35)
  wmean <- sum(p$weights[c("far", "near")] * dr$estimate[1L, c("far", "near")])
  expect_equal(wmean, 0.2)
  expect_gt(abs(dr$estimate[1L, "population"] - wmean), 0.1)

  # The difference scale is the weighted mean of the stratum differences.
  est_d <- pop_estimand(p, scale = "difference")
  dd <- .new_pollin_draws(h, l, est_d, "pods", "Band", "x", "x", 1L, NULL, p)
  expect_equal(unname(dd$estimate[, "population"]), 0.75 * 2 + 0.25 * 50)
})

test_that("population draws from a fit are the area-weighted stratum means", {
  skip_on_cran()
  p <- pollin_population("Band", c(near = 1, far = 3))
  fit <- suppressWarnings(pollin_fit(pop_estimand(p), by = "Band"))
  for (method in c("simulation", "bootstrap")) {
    dr <- suppressWarnings(pollin_draws(
      fit,
      n = if (method == "simulation") 500L else 9L,
      method = method,
      seed = 4L
    ))
    expect_identical(dr$population, p)
    expect_identical(colnames(dr$high), c("far", "near", "population"))
    w <- p$weights[c("far", "near")]
    expect_equal(
      dr$high[, "population"],
      as.vector(dr$high[, c("far", "near")] %*% w)
    )
    expect_equal(
      dr$low[, "population"],
      as.vector(dr$low[, c("far", "near")] %*% w)
    )
    expect_equal(
      dr$estimate[, "population"],
      1 - dr$low[, "population"] / dr$high[, "population"]
    )
    wmean <- as.vector(dr$estimate[, c("far", "near")] %*% w)
    # The near band carries the larger yield and the larger share, so the
    # ratio of population means exceeds the weighted mean of the shares.
    expect_gt(stats::median(dr$estimate[, "population"] - wmean), 0.02)
  }
  s <- summary(dr)
  expect_identical(s$level, c("far", "near", "population"))
})

test_that("pollin_effects reports the trial contrasts of a population fit", {
  skip_on_cran()
  p <- pollin_population("Band", c(near = 1, far = 3))
  fit <- suppressWarnings(pollin_fit(pop_estimand(p), by = "Band"))
  tab <- suppressWarnings(pollin_effects(fit))$table
  expect_identical(tab$Band, c("far", "near"))
  trial <- suppressWarnings(pollin_fit(pop_estimand(), by = "Band"))
  expect_equal(tab, suppressWarnings(pollin_effects(trial))$table)
})

test_that("population means need every stratum and a free column name", {
  p <- pollin_population("Band", c(near = 1, far = 3))
  h <- matrix(c(2, 3), ncol = 2L, dimnames = list(NULL, c("far", "other")))
  expect_error(.population_means(h, h, p), "lack the strata near")
  h2 <- matrix(
    c(2, 3, 4),
    ncol = 3L,
    dimnames = list(NULL, c("far", "near", "population"))
  )
  expect_error(.population_means(h2, h2, p), "already named")
  expect_error(.population_means(h, h, list()), "pollin_population")
})
