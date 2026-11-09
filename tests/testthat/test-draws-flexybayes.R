# Layer L4b: posterior draws of a flexyBayes fit --------------------------------

# A flexyBayes-shaped fit whose latent draws are the estimates of a glmmTMB
# fit with a site slope: site intercepts on `Site`, and the slope as the
# `Open` cell of `pollin_Site_trt` (the `Cage` cell is zero).
glmm_as_flexy <- function() {
  est <- fx_estimand("multi_site")
  g <- suppressWarnings(pollin_fit(est))
  m <- g$fits$pods
  d <- g$data$pods
  re <- glmmTMB::ranef(m)$cond$Site
  sites <- rownames(re)
  d$pollin_Site_trt <- interaction(d$Site, d$trt, drop = TRUE, sep = ":")
  cells <- levels(d$pollin_Site_trt)
  cell_value <- ifelse(
    grepl(":Open$", cells),
    re[sub(":.*$", "", cells), "trtOpen"],
    0
  )
  beta <- glmmTMB::fixef(m)$cond
  eta <- stats::predict(m, type = "link")
  latent <- matrix(c(
    beta,
    re[, "(Intercept)"],
    cell_value,
    eta
  ))
  rownames(latent) <- c(
    paste0(names(beta), ":1"),
    paste0("Site:", seq_along(sites)),
    paste0("pollin_Site_trt:", seq_along(cells)),
    paste0("Predictor:", seq_along(eta))
  )
  fake <- g
  fake$engine <- "flexyBayes"
  fake$data$pods <- d
  fake$flexy <- list(
    fixed = "trt",
    random = c("Site", "pollin_Site_trt"),
    cells = list(list(term = "pollin_Site_trt", key = "Site"))
  )
  terms <- list(
    random_terms = list(
      list(var = "Site", var_levels = sites),
      list(var = "pollin_Site_trt", var_levels = cells)
    )
  )
  list(glmm = g, fake = fake, latent = latent, terms = terms)
}

test_that("the posterior rule is the G-computation rule of the glmmTMB route", {
  skip_on_cran()
  skip_if_not_installed("flexyBayes")
  x <- glmm_as_flexy()
  local_mocked_bindings(.flexy_latent = function(model, n, seed) x$latent)
  local_mocked_bindings(
    fb_model_terms = function(fit) x$terms,
    .package = "flexyBayes"
  )
  m <- .flexy_arm_draws(x$fake, "pods", 1L, NULL, NULL)
  d <- x$glmm$data$pods
  ref <- .arm_means(x$glmm$fits$pods, d, "trt", "Open", "Cage", NULL)
  expect_equal(as.vector(m$high), as.vector(ref$high), tolerance = 1e-10)
  expect_equal(as.vector(m$low), as.vector(ref$low), tolerance = 1e-10)

  d$Region <- factor(d$Region)
  x$fake$data$pods$Region <- d$Region
  m_by <- .flexy_arm_draws(x$fake, "pods", 1L, NULL, "Region")
  ref_by <- .arm_means(x$glmm$fits$pods, d, "trt", "Open", "Cage", "Region")
  expect_equal(m_by$high[1L, ], c(ref_by$high), tolerance = 1e-10)
  expect_equal(m_by$low[1L, ], c(ref_by$low), tolerance = 1e-10)

  bad <- x$latent
  bad["Predictor:1", 1L] <- bad["Predictor:1", 1L] + 0.1
  local_mocked_bindings(.flexy_latent = function(model, n, seed) bad)
  expect_error(
    .flexy_arm_draws(x$fake, "pods", 1L, NULL, NULL),
    "do not rebuild the fit's linear predictor"
  )
})

test_that("posterior draws enter the draws contract, reproducibly", {
  skip_flexy()
  op <- flexy_options()
  on.exit(options(op), add = TRUE)
  est <- fx_estimand("multi_site", unit = c("Site", "Plot"))
  fit <- pollin_fit(est, engine = "flexyBayes")

  set.seed(99)
  before <- .Random.seed
  dr <- pollin_draws(fit, n = 200L, seed = 4L)
  expect_identical(.Random.seed, before)
  expect_s3_class(dr, "pollin_draws")
  expect_identical(dr$method, "posterior")
  expect_identical(dr$engine, "flexyBayes")
  expect_identical(dim(dr$estimate), c(200L, 1L))
  expect_identical(colnames(dr$estimate), "overall")
  expect_identical(pollin_draws(fit, n = 200L, seed = 4L)$estimate, dr$estimate)
  expect_output(print(dr), "joint posterior draws \\(INLA\\)")

  expect_error(
    pollin_draws(fit, n = 10L, method = "simulation"),
    "drawn from its posterior"
  )
  expect_error(
    pollin_draws(fit, n = 10L, method = "bootstrap"),
    "drawn from its posterior"
  )
  g <- suppressWarnings(pollin_fit(est))
  expect_error(
    pollin_draws(g, n = 10L, method = "posterior"),
    "needs a fit with `engine = \"flexyBayes\"`"
  )
  expect_error(
    pollin_effects(g, interval = "posterior"),
    "needs a fit with `engine = \"flexyBayes\"`"
  )
  expect_error(pollin_effects(fit, interval = "delta"), "posterior intervals")
})

test_that("arm means per draw match a row-by-row computation", {
  skip_flexy()
  op <- flexy_options()
  on.exit(options(op), add = TRUE)
  est <- fx_estimand("multi_site", unit = c("Site", "Plot"))
  fit <- pollin_fit(est, engine = "flexyBayes")
  m <- .flexy_arm_draws(fit, "pods", 5L, 8L, NULL)
  latent <- .flexy_latent(fit$fits$pods, 5L, 8L)
  lv <- lapply(
    flexyBayes::fb_model_terms(fit$fits$pods)$random_terms,
    `[[`,
    "var_levels"
  )
  d <- fit$data$pods
  effect <- function(term, level, levels, k) {
    i <- match(level, levels)
    if (is.na(i)) 0 else latent[paste0(term, ":", i), k]
  }
  for (k in c(1L, 5L)) {
    mu <- function(arm) {
      eta <- vapply(
        seq_len(nrow(d)),
        function(i) {
          latent["(Intercept):1", k] +
            (arm == "Open") * latent["trtOpen:1", k] +
            effect("Site", d$Site[i], lv[[1L]], k) +
            effect(
              "pollin_Site_trt",
              paste(d$Site[i], arm, sep = ":"),
              lv[[2L]],
              k
            ) +
            effect(
              "pollin_Site_Plot",
              as.character(d$pollin_Site_Plot[i]),
              lv[[3L]],
              k
            )
        },
        numeric(1)
      )
      mean(exp(eta))
    }
    expect_equal(unname(m$high[k, 1L]), mu("Open"), tolerance = 1e-12)
    expect_equal(unname(m$low[k, 1L]), mu("Cage"), tolerance = 1e-12)
  }
})

test_that("the effects table is the summary of the draws", {
  skip_flexy()
  op <- flexy_options()
  on.exit(options(op), add = TRUE)
  est <- fx_estimand("rcbd")
  fit <- pollin_fit(est, engine = "flexyBayes", by = "Region")
  for (by in list("Region", NULL)) {
    ef <- pollin_effects(fit, by = by, n_draws = 300L, seed = 6L)
    expect_identical(ef$interval, "posterior")
    expect_identical(ef$n_draws, 300L)
    tab <- ef$table
    expect_true(all(is.na(tab$p.value)))
    expect_identical(
      names(tab),
      names(suppressWarnings(pollin_effects(pollin_fit(est, by = by))$table))
    )
    if (!is.null(by)) {
      dr <- pollin_draws(fit, n = 300L, seed = 6L)
      s <- summary(dr)
      expect_identical(tab$Region, s$level)
      expect_identical(tab$estimate, s$median)
      expect_identical(tab$conf.low, s$conf.low)
      expect_identical(tab$conf.high, s$conf.high)
      expect_identical(
        tab$std.error,
        unname(apply(dr$estimate, 2L, stats::sd))
      )
    }
  }
  expect_output(print(ef), "posterior, equal-tailed \\(300 draws\\)")
})

test_that("population, value and the total run on posterior draws", {
  skip_flexy()
  op <- flexy_options()
  on.exit(options(op), add = TRUE)
  pop <- pollin_population("Region", c(north = 3, south = 1))
  est <- fx_estimand("rcbd", population = pop)
  fit <- pollin_fit(est, engine = "flexyBayes", by = "Region")
  dr <- pollin_draws(fit, n = 200L, seed = 2L)
  expect_identical(colnames(dr$estimate), c("north", "south", "population"))
  expect_equal(
    dr$high[, "population"],
    0.75 * dr$high[, "north"] + 0.25 * dr$high[, "south"]
  )
  m <- pollin_market(1000, 500, "t", "AUD", 2024, "Synthetic", crop = "canola")
  v <- pollin_value(dr, m, channel = "dependence")
  expect_s3_class(v, "pollin_value")
  expect_identical(v$summary$method, "posterior")
  expect_equal(
    v$draws,
    (1 - dr$low[, "population"] / dr$high[, "population"]) * 500 * 1000
  )
  m2 <- pollin_market(1, 1, "t", "AUD", 2024, "Synthetic", crop = "faba bean")
  tot <- pollin_guard_doublecount(
    v,
    pollin_value(dr, m2, channel = "dependence")
  )
  expect_s3_class(tot, "pollin_value_total")
})

test_that("the posterior estimate sits near the glmmTMB estimate", {
  skip_flexy()
  op <- flexy_options()
  on.exit(options(op), add = TRUE)
  est <- fx_estimand("rcbd", scale = "lnrr", d = fx_df(seed = 11L, n_site = 8L))
  b <- pollin_effects(
    pollin_fit(est, engine = "flexyBayes"),
    n_draws = 1000L,
    seed = 1L
  )$table
  g <- suppressWarnings(pollin_effects(pollin_fit(est)))$table
  expect_lt(abs(b$estimate - g$estimate), g$std.error)
  expect_gt(b$std.error / g$std.error, 0.7)
  expect_lt(b$std.error / g$std.error, 1.5)
})

test_that("a Poisson fit keeps one linear predictor per data row", {
  skip_flexy()
  op <- flexy_options()
  on.exit(options(op), add = TRUE)
  set.seed(31)
  d <- data.frame(
    trt = rep(rep(c("Cage", "Open"), each = 6L), 8L),
    Block = rep(paste0("b", seq_len(8L)), each = 12L)
  )
  d$visits <- rpois(nrow(d), exp(2 + 0.4 * (d$trt == "Open")))
  pd <- pollin_data(d, "visits", "trt", unit = "Block")
  ds <- pollin_design(pd, c(excluded = "Cage", open = "Open"), "rcbd")
  est <- pollin_estimand(ds, "dependence", "ratio", "loss")
  fit <- pollin_fit(est, engine = "flexyBayes", family = stats::poisson())
  eff <- pollin_effects(fit, n_draws = 200L, seed = 2L)
  expect_identical(nrow(eff$table), 1L)
  expect_true(is.finite(eff$table$estimate))
})
