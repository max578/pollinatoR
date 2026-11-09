# Joint delta interval on the glmmTMB engine -----------------------------------
#
# Oracles fixed before the first run (2026-11-10):
#   * the arm design matrices times the fitted (beta, b), plus the offset,
#     equal predict(type = "link") on the arm-set data, to 1e-10;
#   * the delta estimate equals the G-computation of .arm_means(), to 1e-10;
#   * the analytic gradient equals a central finite difference, to 1e-6;
#   * with the variances held fixed (map), the joint covariance of a Gaussian
#     mixed model equals the inverse of Henderson's mixed-model equations
#     computed by hand, to 1e-6 relative;
#   * interval = "delta_fixed" equals marginaleffects::avg_comparisons().

jt_df <- function(seed = 11L, n_site = 10L) {
  set.seed(seed)
  d <- expand.grid(
    rep = 1:6,
    trt = c("Cage", "Open"),
    Site = sprintf("s%02d", seq_len(n_site)),
    stringsAsFactors = FALSE
  )
  s <- as.integer(factor(d$Site))
  u0 <- stats::rnorm(n_site, 0, 0.3)
  u1 <- stats::rnorm(n_site, 0, 0.25)
  open <- d$trt == "Open"
  d$Region <- ifelse(s <= n_site %/% 2L, "north", "south")
  d$pods <- stats::rpois(nrow(d), exp(3 + u0[s] + (0.3 + u1[s]) * open))
  d$ovules <- 30L
  d$seeds <- stats::rbinom(
    nrow(d),
    d$ovules,
    stats::plogis(0.2 + u0[s] + 0.6 * open)
  )
  d$mass <- 5 + u0[s] + 0.8 * open + stats::rnorm(nrow(d), 0, 0.5)
  d
}
jt_arms <- c(excluded = "Cage", open = "Open")

jt_fit <- function(
  response = "pods",
  scale = "ratio",
  slope = NULL,
  by = NULL,
  denominator = NULL,
  ...
) {
  pd <- pollin_data(
    jt_df(),
    response,
    "trt",
    unit = "Site",
    groups = "Region",
    denominator = denominator
  )
  ds <- pollin_design(pd, jt_arms, "rcbd", slope = slope)
  est <- pollin_estimand(ds, "dependence", scale, "loss")
  pollin_fit(est, by = by, ...)
}

jt_cases <- function() {
  list(
    slope = jt_fit(slope = "Site", family = stats::poisson()),
    by = jt_fit(by = "Region", family = stats::poisson()),
    binomial = jt_fit("seeds", denominator = "ovules")
  )
}

test_that("the arm designs reproduce predict() on the link scale", {
  skip_on_cran()
  for (fit in jt_cases()) {
    r <- names(fit$fits)
    model <- fit$fits[[r]]
    d <- fit$data[[r]]
    jc <- .joint_cov(model, r)
    designs <- .joint_arm_designs(fit, r)
    beta <- jc$par[seq_len(jc$n_beta)]
    b <- jc$par[-seq_len(jc$n_beta)]
    for (arm in c("high", "low")) {
      a <- designs[[arm]]
      nd <- d
      nd$trt <- factor(
        rep(fit$estimand$contrast[[arm]], nrow(d)),
        levels = levels(d$trt)
      )
      eta <- as.vector(a$x %*% beta + a$z %*% b) + a$offset
      want <- stats::predict(model, newdata = nd, type = "link")
      expect_equal(eta, unname(want), tolerance = 1e-10)
    }
  }
})

test_that("the joint delta estimate is the G-computation of .arm_means()", {
  skip_on_cran()
  for (fit in jt_cases()) {
    r <- names(fit$fits)
    ef <- pollin_effects(fit)$table
    e <- fit$estimand
    am <- .arm_means(
      fit$fits[[r]],
      fit$data[[r]],
      "trt",
      "Open",
      "Cage",
      fit$by
    )
    want <- .effect_from_means(am, e$scale, e$basis)
    expect_equal(ef$estimate, unname(as.vector(want)), tolerance = 1e-10)
    expect_identical(unique(ef$interval), "delta")
  }
})

test_that("the analytic gradient matches a finite difference", {
  skip_on_cran()
  fit <- jt_cases()$slope
  model <- fit$fits$pods
  jc <- .joint_cov(model, "pods")
  designs <- .joint_arm_designs(fit, "pods")
  g <- rep("all", nrow(fit$data$pods))
  lnrr <- function(par) {
    h <- .joint_arm_gradient(model, designs$high, par, jc$n_beta, g)
    l <- .joint_arm_gradient(model, designs$low, par, jc$n_beta, g)
    log(h$mean) - log(l$mean)
  }
  h <- .joint_arm_gradient(model, designs$high, jc$par, jc$n_beta, g)
  l <- .joint_arm_gradient(model, designs$low, jc$par, jc$n_beta, g)
  analytic <- as.vector(h$grad / h$mean - l$grad / l$mean)
  step <- 1e-5
  numeric <- vapply(
    seq_along(jc$par),
    function(i) {
      e <- replace(numeric(length(jc$par)), i, step)
      (lnrr(jc$par + e) - lnrr(jc$par - e)) / (2 * step)
    },
    numeric(1L)
  )
  expect_equal(analytic, numeric, tolerance = 1e-6)
})

test_that("with fixed variances the joint covariance is Henderson's inverse", {
  skip_on_cran()
  d <- jt_df()
  sd_site <- 0.4
  sd_res <- 0.5
  model <- glmmTMB::glmmTMB(
    mass ~ trt + (1 | Site),
    data = d,
    family = stats::gaussian(),
    start = list(theta = log(sd_site), betadisp = log(sd_res)),
    map = list(theta = factor(NA), betadisp = factor(NA))
  )
  jc <- .joint_cov(model, "mass")
  x <- stats::model.matrix(~trt, d)
  z <- stats::model.matrix(~ 0 + Site, d)
  mme <- rbind(
    cbind(crossprod(x), crossprod(x, z)),
    cbind(crossprod(z, x), crossprod(z) + diag(sd_res^2 / sd_site^2, ncol(z)))
  )
  want <- solve(mme) * sd_res^2
  expect_equal(unname(jc$v), unname(want), tolerance = 1e-6)
})

test_that("delta_fixed is the marginaleffects delta interval", {
  skip_on_cran()
  fit <- jt_cases()$slope
  ef <- pollin_effects(fit, interval = "delta_fixed")$table
  cmp <- suppressWarnings(
    marginaleffects::avg_comparisons(
      fit$fits$pods,
      variables = list(trt = c("Cage", "Open")),
      comparison = "lnratioavg"
    )
  )
  expect_equal(ef$estimate, 1 - exp(-cmp$estimate), tolerance = 1e-10)
  expect_equal(ef$conf.low, 1 - exp(-cmp$conf.low), tolerance = 1e-10)
  expect_identical(unique(ef$interval), "delta_fixed")
})

test_that("the joint interval is narrower than delta_fixed on a slope design", {
  skip_on_cran()
  fit <- jt_cases()$slope
  joint <- pollin_effects(fit)$table
  fixed <- pollin_effects(fit, interval = "delta_fixed")$table
  expect_lt(joint$std.error, fixed$std.error)
  expect_lt(joint$conf.high - joint$conf.low, fixed$conf.high - fixed$conf.low)
})

test_that("a zero-inflated fit is refused with a pointer to the bootstrap", {
  skip_on_cran()
  fit <- jt_fit(family = stats::poisson(), ziformula = ~1)
  expect_error(pollin_effects(fit), "zero-inflated.*bootstrap")
  expect_error(pollin_draws(fit, n = 10L, seed = 1L), "zero-inflated")
})
