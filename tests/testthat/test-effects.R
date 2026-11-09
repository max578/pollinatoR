# Layer L4: pollin_effects() ---------------------------------------------------

# Balanced: six blocks x both arms x four plants.
l4_df <- function(seed = 5L) {
  set.seed(seed)
  d <- expand.grid(
    plant = 1:4,
    trt = c("Cage", "Open"),
    Block = sprintf("b%d", 1:6),
    stringsAsFactors = FALSE
  )
  blk <- as.integer(factor(d$Block))
  d$Region <- ifelse(blk <= 3L, "north", "south")
  d$pods <- stats::rpois(
    nrow(d),
    exp(3 + 0.4 * (d$trt == "Open") + stats::rnorm(6L, 0, 0.2)[blk])
  )
  d$mass <- 20 +
    5 * (d$trt == "Open") +
    stats::rnorm(6L, 0, 2)[blk] +
    stats::rnorm(nrow(d))
  d$plant <- NULL
  d
}
l4_arms <- c(excluded = "Cage", open = "Open")

l4_fit <- function(
  scale,
  basis = "loss",
  response = "pods",
  family = NULL,
  by = NULL
) {
  pd <- pollin_data(l4_df(), response, "trt", unit = "Block")
  ds <- pollin_design(pd, l4_arms, "rcbd")
  pollin_fit(
    pollin_estimand(ds, "dependence", scale, basis),
    family = family,
    by = by
  )
}

effects_quiet <- function(...) suppressWarnings(pollin_effects(...))
delta_quiet <- function(...) suppressWarnings(.effects_delta(...))

test_that("every scale is labelled with its family, scale and basis", {
  skip_on_cran()
  defs <- list(
    difference = "Open - Cage \\(difference\\)",
    ratio_loss = "1 - Cage/Open \\(loss share\\)",
    ratio_gain = "Open/Cage - 1 \\(gain\\)",
    rii = "\\(Open - Cage\\)/\\(Open \\+ Cage\\) \\(RII\\)",
    lnrr = "log\\(Open/Cage\\) \\(lnRR\\)"
  )
  for (s in c("difference", "ratio", "rii", "lnrr")) {
    for (b in c("loss", "gain")) {
      ef <- effects_quiet(l4_fit(s, b))
      expect_s3_class(ef, "pollin_effects")
      tab <- ef$table
      expect_identical(tab$family, "dependence")
      expect_identical(tab$scale, s)
      expect_identical(tab$basis, b)
      expect_identical(tab$distribution, "nbinom2")
      expect_identical(tab$interval, "delta")
      out <- capture.output(print(ef))
      expect_true(any(grepl(paste0("basis +: ", b), out)))
      key <- if (s == "ratio") paste0("ratio_", b) else s
      expect_true(any(grepl(defs[[key]], out)), info = paste(s, b))
    }
  }
})

test_that("the ratio of averaged predictions equals the ratio of cell means", {
  # Fixed-effects-only Poisson fit on a balanced design: the MLE reproduces
  # the cell means, and var(log r) = 1 / sum(y_open) + 1 / sum(y_cage).
  d <- l4_df()
  d$trt <- stats::relevel(factor(d$trt), ref = "Cage")
  m <- glmmTMB::glmmTMB(pods ~ trt, data = d, family = stats::poisson())
  pd <- pollin_data(l4_df(), "pods", "trt", unit = "Block")
  ds <- pollin_design(pd, l4_arms, "rcbd")
  cell <- tapply(d$pods, d$trt, mean)
  sums <- tapply(d$pods, d$trt, sum)
  r <- cell[["Open"]] / cell[["Cage"]]

  lnrr <- delta_quiet(
    m,
    "trt",
    pollin_estimand(ds, "dependence", "lnrr", "loss"),
    NULL,
    0.95
  )
  expect_equal(lnrr$estimate, log(r), tolerance = 1e-8)
  expect_equal(
    lnrr$std.error,
    sqrt(1 / sums[["Open"]] + 1 / sums[["Cage"]]),
    tolerance = 1e-5
  )
  loss <- delta_quiet(
    m,
    "trt",
    pollin_estimand(ds, "dependence", "ratio", "loss"),
    NULL,
    0.95
  )
  expect_equal(loss$estimate, 1 - 1 / r, tolerance = 1e-8)
  diff <- delta_quiet(
    m,
    "trt",
    pollin_estimand(ds, "dependence", "difference", "loss"),
    NULL,
    0.95
  )
  expect_equal(diff$estimate, cell[["Open"]] - cell[["Cage"]], tolerance = 1e-8)
})

test_that("a balanced Gaussian block design gives the cell-mean ratio", {
  skip_on_cran()
  d <- l4_df()
  cell <- tapply(d$mass, d$trt, mean)
  ef <- effects_quiet(l4_fit(
    "ratio",
    "gain",
    response = "mass",
    family = stats::gaussian()
  ))
  expect_equal(
    ef$table$estimate,
    cell[["Open"]] / cell[["Cage"]] - 1,
    tolerance = 1e-6
  )
})

test_that("the scales are one ratio: RII = D / (2 - D) and its siblings", {
  skip_on_cran()
  fit <- l4_fit("ratio", "loss")
  loss <- effects_quiet(fit)$table
  fit$estimand$scale <- "rii"
  rii <- effects_quiet(fit)$table
  fit$estimand$scale <- "lnrr"
  lnrr <- effects_quiet(fit)$table
  fit$estimand$scale <- "ratio"
  fit$estimand$basis <- "gain"
  gain <- effects_quiet(fit)$table
  for (col in c("estimate", "conf.low", "conf.high")) {
    D <- loss[[col]]
    expect_equal(rii[[col]], D / (2 - D), tolerance = 1e-10)
    expect_equal(lnrr[[col]], -log(1 - D), tolerance = 1e-10)
    expect_equal(gain[[col]], D / (1 - D), tolerance = 1e-10)
  }
  expect_equal(rii$p.value, loss$p.value, tolerance = 1e-12)
  expect_lt(loss$conf.high, 1)
  expect_gt(rii$conf.low, -1)
  expect_lt(rii$conf.high, 1)
})

test_that("bounded scales stay bounded under a wide interval", {
  skip_on_cran()
  fit <- l4_fit("ratio", "loss")
  wide <- effects_quiet(fit, conf_level = 0.999999)$table
  expect_lt(wide$conf.high, 1)
  expect_lt(wide$conf.low, wide$estimate)
})

test_that("by gives one contrast per level from the interaction fit", {
  skip_on_cran()
  fit <- l4_fit("ratio", "loss", by = "Region")
  expect_identical(fit$rhs[1L], "trt * Region")
  ef <- effects_quiet(fit)
  expect_identical(ef$table$Region, c("north", "south"))
  d <- fit$data$pods
  m <- .arm_means(fit$fits$pods, d, "trt", "Open", "Cage", "Region")
  expect_equal(
    ef$table$estimate,
    as.vector(1 - m$low / m$high),
    tolerance = 1e-8
  )
  expect_error(effects_quiet(fit, by = "Zone"), "refit with `pollin_fit")
  expect_error(
    effects_quiet(l4_fit("ratio", "loss"), by = "Region"),
    "refit with `pollin_fit\\(\\.\\.\\., by = \"Region\"\\)`"
  )
  pooled <- effects_quiet(fit, by = NULL)$table
  expect_identical(nrow(pooled), 1L)
  expect_error(l4_fit("ratio", "loss", by = "trt"), "other than the treatment")
  expect_error(l4_fit("ratio", "loss", by = "Zone"), "`by`")
})

test_that("by recovers planted levels that a pooled fit would merge", {
  skip_on_cran()
  set.seed(1)
  d <- expand.grid(
    rep = 1:10,
    trt = c("Cage", "Open"),
    grp = c("a", "b"),
    site = paste0("s", 1:6),
    stringsAsFactors = FALSE
  )
  open_mean <- ifelse(d$grp == "a", 40, 80)
  d$y <- stats::rpois(nrow(d), ifelse(d$trt == "Open", open_mean, 30))
  pd <- pollin_data(d, "y", "trt", unit = "site", groups = "grp")
  ds <- pollin_design(pd, l4_arms, "rcbd")
  est <- pollin_estimand(ds, "dependence", "ratio", "loss")
  fit <- pollin_fit(est, family = stats::poisson(), by = "grp")
  expect_no_warning(ef <- pollin_effects(fit))
  tab <- ef$table
  truth <- c(a = 1 - 30 / 40, b = 1 - 30 / 80)
  expect_identical(tab$grp, c("a", "b"))
  # Poisson cells of 60 plants give loss-share SEs near 0.023 (a) and 0.010
  # (b). Four SEs leave about a 1-in-10,000 chance of a false failure per level.
  expect_true(all(abs(tab$estimate - truth) < 4 * tab$std.error))
  expect_lt(tab$std.error[1L], 0.04)
  expect_lt(tab$std.error[2L], 0.02)
  expect_lt(tab$conf.high[1L], tab$conf.low[2L])
  expect_output(print(ef), "fixed and random effects jointly")
})

test_that("the parametric bootstrap is seeded and reproducible", {
  skip_on_cran()
  fit <- l4_fit("ratio", "loss")
  set.seed(99L)
  before <- stats::runif(1L)
  set.seed(99L)
  expect_warning(
    b1 <- pollin_effects(fit, interval = "bootstrap", n_boot = 20L, seed = 4L),
    "Bootstrap for 'pods': 0 refits failed and 3 gave a warning, of 20"
  )
  expect_identical(stats::runif(1L), before)
  b2 <- effects_quiet(fit, interval = "bootstrap", n_boot = 20L, seed = 4L)
  expect_identical(b1$table, b2$table)
  expect_identical(b1$table$interval, "bootstrap")
  expect_true(is.na(b1$table$p.value))
  expect_lt(b1$table$conf.low, b1$table$estimate)
  expect_lt(b1$table$conf.high, 1)
  delta <- effects_quiet(fit)$table
  expect_equal(b1$table$estimate, delta$estimate, tolerance = 1e-8)
  expect_output(print(b1), "parametric bootstrap \\(20 refits\\)")

  b_by <- effects_quiet(
    l4_fit("ratio", "loss", by = "Region"),
    interval = "bootstrap",
    n_boot = 10L,
    seed = 4L
  )
  expect_identical(b_by$table$Region, c("north", "south"))
})

test_that("a three-arm design reads the supplemented contrast", {
  skip_on_cran()
  d <- l4_df()
  d$trt[seq(2L, nrow(d), by = 4L)] <- "Hand"
  pd <- pollin_data(d, "pods", "trt", unit = "Block")
  ds <- pollin_design(
    pd,
    c(excluded = "Cage", open = "Open", supplemented = "Hand"),
    "rcbd",
    reference = "supplemented"
  )
  ef <- effects_quiet(pollin_fit(pollin_estimand(
    ds,
    "deficit",
    "difference",
    "loss"
  )))
  expect_identical(ef$table$contrast, "Hand vs Open")
  expect_identical(ef$table$family, "deficit")
})

test_that("pollin_effects refuses bad arguments", {
  skip_on_cran()
  fit <- l4_fit("ratio", "loss")
  expect_error(pollin_effects(fit$estimand), "pollin_fit")
  expect_error(pollin_effects(fit, interval = "fieller"), "interval")
  expect_error(pollin_effects(fit, conf_level = 95), "conf_level")
  expect_error(
    pollin_effects(fit, interval = "bootstrap", n_boot = 1L),
    "n_boot"
  )
  expect_error(
    pollin_effects(fit, interval = "bootstrap", seed = "a"),
    "seed"
  )
})
