# Layer L4: pollin_cate() ------------------------------------------------------

cate_slopes <- c(-0.5, -0.3, -0.2, -0.1, 0, 0, 0.1, 0.2, 0.3, 0.5)

# Ten varieties within four sites, each variety with its own treatment slope.
cate_df <- function(seed = 5L) {
  set.seed(seed)
  d <- expand.grid(
    plant = 1:4,
    trt = c("Cage", "Open"),
    Variety = sprintf("v%02d", 1:10),
    Site = sprintf("s%d", 1:4),
    stringsAsFactors = FALSE
  )
  v <- as.integer(factor(d$Variety))
  d$pods <- stats::rnbinom(
    nrow(d),
    mu = exp(
      3 +
        (0.4 + cate_slopes[v]) * (d$trt == "Open") +
        stats::rnorm(4L, 0, 0.2)[as.integer(factor(d$Site))] +
        stats::rnorm(10L, 0, 0.2)[v]
    ),
    size = 20
  )
  d$plant <- NULL
  d
}

test_that("pollin_fit(pooled = ) adds a random treatment slope", {
  est <- fx_estimand("rcbd", d = cate_df())
  inputs <- .pooled_inputs(.estimand_inputs(est), "Variety", NULL)
  expect_identical(inputs$random, c("(1 | Site)", "(trt | Variety)"))
  expect_identical(inputs$pooled, "Variety")
  expect_true(is.factor(inputs$data$Variety))
  expect_error(pollin_fit(est, pooled = "Site"), "declare a slope on a unit")
  expect_error(pollin_fit(est, pooled = "trt"), "other than the treatment")
  expect_error(pollin_fit(est, pooled = "nope"), "one column of the data")
  expect_error(
    pollin_fit(est, pooled = "Variety", by = "Variety"),
    "name the same column"
  )
  skip_on_cran()
  fit <- suppressWarnings(pollin_fit(est, pooled = "Variety"))
  expect_identical(fit$pooled, "Variety")
  expect_identical(fit$rhs, c("trt", "(1 | Site)", "(trt | Variety)"))
  expect_output(print(fit), "pooled    : Variety \\(random treatment slope\\)")
  expect_null(suppressWarnings(pollin_fit(est))$pooled)
})

test_that("pollin_cate resolves the column that carries the slope", {
  skip_on_cran()
  est_r <- fx_estimand("rcbd")
  fit_r <- suppressWarnings(pollin_fit(est_r))
  expect_error(pollin_cate(fit_r), "has no random treatment slope")
  expect_error(
    pollin_cate(fit_r, by = "Region"),
    "carries no random treatment slope.*pooled = \"Region\""
  )
  est_n <- fx_estimand(
    "multi_site",
    unit = c("Site", "Plot"),
    d = fx_df()
  )
  est_n$design <- pollin_design(
    est_n$design$data,
    c(excluded = "Cage", open = "Open"),
    "multi_site",
    slope = "Plot"
  )
  fit_n <- suppressWarnings(pollin_fit(est_n))
  expect_error(pollin_cate(fit_n, by = "Plot"), "nested unit Site:Plot")
  expect_error(pollin_cate(fit_n), "nested unit Site:Plot")
})

test_that("glmmTMB: pooled per-level effects by G-computation", {
  skip_on_cran()
  est <- fx_estimand("rcbd", d = cate_df())
  fit <- suppressWarnings(pollin_fit(est, pooled = "Variety"))
  set.seed(1)
  before <- .Random.seed
  ct <- pollin_cate(fit, n = 500L, seed = 2L)
  expect_identical(.Random.seed, before)
  expect_s3_class(ct, "pollin_cate")
  tab <- ct$table
  expect_identical(tab$level, sprintf("v%02d", 1:10))
  expect_true(all(tab$pooled))
  expect_identical(tab$n, rep(32L, 10L))
  expect_identical(unique(tab$interval), "conditional")
  expect_identical(ct$by, "Variety")
  m <- .arm_means(
    fit$fits$pods,
    fit$data$pods,
    "trt",
    "Open",
    "Cage",
    "Variety"
  )
  expect_equal(tab$estimate, unname(c(1 - m$low / m$high)), tolerance = 1e-12)
  expect_true(all(tab$conf.low < tab$estimate & tab$estimate < tab$conf.high))
  expect_identical(pollin_cate(fit, n = 500L, seed = 2L)$table, tab)
  expect_output(print(ct), "Pooled: partial pooling")
  expect_output(print(ct), "Unpooled: not carried by this fit")
  expect_output(print(ct), "approximate")

  unpooled <- suppressWarnings(pollin_effects(pollin_fit(est, by = "Variety")))
  expect_lt(stats::sd(tab$estimate), stats::sd(unpooled$table$estimate))
})

test_that("glmmTMB: the design's slope unit is read by default", {
  skip_on_cran()
  est <- fx_estimand("multi_site", d = fx_df(seed = 3L, n_site = 8L))
  fit <- suppressWarnings(pollin_fit(est))
  ct <- suppressWarnings(pollin_cate(fit, n = 200L, seed = 1L))
  expect_identical(ct$by, "Site")
  expect_identical(ct$table$level, sprintf("s%02d", 1:8))
})

test_that("flexyBayes: pooled per-level effects from the posterior", {
  skip_flexy()
  op <- flexy_options()
  on.exit(options(op), add = TRUE)
  est <- fx_estimand("rcbd", d = cate_df())
  fit <- pollin_fit(est, engine = "flexyBayes", pooled = "Variety")
  expect_identical(
    fit$rhs,
    c("trt", "Site", "Variety", "pollin_Variety_trt")
  )
  ct <- pollin_cate(fit, n = 300L, seed = 3L)
  tab <- ct$table
  expect_identical(tab$level, sprintf("v%02d", 1:10))
  expect_identical(unique(tab$interval), "posterior")
  m <- .flexy_arm_draws(fit, "pods", 300L, 3L, "Variety")
  reps <- .effect_from_means(m, "ratio", "loss")
  expect_identical(tab$estimate, unname(apply(reps, 2L, stats::median)))
  expect_identical(tab$std.error, unname(apply(reps, 2L, stats::sd)))
  expect_output(print(ct), "posterior, equal-tailed \\(300 draws\\)")

  truth <- 1 - exp(-(0.4 + cate_slopes))
  unpooled <- suppressWarnings(pollin_effects(pollin_fit(est, by = "Variety")))
  expect_lt(stats::sd(tab$estimate), stats::sd(unpooled$table$estimate))
  expect_gt(stats::cor(tab$estimate, truth), 0.9)
})
