# Layer L3: pollin_fit() and the ported model verbs ---------------------------

# Eight sites x two plots x both arms, five plants each; Poisson pods.
l3_df <- function(seed = 3L) {
  set.seed(seed)
  d <- expand.grid(
    plant = 1:5,
    trt = c("Cage", "Open"),
    Plot = c("p1", "p2"),
    Site = sprintf("s%02d", 1:8),
    stringsAsFactors = FALSE
  )
  site <- as.integer(factor(d$Site))
  d$Region <- ifelse(site <= 4L, "north", "south")
  d$bees <- stats::rpois(8L, 8)[site] + stats::rpois(nrow(d), 1)
  d$pods <- stats::rpois(
    nrow(d),
    exp(3 + 0.3 * (d$trt == "Open") + stats::rnorm(8L, 0, 0.2)[site])
  )
  d$seed_wt <- stats::rgamma(nrow(d), shape = 20, rate = 20 / 2)
  d$plant <- NULL
  d
}
l3_arms <- c(excluded = "Cage", open = "Open")

l3_estimand <- function(
  d = l3_df(),
  scale = "difference",
  basis = "loss",
  layout = "multi_site"
) {
  pd <- pollin_data(
    d,
    c("pods", "seed_wt"),
    "trt",
    unit = "Site",
    groups = "Region",
    exposures = "bees"
  )
  ds <- pollin_design(pd, l3_arms, layout)
  pollin_estimand(ds, "dependence", scale, basis)
}

test_that("pollin_fit fits each response with the low arm as reference", {
  skip_on_cran()
  fit <- suppressWarnings(pollin_fit(l3_estimand()))
  expect_s3_class(fit, "pollin_fit")
  expect_identical(names(fit$fits), c("pods", "seed_wt"))
  expect_identical(levels(fit$data$pods$trt), c("Cage", "Open"))
  expect_identical(fit$rhs, c("trt", "(trt | Site)"))
  expect_identical(fit$families$pods$family, "nbinom2")
  expect_identical(fit$families$seed_wt$family, "Gamma")
  expect_output(print(fit), "<pollin_fit>")

  d <- l3_df()
  d$trt[d$trt == "Open" & d$Plot == "p2"] <- "Hand"
  pd <- pollin_data(d, "pods", "trt", unit = "Site")
  ds <- pollin_design(
    pd,
    c(excluded = "Cage", open = "Open", supplemented = "Hand"),
    "rcbd",
    reference = "supplemented"
  )
  fit3 <- pollin_fit(pollin_estimand(ds, "deficit", "ratio", "loss"))
  expect_identical(levels(fit3$data$pods$trt), c("Open", "Cage", "Hand"))
})

test_that("family takes one family or a list named by response", {
  skip_on_cran()
  est <- l3_estimand()
  one <- suppressWarnings(pollin_fit(est, family = stats::poisson()))
  expect_identical(one$families$seed_wt$family, "poisson")
  lst <- suppressWarnings(pollin_fit(
    est,
    family = list(pods = stats::poisson())
  ))
  expect_identical(lst$families$pods$family, "poisson")
  expect_identical(lst$families$seed_wt$family, "Gamma")
  expect_error(pollin_fit(est, family = list(stats::poisson())), "named")
  expect_error(
    pollin_fit(est, family = list(yield = stats::poisson())),
    "not in the data: yield"
  )
  expect_error(pollin_fit(est, family = "poisson"), "family")
})

test_that("pollin_fit refuses other engines and other objects", {
  est <- l3_estimand()
  expect_error(pollin_fit(est, engine = "brms"), "\"glmmTMB\", \"flexyBayes\"")
  expect_error(pollin_fit(est$design), "pollin_estimand")
})

test_that("a response that fails to fit is dropped with a warning", {
  skip_on_cran()
  d <- l3_df()
  d$bad <- NA_real_
  d$bad[1:3] <- 1
  pd <- pollin_data(d, c("pods", "bad"), "trt", unit = "Site")
  ds <- pollin_design(pd, l3_arms, "rcbd", slope = FALSE)
  est <- pollin_estimand(ds, "dependence", "difference", "loss")
  expect_warning(fit <- pollin_fit(est), "'bad' failed to fit")
  expect_identical(names(fit$fits), "pods")
})

test_that("pollin_contribution reports the difference-scale dependence", {
  skip_on_cran()
  d <- l3_df()
  new <- suppressWarnings(pollin_contribution(
    l3_estimand(d),
    outcomes = "pods"
  ))
  expect_s3_class(new, "pollin_contribution")
  tab <- new$table
  eff <- suppressWarnings(pollin_effects(pollin_fit(l3_estimand(d))))$table
  eff <- eff[eff$response == "pods", ]
  expect_equal(tab$estimate, eff$estimate, tolerance = 1e-10)
  expect_equal(tab$std.error, eff$std.error, tolerance = 1e-10)
  expect_true(is.finite(tab$slope_sd))
  expect_identical(tab$scale, "difference")
  expect_identical(tab$contrast, "Open vs Cage")
  expect_true("p.fdr" %in% names(tab))

  by_tab <- suppressWarnings(pollin_contribution(
    l3_estimand(d),
    outcomes = "pods",
    by = "Region"
  ))$table
  by_fit <- suppressWarnings(pollin_fit(l3_estimand(d), by = "Region"))
  by_eff <- suppressWarnings(pollin_effects(by_fit))$table
  by_eff <- by_eff[by_eff$response == "pods", ]
  expect_equal(by_tab$estimate, by_eff$estimate, tolerance = 1e-10)
  expect_identical(by_tab$Region, by_eff$Region)

  expect_error(
    pollin_contribution(l3_estimand(d, scale = "ratio")),
    "difference scale"
  )
  expect_error(pollin_contribution(l3_estimand(d), outcomes = "x"), "x")
})

test_that("pollin_moderation from a design splits within and between", {
  skip_on_cran()
  d <- l3_df()
  ds <- l3_estimand(d)$design
  new <- suppressWarnings(pollin_moderation(ds, outcome = "pods"))
  expect_s3_class(new, "pollin_moderation")
  expect_setequal(
    new$table$component,
    c("within (robust)", "between (vulnerable)")
  )
  expect_identical(new$design, ds)
  expect_error(pollin_moderation(ds$data, outcome = "pods"), "pollin_design")
  expect_error(pollin_moderation(ds, outcome = "yield"), "yield")
  pd <- pollin_data(d, "pods", "trt", unit = "Site")
  expect_error(
    pollin_moderation(pollin_design(pd, l3_arms, "rcbd"), outcome = "pods"),
    "No exposures"
  )
})

test_that("moderation splits on the full key of an inner slope unit", {
  d <- l3_df()
  pd <- pollin_data(d, "pods", "trt", unit = c("Site", "Plot"))
  ds <- pollin_design(pd, l3_arms, "split_plot", slope = "Plot")
  inputs <- .design_inputs(ds, "Cage")
  expect_identical(nlevels(inputs$wb_group), 16L)
  expect_identical(inputs$slope_key, "Site:Plot")
})
