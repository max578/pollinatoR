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

test_that("a direct effect adds its adjustment as main effects", {
  skip_on_cran()
  tot <- pollin_fit(fd_estimand())
  dir <- pollin_fit(fd_estimand("direct"))
  expect_identical(tot$rhs, c("trt", "(1 | Site)"))
  expect_identical(dir$rhs, c("trt", "veg", "(1 | Site)"))
  # On the log link with main effects only, the ratio of the averaged arm
  # predictions is exp(b_trt) exactly, so the loss share is 1 - exp(-b_trt)
  b <- glmmTMB::fixef(dir$fits$pods)$cond[["trtOpen"]]
  e_dir <- pollin_effects(dir)
  expect_equal(e_dir$table$estimate, 1 - exp(-b), tolerance = 1e-8)
  expect_identical(e_dir$table$effect, "direct")
  expect_identical(pollin_effects(tot)$table$effect, "total")
  # The cage lowers veg, which lowers pods, so the total exceeds the direct
  expect_gt(pollin_effects(tot)$table$estimate, e_dir$table$estimate)
  expect_output(print(e_dir), "effect   : direct, holding veg fixed")
})

test_that("a direct effect drops rows missing an adjustment column", {
  skip_on_cran()
  est <- fd_estimand("direct")
  est$design$data$data$veg[1:3] <- NA
  expect_identical(nrow(pollin_fit(est)$data$pods), 157L)
  est_t <- fd_estimand()
  est_t$design$data$data$veg[1:3] <- NA
  expect_identical(nrow(pollin_fit(est_t)$data$pods), 160L)
})

test_that("a site-level exposure is fitted between units only", {
  skip_on_cran()
  ds <- fd_moderation_design(3L, site_only = TRUE)
  expect_message(
    m <- pollin_moderation(ds, "y"),
    "`bees` does not vary within units"
  )
  expect_identical(m$table$component, "between (vulnerable)")
  expect_false(anyNA(m$table$estimate))
  expect_identical(m$between_only, "bees")
  a <- m$assumptions
  expect_match(
    a$statement[a$assumption == "within-unit moderation"],
    "Not estimated for bees"
  )
  expect_output(print(m), "between units only: bees")
})

test_that("the excluded-arm slope is reported as the negative control", {
  skip_on_cran()
  # Fertility raises both arms alike and the bee count: no true moderation
  m <- pollin_moderation(fd_moderation_design(2L, conf = 0.5), "y")
  nc <- m$negative_control
  expect_identical(nc$term, c("bees_within", "bees_between"))
  co <- summary(m$fit)$coefficients$cond
  expect_identical(nc$estimate, unname(co[nc$term, 1L]))
  expect_lt(nc$p.value[nc$term == "bees_between"], 0.001)
  expect_output(print(m), "negative control")
  # Without an excluded arm there is no negative control
  pd <- m$design$data
  pd$data$trt[pd$data$trt == "Cage"] <- "Hand"
  ds2 <- pollin_design(pd, c(open = "Open", supplemented = "Hand"),
    "multi_site")
  m2 <- pollin_moderation(ds2, "y")
  expect_null(m2$negative_control)
  expect_match(
    m2$assumptions$statement[m2$assumptions$assumption == "negative control"],
    "no excluded arm"
  )
})

test_that("a graph's exposure set enters as main effects and interactions", {
  skip_on_cran()
  skip_if_not_installed("dagitty")
  ds <- fd_moderation_design(4L)
  g <- pollin_dag(
    ds,
    "dag { forest -> bees; forest -> y; bees -> y; trt -> y }",
    "y"
  )
  m <- pollin_moderation(ds, "y", dag = g)
  expect_identical(m$adjustment, "forest")
  expect_true("trtOpen:forest" %in% m$table$term)
  expect_identical(
    m$table$component[m$table$term == "trtOpen:forest"],
    "adjustment"
  )
  expect_match(
    m$assumptions$statement[m$assumptions$assumption == "adjustment"],
    "cannot draw effect modification"
  )
  # Without the graph the declared adjust column is a main effect only
  ds_a <- fd_moderation_design(4L, adjust = "forest")
  m_a <- pollin_moderation(ds_a, "y")
  expect_identical(m_a$adjustment, "forest")
  expect_false(any(m_a$table$component == "adjustment"))
  # A latent common cause leaves no set
  g_lat <- pollin_dag(
    ds,
    "dag { u [latent]; u -> bees; u -> y; bees -> y; trt -> y }",
    "y"
  )
  expect_error(
    pollin_moderation(ds, "y", dag = g_lat),
    "identifies no adjustment set for the effect of `bees`"
  )
  expect_error(pollin_moderation(ds_a, "y", dag = g), "different design")
})

test_that("a fit without a positive-definite Hessian warns by name", {
  skip_on_cran()
  ds <- fd_moderation_design(1L, site_sd = 0, slope_sd = 0)
  seen <- character(0)
  withCallingHandlers(
    suppressMessages(pollin_moderation(ds, "y")),
    warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("moderation fit for `y` .*not reliable", seen)))
})
