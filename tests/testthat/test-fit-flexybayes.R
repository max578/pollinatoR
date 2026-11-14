# Layer L3: the flexyBayes engine ----------------------------------------------

test_that("glmmTMB families map to flexyBayes families by name", {
  expect_identical(
    .flexy_family(glmmTMB::nbinom2(), "y"),
    list(family = "negative_binomial", link = "log")
  )
  expect_identical(.flexy_family(stats::poisson(), "y")$family, "poisson")
  expect_identical(.flexy_family(stats::binomial(), "y")$family, "binomial")
  expect_identical(
    .flexy_family(stats::Gamma(link = "log"), "y")$family,
    "gamma"
  )
  expect_identical(.flexy_family(glmmTMB::beta_family(), "y")$family, "beta")
  expect_error(
    .flexy_family(glmmTMB::tweedie(), "y"),
    "Family 'tweedie' of 'y' has no flexyBayes counterpart"
  )
  expect_error(.flexy_family(stats::Gamma(), "y"), "link 'inverse'")
})

test_that("the engine route refuses brms and what it sets itself", {
  expect_error(
    .flexy_check_dots(list(backend = "brms")),
    "posterior draws of a brms fit are not in the flexyBayes extension API"
  )
  expect_error(.flexy_check_dots(list(backend = "auto")), "INLA|inla")
  expect_error(.flexy_check_dots(list(data = 1)), "`data` is set by")
  expect_true(.flexy_check_dots(list(backend = "inla", seed = 1L)))
})

test_that("random terms mirror the design, with the slope as a cell term", {
  d <- fx_df()
  r <- .flexy_random(d, c("Site", "Plot"), "trt", "Site")
  expect_identical(r$terms, c("Site", "pollin_Site_trt", "pollin_Site_Plot"))
  expect_identical(r$cells, list(list(term = "pollin_Site_trt", key = "Site")))
  expect_identical(
    as.character(r$data$pollin_Site_trt[1L]),
    paste(d$Site[1L], d$trt[1L], sep = ":")
  )
  expect_identical(.flexy_random(d, "Site", "trt", NULL)$terms, "Site")
  d$pollin_Site_Plot <- 1
  expect_error(
    .flexy_random(d, c("Site", "Plot"), "trt", NULL),
    "already hold a column 'pollin_Site_Plot'"
  )
})

test_that("R/ calls flexyBayes through its extension API only (C5)", {
  allow <- c(
    "flexybayes",
    "fb",
    "fb_extension_api",
    "fb_model_terms",
    "fb_response_values",
    "fb_offset_values",
    "fb_model_weights",
    "fb_known_variance",
    "fb_engine_settings",
    "fb_convergence",
    "fb_inla_posterior_sample",
    "fb_check_refit",
    "fb_fixed_design",
    "fb_refit",
    "fb_marginal_cov",
    "fb_isolated_lapply",
    "fb_register_refusals",
    "fb_refuse",
    "fb_refusals"
  )
  dir <- test_path("..", "..", "R")
  skip_if_not(dir.exists(dir), "package sources not available")
  src <- unlist(lapply(
    list.files(dir, "\\.R$", full.names = TRUE),
    readLines,
    warn = FALSE
  ))
  expect_false(any(grepl("flexyBayes:::", src, fixed = TRUE)))
  used <- regmatches(src, gregexpr("flexyBayes::[A-Za-z_.][A-Za-z0-9_.]*", src))
  used <- unique(sub("flexyBayes::", "", unlist(used), fixed = TRUE))
  expect_gt(length(used), 0L)
  expect_identical(setdiff(used, allow), character(0))
})

test_that("a flexyBayes fit carries the engine and the design's terms", {
  skip_flexy()
  op <- flexy_options()
  on.exit(options(op), add = TRUE)
  est <- fx_estimand("multi_site", unit = c("Site", "Plot"))
  fit <- pollin_fit(est, engine = "flexyBayes", by = "Region")
  expect_s3_class(fit, "pollin_fit")
  expect_identical(fit$engine, "flexyBayes")
  expect_identical(
    fit$rhs,
    c("trt * Region", "Site", "pollin_Site_trt", "pollin_Site_Plot")
  )
  expect_s3_class(fit$fits$pods, "flexybayes")
  expect_output(print(fit), "engine    : flexyBayes")
  expect_error(
    pollin_fit(est, engine = "flexyBayes", backend = "brms"),
    "extension API version 1"
  )
  expect_error(
    pollin_fit(est, engine = "flexyBayes", family = glmmTMB::tweedie()),
    "no flexyBayes counterpart"
  )
})

test_that("the API version and the convergence verdict are enforced", {
  skip_flexy()
  op <- flexy_options()
  on.exit(options(op), add = TRUE)
  est <- fx_estimand()
  local_mocked_bindings(
    fb_extension_api = function() 2L,
    .package = "flexyBayes"
  )
  expect_error(
    pollin_fit(est, engine = "flexyBayes"),
    "extension API version 1; the installed flexyBayes offers version 2"
  )
})

test_that("a FALSE convergence verdict stops and NA warns", {
  skip_flexy()
  op <- flexy_options()
  on.exit(options(op), add = TRUE)
  est <- fx_estimand()
  local_mocked_bindings(
    fb_convergence = function(fit) list(ok = FALSE, reasons = "mode not found"),
    .package = "flexyBayes"
  )
  expect_error(
    pollin_fit(est, engine = "flexyBayes"),
    "did not converge: mode not found"
  )
  local_mocked_bindings(
    fb_convergence = function(fit) list(ok = NA, reasons = character(0)),
    .package = "flexyBayes"
  )
  expect_warning(
    pollin_fit(est, engine = "flexyBayes"),
    "no convergence verdict for 'pods'"
  )
})

test_that("a direct effect carries its adjustment on the flexyBayes engine", {
  skip_flexy()
  op <- flexy_options()
  on.exit(options(op), add = TRUE)
  fit <- pollin_fit(fd_estimand("direct"), engine = "flexyBayes")
  expect_identical(fit$flexy$fixed, c("trt", "veg"))
  e <- pollin_effects(fit, n_draws = 1000L, seed = 1L)
  expect_identical(e$table$effect, "direct")
  ref <- pollin_effects(pollin_fit(fd_estimand("direct")))$table$estimate
  expect_lt(abs(e$table$estimate - ref), 0.02)
})
