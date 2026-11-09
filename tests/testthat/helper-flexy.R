# Live flexyBayes fits: skips, options and a small synthetic trial -------------

skip_flexy <- function() {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("flexyBayes")
  testthat::skip_if_not_installed("INLA")
}

# One INLA fit per call instead of the best of three; returns the old options.
flexy_options <- function() {
  options(
    flexyBayes.inla_best_of = 1L,
    flexyBayes.silence_default_prior_note = TRUE
  )
}

fx_df <- function(seed = 7L, n_site = 6L, n_plant = 6L) {
  set.seed(seed)
  d <- expand.grid(
    plant = seq_len(n_plant),
    trt = c("Cage", "Open"),
    Site = sprintf("s%02d", seq_len(n_site)),
    stringsAsFactors = FALSE
  )
  s <- as.integer(factor(d$Site))
  d$Region <- ifelse(s <= n_site / 2L, "north", "south")
  d$Plot <- ifelse(d$plant <= n_plant / 2L, "p1", "p2")
  site <- stats::rnorm(n_site, 0, 0.3)[s]
  slope <- stats::rnorm(n_site, 0, 0.15)[s]
  d$pods <- stats::rnbinom(
    nrow(d),
    mu = exp(3 + (0.4 + slope) * (d$trt == "Open") + site),
    size = 10
  )
  d$plant <- NULL
  d
}

fx_estimand <- function(
  layout = "rcbd",
  unit = "Site",
  scale = "ratio",
  basis = "loss",
  population = "trial",
  d = fx_df()
) {
  pd <- pollin_data(d, "pods", "trt", unit = unit)
  ds <- pollin_design(pd, c(excluded = "Cage", open = "Open"), layout)
  pollin_estimand(ds, "dependence", scale, basis, population = population)
}
