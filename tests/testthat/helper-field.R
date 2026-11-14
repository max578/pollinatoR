# Field-trial fixtures shared across layers ----------------------------------

# A field trial where the cage lowers vegetative weight and both drive pods
fd_estimand <- function(effect = "total") {
  set.seed(11)
  d <- expand.grid(
    plant = seq_len(8L),
    trt = c("Cage", "Open"),
    Site = sprintf("s%02d", seq_len(10L)),
    stringsAsFactors = FALSE
  )
  open <- d$trt == "Open"
  site <- stats::rnorm(10L, 0, 0.2)[as.integer(factor(d$Site))]
  d$veg <- stats::rgamma(nrow(d), 10, 10 / (40 + 10 * open))
  d$pods <- stats::rpois(nrow(d), exp(2 + site + 0.3 * open + 0.02 * d$veg))
  pd <- pollin_data(d, "pods", "trt", unit = "Site", adjust = "veg")
  ds <- pollin_design(pd, c(excluded = "Cage", open = "Open"), "rcbd")
  pollin_estimand(
    ds,
    "dependence",
    "ratio",
    "loss",
    effect = effect,
    mediators = if (effect == "direct") "veg"
  )
}

# A multi-site study with one exposure. `conf` lets site fertility raise both
# arms alike and the bee count; `site_only` makes the count constant in a site.
fd_moderation_design <- function(
  seed,
  site_sd = 0.3,
  slope_sd = 0.2,
  conf = 0,
  site_only = FALSE,
  adjust = NULL
) {
  set.seed(seed)
  n_site <- 30L
  s <- rep(seq_len(n_site), each = 12L)
  trt <- rep(rep(c("Cage", "Open"), each = 6L), n_site)
  fert <- stats::rnorm(n_site)
  bees_site <- stats::rpois(n_site, exp(2 + conf * fert))
  bees <- bees_site[s]
  if (!site_only) {
    bees <- bees + stats::rpois(length(s), 2)
  }
  shift <- stats::rnorm(n_site, 0, site_sd)[s] + 0.8 * conf * fert[s]
  slope <- 0.3 + stats::rnorm(n_site, 0, slope_sd)[s]
  d <- data.frame(
    y = stats::rpois(length(s), exp(2 + shift + slope * (trt == "Open"))),
    trt = trt,
    Site = sprintf("s%02d", s),
    bees = bees,
    forest = round(fert + stats::rnorm(n_site, 0, 0.5), 3)[s]
  )
  pd <- pollin_data(d, "y", "trt", unit = "Site", exposures = "bees",
    adjust = adjust)
  pollin_design(pd, c(excluded = "Cage", open = "Open"), "multi_site")
}
