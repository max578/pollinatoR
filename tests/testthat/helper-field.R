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
