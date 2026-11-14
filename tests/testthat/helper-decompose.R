# Fixtures for the whole-yield decomposition ----------------------------------

# Yield is the product of three components, each with its own site effect and
# plant noise, independent within an arm, so the planted log ratio of yield is
# the sum of the planted component effects `b` and the interaction term is 0.
decompose_fixture <- function(
  n_sites,
  n_plants,
  seed,
  b = c(pods = 0.10, seeds = 0.30, weight = -0.10)
) {
  set.seed(seed)
  d <- data.frame(
    trt = rep(rep(c("Cage", "Open"), each = n_plants), n_sites),
    Site = rep(sprintf("s%02d", seq_len(n_sites)), each = 2L * n_plants)
  )
  open <- d$trt == "Open"
  site <- as.integer(factor(d$Site))
  u <- function() stats::rnorm(n_sites, 0, 0.2)[site]
  d$pods <- stats::rnbinom(
    nrow(d),
    mu = exp(4 + u() + b[["pods"]] * open),
    size = 20
  )
  d$seeds <- stats::rgamma(
    nrow(d),
    shape = 20,
    rate = 20 / exp(2.5 + u() + b[["seeds"]] * open)
  )
  d$weight <- stats::rgamma(
    nrow(d),
    shape = 30,
    rate = 30 / exp(1 + u() + b[["weight"]] * open)
  )
  d$yield <- d$pods * d$seeds * d$weight / 1000
  d
}

decompose_fit <- function(d) {
  pd <- pollin_data(d, c("yield", "pods", "seeds", "weight"), "trt",
    unit = "Site")
  ds <- pollin_design(pd, c(excluded = "Cage", open = "Open"), "rcbd")
  est <- pollin_estimand(ds, "dependence", "ratio", "loss")
  pollin_fit(est, family = list(
    yield = stats::Gamma(link = "log"),
    pods = glmmTMB::nbinom2(),
    seeds = stats::Gamma(link = "log"),
    weight = stats::Gamma(link = "log")
  ))
}
