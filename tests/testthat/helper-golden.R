# Fixed inputs for the golden-value tests (outputs at 0.0.0.9000) -------------

golden_inputs <- function() {
  set.seed(20261108L)
  n_site <- 12L
  n_plot <- 2L
  n_plant <- 6L
  grid <- expand.grid(
    plant = seq_len(n_plant),
    trt = c("Cage", "Open"),
    Plot = paste0("p", seq_len(n_plot)),
    Site = sprintf("s%02d", seq_len(n_site)),
    stringsAsFactors = FALSE
  )
  site_id <- as.integer(factor(grid$Site))
  grid$Region <- ifelse(site_id <= n_site / 2L, "north", "south")
  site_eff <- stats::rnorm(n_site, 0, 0.25)[site_id]
  site_slope <- stats::rnorm(n_site, 0, 0.10)[site_id]
  bees_site <- stats::rpois(n_site, 8)[site_id]
  grid$bees <- bees_site + stats::rpois(nrow(grid), 1)
  grid$wild <- stats::rpois(n_site, 4)[site_id] + stats::rpois(nrow(grid), 1)
  grid$rain <- stats::rnorm(n_site, 300, 40)[site_id]
  open <- grid$trt == "Open"
  mu_pods <- exp(
    3 + (0.30 + site_slope + 0.02 * grid$bees) * open + site_eff
  )
  grid$pods <- stats::rpois(nrow(grid), mu_pods)
  grid$seed_wt <- stats::rgamma(
    nrow(grid),
    shape = 20,
    rate = 20 / exp(1 + 0.12 * open + site_eff / 2)
  )
  grid$seeds_per_pod <- stats::rpois(nrow(grid), exp(2.8 - 0.002 * grid$pods))
  grid$plant <- NULL

  units <- stats::aggregate(
    cbind(pods, seed_wt) ~ Site + trt,
    data = grid,
    FUN = mean
  )
  open_u <- units$pods[units$trt == "Open"]
  cage_u <- units$pods[units$trt == "Cage"]
  vigour <- stats::rnorm(length(open_u), 50, 5)

  list(
    plants = grid,
    open = open_u,
    cage = cage_u,
    vigour = vigour,
    components = grid[, c("pods", "seeds_per_pod", "seed_wt")],
    p = c(0.001, 0.008, 0.02, 0.04, 0.3, 0.8)
  )
}

# Drop heavy or environment-bound members (fits, data copies) so the comparison
# is on the numbers an analyst reads.
golden_strip <- function(x) {
  x <- unclass(x)
  x[setdiff(names(x), c("fits", "fit", "spec", "data"))]
}

golden_path <- function() {
  testthat::test_path("_golden", "golden_v0.0.0.9000.rds")
}

# Fixed inputs for the golden-value tests (outputs at 0.1.0) ------------------

# A three-arm randomised-block trial with two regions: eight blocks x three
# arms x four plants, a count response and a positive continuous response.
golden_inputs_v010 <- function() {
  set.seed(20261109L)
  n_block <- 8L
  d <- expand.grid(
    plant = seq_len(4L),
    trt = c("Cage", "Open", "Hand"),
    Block = sprintf("b%02d", seq_len(n_block)),
    stringsAsFactors = FALSE
  )
  blk <- as.integer(factor(d$Block))
  d$Region <- ifelse(blk <= n_block / 2L, "north", "south")
  blk_eff <- stats::rnorm(n_block, 0, 0.3)[blk]
  beta <- c(Cage = 0, Open = 0.4, Hand = 0.55)[d$trt]
  d$pods <- stats::rnbinom(nrow(d), mu = exp(3 + beta + blk_eff), size = 8)
  d$seed_wt <- stats::rgamma(
    nrow(d),
    shape = 15,
    rate = 15 / exp(1 + beta / 3 + blk_eff / 2)
  )
  d$bees <- stats::rpois(n_block, 8)[blk] + stats::rpois(nrow(d), 1)
  d$plant <- NULL
  d
}

golden_arms_v010 <- c(excluded = "Cage", open = "Open", supplemented = "Hand")

# The scale-by-basis grid of the 0.1.0 effects capture.
golden_grid_v010 <- function() {
  data.frame(
    scale = c("difference", "ratio", "ratio", "rii", "lnrr"),
    basis = c("loss", "loss", "gain", "loss", "loss"),
    stringsAsFactors = FALSE
  )
}

# The numbers an analyst reads from a pollin_effects object.
golden_effects <- function(x) {
  x <- unclass(x)
  x[c("table", "interval", "n_boot", "seed", "conf_level", "by")]
}

golden_path_v010 <- function() {
  testthat::test_path("_golden", "golden_v0.1.0.rds")
}
