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

# The delta route of 0.1.0 to 0.3.1 is `interval = "delta_fixed"` since 0.4.0;
# its outputs are compared to the 0.1.0 oracle under the old label.
golden_fixed_as_delta <- function(x) {
  relabel <- function(v) replace(v, v == "delta_fixed", "delta")
  if (is.data.frame(x)) {
    x$interval <- relabel(x$interval)
    return(x)
  }
  x$table$interval <- relabel(x$table$interval)
  x$interval <- relabel(x$interval)
  x
}

golden_path_v010 <- function() {
  testthat::test_path("_golden", "golden_v0.1.0.rds")
}

# Outputs at 0.2.0: draws, population, value and total ------------------------

golden_path_v020 <- function() {
  testthat::test_path("_golden", "golden_v0.2.0.rds")
}

# The numbers an analyst reads from a pollin_draws object.
golden_draws <- function(x) {
  list(
    summary = summary(x),
    high = x$high,
    low = x$low,
    estimate = x$estimate,
    method = x$method,
    engine = x$engine,
    n = x$n,
    seed = x$seed
  )
}

# The numbers an analyst reads from a pollin_value or pollin_value_total.
golden_value <- function(x) {
  list(draws = x$draws, summary = x$summary)
}

# Every 0.2.0 output the Stage 3 golden test compares, on the 0.1.0 inputs.
golden_capture_v020 <- function() {
  d <- golden_inputs_v010()
  pd <- pollin_data(d, "pods", "trt", unit = "Block")
  ds <- pollin_design(pd, golden_arms_v010, "rcbd")
  e_dep <- pollin_estimand(ds, "dependence", "ratio", "loss")
  e_def <- pollin_estimand(ds, "deficit", "ratio", "gain")
  pop <- pollin_population(
    "Region",
    c(north = 1200, south = 800),
    label = "two regions"
  )
  e_pop <- pollin_estimand(ds, "dependence", "ratio", "loss", population = pop)

  fit_dep <- suppressWarnings(pollin_fit(e_dep))
  fit_dep_by <- suppressWarnings(pollin_fit(e_dep, by = "Region"))
  fit_def <- suppressWarnings(pollin_fit(e_def))
  fit_pop <- suppressWarnings(pollin_fit(e_pop, by = "Region"))

  dr_sim <- pollin_draws(fit_dep, n = 1000L, seed = 41L)
  dr_by <- pollin_draws(fit_dep_by, n = 1000L, seed = 42L)
  dr_boot <- suppressWarnings(
    pollin_draws(fit_dep, n = 19L, method = "bootstrap", seed = 43L)
  )
  dr_pop <- pollin_draws(fit_pop, n = 1000L, seed = 44L)
  dr_def <- pollin_draws(fit_def, n = 1000L, seed = 45L)

  canola <- pollin_market(
    1000,
    500,
    "t",
    "AUD",
    2024,
    "Synthetic inputs",
    crop = "canola"
  )
  faba <- pollin_market(
    400,
    650,
    "t",
    "AUD",
    2024,
    "Synthetic inputs",
    price_draws = stats::setNames(
      650 * exp(stats::qnorm(stats::ppoints(1000L), 0, 0.1)),
      NULL
    ),
    crop = "faba bean"
  )
  v_dep <- pollin_value(dr_pop, canola, channel = "dependence")
  v_def <- pollin_value(dr_def, faba, channel = "deficit")
  tot <- pollin_guard_doublecount(v_dep, v_def)

  list(
    draws_simulation = golden_draws(dr_sim),
    draws_by = golden_draws(dr_by),
    draws_bootstrap = golden_draws(dr_boot),
    draws_population = golden_draws(dr_pop),
    draws_deficit = golden_draws(dr_def),
    value_dependence = golden_value(v_dep),
    value_deficit = golden_value(v_def),
    total = golden_value(tot)
  )
}

# Outputs at 0.4.0: the joint delta interval and the joint simulation draws ----

golden_path_v040 <- function() {
  testthat::test_path("_golden", "golden_v0.4.0.rds")
}

# The joint delta effects on the 0.1.0 inputs, every scale, overall and by
# region, and the contributions; then every 0.2.0 output.
golden_capture_v040 <- function() {
  d <- golden_inputs_v010()
  pd <- pollin_data(
    d,
    c("pods", "seed_wt"),
    "trt",
    unit = "Block",
    groups = "Region",
    exposures = "bees"
  )
  ds <- pollin_design(pd, golden_arms_v010, "rcbd")
  grid <- golden_grid_v010()
  delta <- list()
  for (fam in c("deficit", "dependence")) {
    for (i in seq_len(nrow(grid))) {
      key <- paste(fam, grid$scale[i], grid$basis[i], sep = "_")
      est <- pollin_estimand(ds, fam, grid$scale[i], grid$basis[i])
      delta[[key]] <- list(
        delta = golden_effects(pollin_effects(pollin_fit(est))),
        delta_by = golden_effects(
          pollin_effects(pollin_fit(est, by = "Region"))
        )
      )
    }
  }
  e_diff <- pollin_estimand(ds, "dependence", "difference", "loss")
  c(
    list(
      delta = delta,
      contribution = pollin_contribution(e_diff)$table,
      contribution_by = pollin_contribution(e_diff, by = "Region")$table
    ),
    golden_capture_v020()
  )
}

# Outputs at 0.4.0: moderation ------------------------------------------------

golden_path_v040_moderation <- function() {
  testthat::test_path("_golden", "golden_v0.4.0_moderation.rds")
}

# The numbers an analyst reads from a pollin_moderation object. A fit without a
# positive-definite Hessian or without convergence is never an oracle.
golden_moderation <- function(x) {
  if (!isTRUE(x$fit$sdr$pdHess) || x$fit$fit$convergence != 0L) {
    stop("Moderation fit for `", x$outcome, "` did not converge cleanly.",
      call. = FALSE)
  }
  list(table = x$table, outcome = x$outcome)
}

# Moderation on two designs that Stage 4 must leave unchanged: the shipped
# example (exposure varying within site) and the three-arm block trial of
# 0.1.0. Not captured: the nested design of 0.0.0.9000 (its fits lack a
# positive-definite Hessian; `pods` stays frozen by the older oracles) and a
# site-level-only exposure (its 0.4.0 output is the defect Stage 4 fixes).
golden_capture_v040_moderation <- function() {
  set.seed(1)
  site <- rep(paste0("s", seq_len(14L)), each = 20L)
  trt <- rep(rep(c("Cage", "Open"), each = 10L), 14L)
  bees <- rep(stats::rpois(14L, 8), each = 20L) + stats::rpois(280L, 1)
  y <- stats::rpois(
    length(trt),
    exp(1.6 + 0.3 * (trt == "Open") + 0.02 * bees * (trt == "Open"))
  )
  d_ex <- data.frame(y = y, trt = trt, Site = site, bees = bees)
  pd_ex <- pollin_data(d_ex, "y", "trt", unit = "Site", exposures = "bees")
  ds_ex <- pollin_design(
    pd_ex,
    c(excluded = "Cage", open = "Open"),
    "multi_site"
  )

  d3 <- golden_inputs_v010()
  pd3 <- pollin_data(
    d3,
    c("pods", "seed_wt"),
    "trt",
    unit = "Block",
    groups = "Region",
    exposures = "bees"
  )
  ds3 <- pollin_design(pd3, golden_arms_v010, "rcbd")

  list(
    example = golden_moderation(pollin_moderation(ds_ex, outcome = "y")),
    rcbd_pods = golden_moderation(pollin_moderation(ds3, outcome = "pods")),
    rcbd_seed_wt = golden_moderation(
      pollin_moderation(ds3, outcome = "seed_wt")
    )
  )
}

# Outputs at 0.4.0: CATE -------------------------------------------------------

golden_path_v040_cate <- function() {
  testthat::test_path("_golden", "golden_v0.4.0_cate.rds")
}

# The shipped `pollin_cate()` example, glmmTMB engine, on two scales.
golden_capture_v040_cate <- function() {
  set.seed(3)
  d <- data.frame(
    trt = rep(rep(c("Cage", "Open"), each = 8L), 12L),
    Site = rep(sprintf("s%02d", seq_len(12L)), each = 16L)
  )
  slope <- stats::rnorm(12L, 0, 0.3)[as.integer(factor(d$Site))]
  d$pods <- stats::rnbinom(
    nrow(d),
    mu = exp(3 + (0.3 + slope) * (d$trt == "Open")),
    size = 20
  )
  pd <- pollin_data(d, "pods", "trt", unit = "Site")
  ds <- pollin_design(pd, c(excluded = "Cage", open = "Open"), "multi_site")
  est_r <- pollin_estimand(ds, "dependence", "ratio", "loss")
  est_d <- pollin_estimand(ds, "dependence", "difference", "loss")
  list(
    ratio_loss = pollin_cate(pollin_fit(est_r), n = 1000L, seed = 1L)$table,
    difference = pollin_cate(pollin_fit(est_d), n = 1000L, seed = 1L)$table
  )
}
