# Planted-truth recovery suite, Stage 3 (criteria C2 engine agreement, C3
# Bayesian calibration, C4 shrinkage recovery).
#
# Implements validation/recovery_rule_S3.md exactly. Run from the package root:
#
#   Rscript --vanilla validation/planted_truth_S3.R [n_seeds] [n_seeds_c4] [cores] [out_dir]
#
# Defaults: 200 seeds per C2/C3 cell, 100 per C4 cell, 8 cores,
# out_dir = "validation". A smoke run uses smaller seed counts and a scratch
# out_dir; only the full run writes the results named in the rule.

t_start <- Sys.time()
args <- commandArgs(trailingOnly = TRUE)
n_seeds <- if (length(args) >= 1L) as.integer(args[1L]) else 200L
n_seeds_c4 <- if (length(args) >= 2L) as.integer(args[2L]) else 100L
cores <- if (length(args) >= 3L) as.integer(args[3L]) else 8L
out_dir <- if (length(args) >= 4L) args[4L] else "validation"
stopifnot(n_seeds >= 1L, n_seeds_c4 >= 0L, cores >= 1L, dir.exists(out_dir))

suppressMessages(devtools::load_all(".", quiet = TRUE))
options(
  flexyBayes.inla_best_of = 1L,
  flexyBayes.silence_default_prior_note = TRUE
)
INLA::inla.setOption(num.threads = "1:1")

n_draws <- 1000L
conf_level <- 0.95
me_advisory <- "only takes into account the uncertainty in fixed-effect"

# Designs (rule Section 2) -----------------------------------------------------
# Each simulator returns the data, the linear predictor of every row under
# each arm at the realised random effects (`eta`, a matrix with one column per
# arm), and the response drawn from the row's own arm.

rnb <- function(mu, theta) stats::rnbinom(length(mu), mu = mu, size = theta)
rgam <- function(mu, shape) stats::rgamma(length(mu), shape = shape, rate = shape / mu)
rbern <- function(p) stats::rbinom(length(p), 1L, p)

# Multi-site layout: site intercept plus an independent site-by-arm effect
sim_multi_site <- function(n_site, n_rep, arms, beta, sd_site, sd_cell, draw, linkinv) {
  sites <- sprintf("s%02d", seq_len(n_site))
  g <- expand.grid(
    rep = seq_len(n_rep),
    trt = arms,
    site = sites,
    stringsAsFactors = FALSE
  )
  u_site <- stats::rnorm(n_site, 0, sd_site)
  names(u_site) <- sites
  u_cell <- matrix(
    stats::rnorm(n_site * length(arms), 0, sd_cell),
    n_site,
    dimnames = list(sites, arms)
  )
  eta <- sapply(arms, function(a) {
    beta[[a]] + u_site[g$site] + u_cell[cbind(g$site, a)]
  })
  own <- eta[cbind(seq_len(nrow(g)), match(g$trt, arms))]
  g$y <- draw(linkinv(own))
  list(data = g, eta = eta)
}

# Randomised blocks: block intercept only
sim_rcbd <- function(n_block, n_rep, arms, beta, sd_block, draw, linkinv) {
  blocks <- sprintf("b%02d", seq_len(n_block))
  g <- expand.grid(
    rep = seq_len(n_rep),
    trt = arms,
    block = blocks,
    stringsAsFactors = FALSE
  )
  u_block <- stats::rnorm(n_block, 0, sd_block)
  names(u_block) <- blocks
  eta <- sapply(arms, function(a) beta[[a]] + u_block[g$block])
  own <- eta[cbind(seq_len(nrow(g)), match(g$trt, arms))]
  g$y <- draw(linkinv(own))
  list(data = g, eta = eta)
}

# Varieties within blocks: block intercept, variety intercept and an
# independent variety-by-arm effect
sim_variety <- function(n_variety, n_block, n_rep, arms, beta, sd_block, sd_var, sd_cell, theta) {
  blocks <- sprintf("b%02d", seq_len(n_block))
  vars <- sprintf("v%02d", seq_len(n_variety))
  g <- expand.grid(
    rep = seq_len(n_rep),
    trt = arms,
    variety = vars,
    block = blocks,
    stringsAsFactors = FALSE
  )
  u_block <- stats::rnorm(n_block, 0, sd_block)
  names(u_block) <- blocks
  u_var <- stats::rnorm(n_variety, 0, sd_var)
  names(u_var) <- vars
  u_cell <- matrix(
    stats::rnorm(n_variety * length(arms), 0, sd_cell),
    n_variety,
    dimnames = list(vars, arms)
  )
  eta <- sapply(arms, function(a) {
    beta[[a]] + u_block[g$block] + u_var[g$variety] + u_cell[cbind(g$variety, a)]
  })
  own <- eta[cbind(seq_len(nrow(g)), match(g$trt, arms))]
  g$y <- rnb(exp(own), theta)
  list(data = g, eta = eta)
}

two <- c("Cage", "Open")
three <- c("Cage", "Open", "Hand")
arms_two <- c(excluded = "Cage", open = "Open")
arms_three <- c(excluded = "Cage", open = "Open", supplemented = "Hand")

# C2/C3 data designs: k is the seed block; A6 reuses the data of A1
data_designs <- list(
  A1 = list(
    k = 1L,
    simulate = function() {
      sim_multi_site(12L, 8L, two, log(c(Cage = 20, Open = 30)), 0.40, 0.20,
                     function(mu) rnb(mu, 8), exp)
    },
    unit = "site", layout = "multi_site", arms = arms_two
  ),
  A2 = list(
    k = 2L,
    simulate = function() {
      sim_rcbd(8L, 6L, two, log(c(Cage = 20, Open = 30)), 0.40,
               function(mu) stats::rpois(length(mu), mu), exp)
    },
    unit = "block", layout = "rcbd", arms = arms_two
  ),
  A3 = list(
    k = 3L,
    simulate = function() {
      sim_multi_site(12L, 8L, two, log(c(Cage = 2.0, Open = 2.6)), 0.30, 0.20,
                     function(mu) rgam(mu, 10), exp)
    },
    unit = "site", layout = "multi_site", arms = arms_two
  ),
  A4 = list(
    k = 4L,
    simulate = function() {
      sim_rcbd(8L, 24L, two, stats::qlogis(c(Cage = 0.55, Open = 0.75)), 0.40,
               rbern, stats::plogis)
    },
    unit = "block", layout = "rcbd", arms = arms_two
  ),
  A5 = list(
    k = 5L,
    simulate = function() {
      sim_rcbd(8L, 4L, three, log(c(Cage = 20, Open = 30, Hand = 36)), 0.40,
               function(mu) rnb(mu, 8), exp)
    },
    unit = "block", layout = "rcbd", arms = arms_three
  )
)

# C2/C3 cells (rule Section 3)
cells <- list(
  A1 = list(data = "A1", family = glmmTMB::nbinom2(), linkinv = exp,
            channel = "dependence", scale = "ratio", basis = "loss"),
  A2 = list(data = "A2", family = stats::poisson(), linkinv = exp,
            channel = "dependence", scale = "ratio", basis = "loss"),
  A3 = list(data = "A3", family = stats::Gamma(link = "log"), linkinv = exp,
            channel = "dependence", scale = "ratio", basis = "loss"),
  A4 = list(data = "A4", family = stats::binomial(), linkinv = stats::plogis,
            channel = "dependence", scale = "ratio", basis = "loss"),
  A5 = list(data = "A5", family = glmmTMB::nbinom2(), linkinv = exp,
            channel = "deficit", scale = "ratio", basis = "gain"),
  A6 = list(data = "A1", family = glmmTMB::nbinom2(), linkinv = exp,
            channel = "dependence", scale = "difference", basis = "loss"),
  P8 = list(data = "S8", family = glmmTMB::nbinom2(), linkinv = exp,
            channel = "dependence", scale = "ratio", basis = "loss"),
  P12 = list(data = "S12", family = glmmTMB::nbinom2(), linkinv = exp,
             channel = "dependence", scale = "ratio", basis = "loss")
)
# Cells whose fit carries a random treatment slope (rule Section 3): the
# width ratio is reported there, not gated
slope_cells <- c("A1", "A3", "A6", "P8", "P12")

# C4 cells (rule Section 4)
c4_designs <- list(
  S8 = list(k = 7L, n_variety = 8L, overall_cell = "P8"),
  S12 = list(k = 8L, n_variety = 12L, overall_cell = "P12")
)
c4_sim <- function(des) {
  sim_variety(des$n_variety, 4L, 3L, two, log(c(Cage = 20, Open = 30)),
              0.30, 0.30, 0.20, 8)
}

contrast_of <- function(arms, channel) {
  if (channel == "deficit") {
    c(high = arms[["supplemented"]], low = arms[["open"]])
  } else {
    c(high = arms[["open"]], low = arms[["excluded"]])
  }
}

effect_of <- function(H, L, scale, basis) {
  if (scale == "difference") {
    return(H - L)
  }
  if (basis == "loss") 1 - L / H else H / L - 1
}

# Truth: arm means at the realised random effects, averaged over the rows
# (overall, or per level of `g`)
truth_of <- function(sim, hl, linkinv, scale, basis, g = NULL) {
  h <- linkinv(sim$eta[, hl[["high"]]])
  l <- linkinv(sim$eta[, hl[["low"]]])
  if (is.null(g)) {
    return(effect_of(mean(h), mean(l), scale, basis))
  }
  effect_of(tapply(h, g, mean), tapply(l, g, mean), scale, basis)
}

seed_of <- function(k, s, offset) 3000000L + offset + 100000L * k + s

# Running one engine, recording errors, warnings and convergence -----------------

run_engine <- function(expr_fun) {
  warns <- character(0)
  res <- tryCatch(
    withCallingHandlers(
      expr_fun(),
      warning = function(w) {
        msg <- conditionMessage(w)
        if (!grepl(me_advisory, msg, fixed = TRUE)) {
          warns <<- c(warns, msg)
        }
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) structure(conditionMessage(e), class = "run_error")
  )
  err <- if (inherits(res, "run_error")) as.character(unclass(res)) else NA_character_
  # A flexyBayes fit with no convergence verdict is counted as not converged
  if (is.na(err) && any(grepl("gives no convergence verdict", warns, fixed = TRUE))) {
    err <- "flexyBayes: no convergence verdict"
  }
  list(
    value = if (is.na(err)) res else NULL,
    error = err,
    n_warn = length(warns),
    first_warn = if (length(warns)) warns[1L] else NA_character_
  )
}

# glmmTMB convergence: optimiser code 0 and a positive-definite Hessian
glmm_converged <- function(fit) {
  m <- fit$fits[[1L]]
  isTRUE(m$fit$convergence == 0L) && isTRUE(m$sdr$pdHess)
}

build_estimand <- function(d, unit, layout, arms, channel, scale, basis) {
  pd <- pollin_data(d, "y", "trt", unit = unit)
  ds <- pollin_design(pd, arms, layout)
  pollin_estimand(ds, channel, scale, basis)
}

# C2/C3 task: one data set, its cells, both engines ----------------------------

run_c2_task <- function(task) {
  dd <- data_designs[[task$data]]
  s <- task$seed
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed_of(dd$k, s, 0L))
  sim <- dd$simulate()
  ids <- names(cells)[vapply(cells, function(z) z$data == task$data, logical(1))]
  rows <- lapply(ids, function(id) {
    cl <- cells[[id]]
    hl <- contrast_of(dd$arms, cl$channel)
    truth <- truth_of(sim, hl, cl$linkinv, cl$scale, cl$basis)
    est <- build_estimand(sim$data, dd$unit, dd$layout, dd$arms,
                          cl$channel, cl$scale, cl$basis)
    t0 <- Sys.time()
    g <- run_engine(function() {
      fit <- pollin_fit(est, engine = "glmmTMB", family = cl$family)
      if (!glmm_converged(fit)) stop("glmmTMB: not converged")
      pollin_effects(fit, interval = "delta", conf_level = conf_level)$table
    })
    t1 <- Sys.time()
    b <- run_engine(function() {
      fit <- pollin_fit(est, engine = "flexyBayes", family = cl$family,
                        aggregate = FALSE)
      pollin_effects(fit, interval = "posterior", n_draws = n_draws,
                     seed = seed_of(dd$k, s, 20000000L),
                     conf_level = conf_level)$table
    })
    t2 <- Sys.time()
    pick <- function(r, col) if (is.null(r$value)) NA_real_ else r$value[[col]][1L]
    data.frame(
      part = "C2C3", cell = id, seed = s, level = NA_character_,
      truth = truth,
      g_estimate = pick(g, "estimate"), g_se = pick(g, "std.error"),
      g_low = pick(g, "conf.low"), g_high = pick(g, "conf.high"),
      b_estimate = pick(b, "estimate"), b_se = pick(b, "std.error"),
      b_low = pick(b, "conf.low"), b_high = pick(b, "conf.high"),
      g_error = g$error, g_n_warn = g$n_warn, g_first_warn = g$first_warn,
      b_error = b$error, b_n_warn = b$n_warn, b_first_warn = b$first_warn,
      g_secs = as.numeric(difftime(t1, t0, units = "secs")),
      b_secs = as.numeric(difftime(t2, t1, units = "secs")),
      stringsAsFactors = FALSE
    )
  })
  list(c2 = do.call(rbind, rows), c4 = NULL)
}

# C4 task: pooled and unpooled fits, both engines -------------------------------

run_c4_task <- function(task) {
  des <- c4_designs[[task$cell]]
  s <- task$seed
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed_of(des$k, s, 0L))
  sim <- c4_sim(des)
  d <- sim$data
  hl <- contrast_of(arms_two, "dependence")
  truth <- truth_of(sim, hl, exp, "ratio", "loss", g = d$variety)
  est <- build_estimand(d, "block", "rcbd", arms_two, "dependence", "ratio", "loss")
  fam <- glmmTMB::nbinom2()
  draw_seed <- seed_of(des$k, s, 20000000L)
  tabs <- list()
  t0 <- Sys.time()
  tabs$gp <- run_engine(function() {
    fit <- pollin_fit(est, engine = "glmmTMB", family = fam, pooled = "variety")
    if (!glmm_converged(fit)) stop("glmmTMB: not converged")
    list(
      cate = pollin_cate(fit, n = n_draws, seed = draw_seed, conf_level = conf_level)$table,
      overall = pollin_effects(fit, interval = "delta", conf_level = conf_level)$table
    )
  })
  tabs$gu <- run_engine(function() {
    fit <- pollin_fit(est, engine = "glmmTMB", family = fam, by = "variety")
    if (!glmm_converged(fit)) stop("glmmTMB: not converged")
    pollin_effects(fit, interval = "delta", conf_level = conf_level)$table
  })
  t1 <- Sys.time()
  tabs$bp <- run_engine(function() {
    fit <- pollin_fit(est, engine = "flexyBayes", family = fam, pooled = "variety",
                      aggregate = FALSE)
    list(
      cate = pollin_cate(fit, n = n_draws, seed = draw_seed, conf_level = conf_level)$table,
      overall = pollin_effects(fit, interval = "posterior", n_draws = n_draws,
                               seed = draw_seed, conf_level = conf_level)$table
    )
  })
  tabs$bu <- run_engine(function() {
    fit <- pollin_fit(est, engine = "flexyBayes", family = fam, by = "variety",
                      aggregate = FALSE)
    pollin_effects(fit, interval = "posterior", n_draws = n_draws,
                   seed = draw_seed, conf_level = conf_level)$table
  })
  t2 <- Sys.time()
  lv <- names(truth)
  out <- data.frame(
    part = "C4", cell = task$cell, seed = s, level = lv,
    truth = unname(truth), stringsAsFactors = FALSE
  )
  # Overall estimand on the pooled fits: the C2/C3 cells P8 and P12
  pick <- function(r, col) {
    if (is.null(r$value)) NA_real_ else r$value$overall[[col]][1L]
  }
  g <- tabs$gp
  b <- tabs$bp
  c2 <- data.frame(
    part = "C2C3", cell = des$overall_cell, seed = s, level = NA_character_,
    truth = truth_of(sim, hl, exp, "ratio", "loss"),
    g_estimate = pick(g, "estimate"), g_se = pick(g, "std.error"),
    g_low = pick(g, "conf.low"), g_high = pick(g, "conf.high"),
    b_estimate = pick(b, "estimate"), b_se = pick(b, "std.error"),
    b_low = pick(b, "conf.low"), b_high = pick(b, "conf.high"),
    g_error = g$error, g_n_warn = g$n_warn, g_first_warn = g$first_warn,
    b_error = b$error, b_n_warn = b$n_warn, b_first_warn = b$first_warn,
    g_secs = NA_real_, b_secs = NA_real_,
    stringsAsFactors = FALSE
  )
  for (nm in c("gp", "bp")) {
    if (!is.null(tabs[[nm]]$value)) tabs[[nm]]$value <- tabs[[nm]]$value$cate
  }
  for (nm in names(tabs)) {
    r <- tabs[[nm]]
    key <- if (nm %in% c("gp", "bp")) "level" else "variety"
    m <- if (is.null(r$value)) NA_integer_ else match(lv, as.character(r$value[[key]]))
    col <- function(x) if (is.null(r$value)) NA_real_ else r$value[[x]][m]
    out[[paste0(nm, "_estimate")]] <- col("estimate")
    out[[paste0(nm, "_low")]] <- col("conf.low")
    out[[paste0(nm, "_high")]] <- col("conf.high")
    out[[paste0(nm, "_error")]] <- r$error
    out[[paste0(nm, "_n_warn")]] <- r$n_warn
    out[[paste0(nm, "_first_warn")]] <- r$first_warn
  }
  out$g_secs <- as.numeric(difftime(t1, t0, units = "secs"))
  out$b_secs <- as.numeric(difftime(t2, t1, units = "secs"))
  list(c2 = c2, c4 = out)
}

run_task <- function(task) {
  if (task$part == "C2C3") run_c2_task(task) else run_c4_task(task)
}

# Tasks (rule Section 6) -------------------------------------------------------

tasks <- c(
  lapply(names(c4_designs), function(k) {
    lapply(seq_len(n_seeds_c4), function(s) list(part = "C4", cell = k, seed = s))
  }),
  lapply(names(data_designs), function(k) {
    lapply(seq_len(n_seeds), function(s) list(part = "C2C3", data = k, seed = s))
  })
)
tasks <- unlist(tasks, recursive = FALSE)

message("Running ", length(tasks), " tasks on ", cores, " cores ...")
res <- parallel::mclapply(
  tasks,
  run_task,
  mc.cores = cores,
  mc.preschedule = FALSE
)
bad <- vapply(res, inherits, logical(1), what = "try-error")
if (any(bad)) {
  stop(
    sum(bad), " task(s) crashed outside the cell handler; first: ",
    as.character(res[[which(bad)[1L]]])
  )
}
c2_rows <- do.call(rbind, lapply(res, `[[`, "c2"))
c4_rows <- do.call(rbind, lapply(res, `[[`, "c4"))

# Summaries (rule Section 7) ---------------------------------------------------

wilson <- function(x, n, z = stats::qnorm(0.975)) {
  p <- x / n
  centre <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  c(centre - half, centre + half)
}

# PASS inside the band; BORDERLINE outside it by at most 2 MC SE; else FAIL
band_verdict <- function(x, mcse, lo, hi) {
  if (!is.finite(x) || !is.finite(mcse)) {
    return("FAIL")
  }
  if (x >= lo && x <= hi) {
    return("PASS")
  }
  gap <- if (x < lo) lo - x else x - hi
  if (gap <= 2 * mcse) "BORDERLINE" else "FAIL"
}

fail_limit <- function(n) floor(0.02 * n)

c2_summary <- do.call(rbind, lapply(names(cells), function(id) {
  z <- c2_rows[c2_rows$cell == id, ]
  g_ok <- is.na(z$g_error) & is.finite(z$g_estimate) & is.finite(z$g_se) &
    is.finite(z$g_low) & is.finite(z$g_high)
  b_ok <- is.na(z$b_error) & is.finite(z$b_estimate) & is.finite(z$b_se) &
    is.finite(z$b_low) & is.finite(z$b_high)
  ok <- g_ok & b_ok
  zo <- z[ok, ]
  n_ok <- nrow(zo)
  n_g_failed <- sum(!g_ok)
  n_b_failed <- sum(!b_ok)
  lim <- fail_limit(nrow(z))
  fits_ok <- n_g_failed <= lim && n_b_failed <= lim
  agree <- abs(zo$b_estimate - zo$g_estimate) <= 0.25 * zo$g_se
  agree_share <- mean(agree)
  width_ratio <- stats::median((zo$b_high - zo$b_low) / (zo$g_high - zo$g_low))
  agree_met <- isTRUE(agree_share >= 0.95)
  width_gated <- !id %in% slope_cells
  width_met <- !width_gated || isTRUE(width_ratio >= 0.8 && width_ratio <= 1.25)
  c2_pass <- fits_ok && agree_met && width_met
  c2_missed <- c(
    if (!fits_ok) "fits",
    if (!agree_met) "agreement",
    if (!width_met) "width ratio"
  )
  b_cov_ev <- zo$truth >= zo$b_low & zo$truth <= zo$b_high
  g_cov_ev <- zo$truth >= zo$g_low & zo$truth <= zo$g_high
  b_cov <- mean(b_cov_ev)
  b_mcse <- sqrt(b_cov * (1 - b_cov) / n_ok)
  ci <- wilson(sum(b_cov_ev), n_ok)
  c3 <- if (!fits_ok) "FAIL" else band_verdict(b_cov, b_mcse, 0.93, 0.97)
  data.frame(
    cell = id,
    family = cells[[id]]$family$family,
    channel = cells[[id]]$channel,
    scale = cells[[id]]$scale,
    basis = cells[[id]]$basis,
    data = cells[[id]]$data,
    n_seeds = nrow(z),
    n_ok = n_ok,
    n_failed_glmmTMB = n_g_failed,
    n_failed_flexyBayes = n_b_failed,
    n_warned_glmmTMB = sum(z$g_n_warn > 0L),
    n_warned_flexyBayes = sum(z$b_n_warn > 0L),
    truth_mean = mean(z$truth),
    agreement_share = agree_share,
    median_abs_diff_over_se = stats::median(abs(zo$b_estimate - zo$g_estimate) / zo$g_se),
    median_width_ratio = width_ratio,
    width_ratio_gated = width_gated,
    c2_verdict = if (c2_pass) "PASS" else paste0("FAIL (", paste(c2_missed, collapse = ", "), ")"),
    coverage_flexyBayes = b_cov,
    coverage_mcse = b_mcse,
    coverage_lo = ci[1L],
    coverage_hi = ci[2L],
    c3_verdict = c3,
    bias_flexyBayes = mean(zo$b_estimate - zo$truth),
    bias_over_sd_flexyBayes = mean(zo$b_estimate - zo$truth) / stats::sd(zo$b_estimate),
    n_near_coverage = sum(pmin(abs(zo$truth - zo$b_low), abs(zo$truth - zo$b_high)) < 0.15 * zo$b_se),
    n_near_agreement = sum(abs(abs(zo$b_estimate - zo$g_estimate) - 0.25 * zo$g_se) < 0.15 * zo$b_se),
    coverage_glmmTMB = mean(g_cov_ev),
    bias_glmmTMB = mean(zo$g_estimate - zo$truth),
    median_secs_glmmTMB = stats::median(z$g_secs),
    median_secs_flexyBayes = stats::median(z$b_secs),
    stringsAsFactors = FALSE
  )
}))

c4_summary <- if (is.null(c4_rows)) NULL else do.call(rbind, lapply(names(c4_designs), function(id) {
  z <- c4_rows[c4_rows$cell == id, ]
  ok_of <- function(p) {
    e <- tapply(
      is.na(z[[paste0(p, "_error")]]) & is.finite(z[[paste0(p, "_estimate")]]) &
        is.finite(z[[paste0(p, "_low")]]) & is.finite(z[[paste0(p, "_high")]]),
      z$seed,
      all
    )
    e
  }
  okp <- lapply(c(gp = "gp", gu = "gu", bp = "bp", bu = "bu"), ok_of)
  n_seeds_cell <- length(unique(z$seed))
  lim <- fail_limit(n_seeds_cell)
  rmse <- function(p) {
    sq <- (z[[paste0(p, "_estimate")]] - z$truth)^2
    sqrt(tapply(sq, z$seed, mean))
  }
  eng <- function(pp, pu) {
    ok <- okp[[pp]] & okp[[pu]]
    rp <- rmse(pp)[ok]
    ru <- rmse(pu)[ok]
    list(
      n_ok = sum(ok),
      n_failed_pooled = sum(!okp[[pp]]),
      n_failed_unpooled = sum(!okp[[pu]]),
      win = mean(rp < ru),
      rmse_pooled = mean(rp),
      rmse_unpooled = mean(ru),
      fits_ok = sum(!okp[[pp]]) <= lim && sum(!okp[[pu]]) <= lim
    )
  }
  cov_of <- function(p) {
    keep <- okp[[p]][as.character(z$seed)]
    zz <- z[keep, ]
    ev <- zz$truth >= zz[[paste0(p, "_low")]] & zz$truth <= zz[[paste0(p, "_high")]]
    sdl <- (zz[[paste0(p, "_high")]] - zz[[paste0(p, "_low")]]) / (2 * stats::qnorm(0.975))
    near <- pmin(abs(zz$truth - zz[[paste0(p, "_low")]]), abs(zz$truth - zz[[paste0(p, "_high")]])) < 0.15 * sdl
    c(coverage = mean(ev), n = length(ev), x = sum(ev), near = sum(near))
  }
  g <- eng("gp", "gu")
  b <- eng("bp", "bu")
  cb <- cov_of("bp")
  cb_mcse <- sqrt(cb[["coverage"]] * (1 - cb[["coverage"]]) / cb[["n"]])
  ci <- wilson(cb[["x"]], cb[["n"]])
  g_v <- if (g$fits_ok && isTRUE(g$win >= 0.90)) "PASS" else paste0("FAIL (", paste(c(if (!g$fits_ok) "fits", if (!isTRUE(g$win >= 0.90)) "shrinkage"), collapse = ", "), ")")
  b_v <- if (b$fits_ok && isTRUE(b$win >= 0.90)) "PASS" else paste0("FAIL (", paste(c(if (!b$fits_ok) "fits", if (!isTRUE(b$win >= 0.90)) "shrinkage"), collapse = ", "), ")")
  cov_v <- if (!b$fits_ok) "FAIL" else band_verdict(cb[["coverage"]], cb_mcse, 0.92, 0.98)
  data.frame(
    cell = id,
    n_levels = c4_designs[[id]]$n_variety,
    n_seeds = n_seeds_cell,
    n_ok_glmmTMB = g$n_ok,
    n_failed_glmmTMB_pooled = g$n_failed_pooled,
    n_failed_glmmTMB_unpooled = g$n_failed_unpooled,
    n_ok_flexyBayes = b$n_ok,
    n_failed_flexyBayes_pooled = b$n_failed_pooled,
    n_failed_flexyBayes_unpooled = b$n_failed_unpooled,
    win_share_glmmTMB = g$win,
    rmse_pooled_glmmTMB = g$rmse_pooled,
    rmse_unpooled_glmmTMB = g$rmse_unpooled,
    shrinkage_verdict_glmmTMB = g_v,
    win_share_flexyBayes = b$win,
    rmse_pooled_flexyBayes = b$rmse_pooled,
    rmse_unpooled_flexyBayes = b$rmse_unpooled,
    shrinkage_verdict_flexyBayes = b_v,
    level_coverage_flexyBayes = cb[["coverage"]],
    level_coverage_mcse = cb_mcse,
    level_coverage_lo = ci[1L],
    level_coverage_hi = ci[2L],
    level_coverage_verdict = cov_v,
    n_near_level_coverage = cb[["near"]],
    level_coverage_glmmTMB = cov_of("gp")[["coverage"]],
    level_coverage_unpooled_glmmTMB = cov_of("gu")[["coverage"]],
    level_coverage_unpooled_flexyBayes = cov_of("bu")[["coverage"]],
    median_secs_glmmTMB = stats::median(tapply(z$g_secs, z$seed, `[`, 1L)),
    median_secs_flexyBayes = stats::median(tapply(z$b_secs, z$seed, `[`, 1L)),
    stringsAsFactors = FALSE
  )
}))

c2_met <- all(c2_summary$c2_verdict == "PASS")
c3_met <- all(c2_summary$c3_verdict == "PASS")
c4_met <- !is.null(c4_summary) &&
  all(c4_summary$shrinkage_verdict_glmmTMB == "PASS") &&
  all(c4_summary$shrinkage_verdict_flexyBayes == "PASS") &&
  all(c4_summary$level_coverage_verdict == "PASS")

t_end <- Sys.time()
wall <- difftime(t_end, t_start, units = "mins")

# Write ------------------------------------------------------------------------

stack_rows <- function(a, b) {
  if (is.null(b)) {
    return(a)
  }
  cols <- union(names(a), names(b))
  for (x in setdiff(cols, names(a))) a[[x]] <- NA
  for (x in setdiff(cols, names(b))) b[[x]] <- NA
  rbind(a[cols], b[cols])
}
seeds_out <- stack_rows(c2_rows, c4_rows)
seeds_out <- seeds_out[order(seeds_out$part, seeds_out$cell, seeds_out$seed, seeds_out$level), ]
rownames(seeds_out) <- NULL

results <- stack_rows(
  cbind(data.frame(criterion = "C2/C3", stringsAsFactors = FALSE), c2_summary),
  if (is.null(c4_summary)) NULL else cbind(data.frame(criterion = "C4", stringsAsFactors = FALSE), c4_summary)
)
utils::write.csv(results, file.path(out_dir, "planted_truth_S3_results.csv"), row.names = FALSE)
gz <- gzfile(file.path(out_dir, "planted_truth_S3_seeds.csv.gz"), "w")
utils::write.csv(seeds_out, gz, row.names = FALSE)
close(gz)

err_cols <- grep("_error$", names(seeds_out), value = TRUE)
warn_cols <- grep("_first_warn$", names(seeds_out), value = TRUE)
sink(file.path(out_dir, "planted_truth_S3_session.txt"))
cat("Planted-truth suite S3 (rule: validation/recovery_rule_S3.md)\n")
cat("Started :", format(t_start, "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Finished:", format(t_end, "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Wall time (min):", round(as.numeric(wall), 2), "\n")
cat(
  "Seeds per C2/C3 cell:", n_seeds, "; per C4 cell:", n_seeds_c4,
  "; cores:", cores, "; tasks:", length(tasks), "\n"
)
cat("Posterior and interval draws per fit:", n_draws, "\n")
cat("C2 met:", c2_met, "; C3 met:", c3_met, "; C4 met:", c4_met, "\n")
cat("pollinatoR git HEAD:", system("git rev-parse HEAD", intern = TRUE), "\n")
cat("INLA num.threads:", INLA::inla.getOption("num.threads"),
    "; flexyBayes.inla_best_of:", getOption("flexyBayes.inla_best_of"), "\n")
cat("\nPackage versions:\n")
for (p in c("pollinatoR", "flexyBayes", "INLA", "glmmTMB", "TMB", "marginaleffects", "Matrix")) {
  cat(sprintf("  %-16s %s\n", p, as.character(utils::packageVersion(p))))
}
for (x in err_cols) {
  v <- seeds_out[[x]]
  cat("\nErrors in", x, "(message, rows):\n")
  if (any(!is.na(v))) print(table(v[!is.na(v)])) else cat("  none\n")
}
for (x in warn_cols) {
  v <- seeds_out[[x]]
  cat("\nWarnings in", x, "(first per row, rows; marginaleffects advisory excluded):\n")
  if (any(!is.na(v))) print(table(v[!is.na(v)])) else cat("  none\n")
}
cat("\n")
print(utils::sessionInfo())
sink()

print(c2_summary[c("cell", "n_ok", "n_failed_glmmTMB", "n_failed_flexyBayes",
                   "agreement_share", "median_width_ratio", "c2_verdict",
                   "coverage_flexyBayes", "coverage_mcse", "c3_verdict")], digits = 3)
if (!is.null(c4_summary)) print(c4_summary, digits = 3)
message("C2 met: ", c2_met, "; C3 met: ", c3_met, "; C4 met: ", c4_met)
message("Wall time: ", round(as.numeric(wall), 2), " min")
