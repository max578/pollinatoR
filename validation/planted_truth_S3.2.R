# Planted-truth recovery suite, Stage 3.2 (the joint delta interval of 0.4.0:
# gates G3 on B1, A2, A5 and G4 on the slope cells A1, A3, A1W).
#
# Implements validation/recovery_rule_S3.2.md exactly. Run from the package root:
#
#   Rscript --vanilla validation/planted_truth_S3.2.R [scale] [cores] [out_dir]
#
# `scale` multiplies every cell's seed count (1 for the full run; a smoke run
# uses a small fraction and a scratch out_dir). Defaults: 1, 8 cores,
# out_dir = "validation".

t_start <- Sys.time()
args <- commandArgs(trailingOnly = TRUE)
scale <- if (length(args) >= 1L) as.numeric(args[1L]) else 1
cores <- if (length(args) >= 2L) as.integer(args[2L]) else 8L
out_dir <- if (length(args) >= 3L) args[3L] else "validation"
stopifnot(scale > 0, cores >= 1L, dir.exists(out_dir))

suppressMessages(devtools::load_all(".", quiet = TRUE))
options(
  flexyBayes.inla_best_of = 1L,
  flexyBayes.silence_default_prior_note = TRUE
)
INLA::inla.setOption(num.threads = "1:1")

n_draws <- 1000L
conf_level <- 0.95
me_advisory <- "only takes into account the uncertainty in fixed-effect"
fallback_note <- "holds the variance parameters at their estimates"

# Designs (rule Section 1) -----------------------------------------------------
# Each simulator returns the data and the linear predictor of every row under
# each arm at the realised random effects (`eta`, one column per arm).

rnb <- function(mu, theta) stats::rnbinom(length(mu), mu = mu, size = theta)

# Randomised blocks: block intercept only (A2, A5; as S3)
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

# Randomised blocks, seeds counted out of ovules per plant (B1): block
# intercept, ovules drawn uniformly on 20 to 40, seeds binomial out of them
sim_rcbd_trials <- function(n_block, n_rep, arms, beta, sd_block, ovules) {
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
  g$n <- sample(ovules, nrow(g), replace = TRUE)
  g$y <- stats::rbinom(nrow(g), g$n, stats::plogis(own))
  list(data = g, eta = eta)
}

two <- c("Cage", "Open")
three <- c("Cage", "Open", "Hand")
arms_two <- c(excluded = "Cage", open = "Open")
arms_three <- c(excluded = "Cage", open = "Open", supplemented = "Hand")

# Multi-site layout: site intercept plus an independent site-by-arm effect (as S3)
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
rgam <- function(mu, shape) stats::rgamma(length(mu), shape = shape, rate = shape / mu)

n_of <- function(n) max(1L, as.integer(round(scale * n)))

# Cells (rule Section 1): k is the seed index
cells <- list(
  B1 = list(
    k = 1L, n_seeds = n_of(200L),
    simulate = function() {
      sim_rcbd_trials(8L, 6L, two, stats::qlogis(c(Cage = 0.55, Open = 0.75)),
                      0.40, 20:40)
    },
    unit = "block", layout = "rcbd", arms = arms_two, denominator = "n",
    family = stats::binomial(), linkinv = stats::plogis,
    channel = "dependence", scale = "ratio", basis = "loss", slope = FALSE
  ),
  A2 = list(
    k = 2L, n_seeds = n_of(100L),
    simulate = function() {
      sim_rcbd(8L, 6L, two, log(c(Cage = 20, Open = 30)), 0.40,
               function(mu) stats::rpois(length(mu), mu), exp)
    },
    unit = "block", layout = "rcbd", arms = arms_two, denominator = NULL,
    family = stats::poisson(), linkinv = exp,
    channel = "dependence", scale = "ratio", basis = "loss", slope = FALSE
  ),
  A5 = list(
    k = 3L, n_seeds = n_of(100L),
    simulate = function() {
      sim_rcbd(8L, 4L, three, log(c(Cage = 20, Open = 30, Hand = 36)), 0.40,
               function(mu) rnb(mu, 8), exp)
    },
    unit = "block", layout = "rcbd", arms = arms_three, denominator = NULL,
    family = glmmTMB::nbinom2(), linkinv = exp,
    channel = "deficit", scale = "ratio", basis = "gain", slope = FALSE
  ),
  A1 = list(
    k = 4L, n_seeds = n_of(100L),
    simulate = function() {
      sim_multi_site(12L, 8L, two, log(c(Cage = 20, Open = 30)), 0.40, 0.20,
                     function(mu) rnb(mu, 8), exp)
    },
    unit = "site", layout = "multi_site", arms = arms_two, denominator = NULL,
    family = glmmTMB::nbinom2(), linkinv = exp,
    channel = "dependence", scale = "ratio", basis = "loss", slope = TRUE
  ),
  A3 = list(
    k = 5L, n_seeds = n_of(100L),
    simulate = function() {
      sim_multi_site(12L, 8L, two, log(c(Cage = 2.0, Open = 2.6)), 0.30, 0.20,
                     function(mu) rgam(mu, 10), exp)
    },
    unit = "site", layout = "multi_site", arms = arms_two, denominator = NULL,
    family = stats::Gamma(link = "log"), linkinv = exp,
    channel = "dependence", scale = "ratio", basis = "loss", slope = TRUE
  ),
  A1W = list(
    k = 6L, n_seeds = n_of(200L),
    simulate = function() {
      sim_multi_site(24L, 8L, two, log(c(Cage = 20, Open = 30)), 0.40, 0.20,
                     function(mu) rnb(mu, 8), exp)
    },
    unit = "site", layout = "multi_site", arms = arms_two, denominator = NULL,
    family = glmmTMB::nbinom2(), linkinv = exp,
    channel = "dependence", scale = "ratio", basis = "loss", slope = TRUE
  )
)

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

truth_of <- function(sim, hl, linkinv, scale, basis) {
  h <- linkinv(sim$eta[, hl[["high"]]])
  l <- linkinv(sim$eta[, hl[["low"]]])
  effect_of(mean(h), mean(l), scale, basis)
}

# Seeds (rule Section 3)
seed_data <- function(k, s) 51000000L + 100000L * k + s
seed_draw <- function(k, s) 71000000L + 100000L * k + s

# Running one step, recording errors and warnings -------------------------------

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
  if (is.na(err) && any(grepl("gives no convergence verdict", warns, fixed = TRUE))) {
    err <- "flexyBayes: no convergence verdict"
  }
  list(
    value = if (is.na(err)) res else NULL,
    error = err,
    n_warn = length(warns),
    first_warn = if (length(warns)) warns[1L] else NA_character_,
    fallback = any(grepl(fallback_note, warns, fixed = TRUE))
  )
}

glmm_converged <- function(fit) {
  m <- fit$fits[[1L]]
  isTRUE(m$fit$convergence == 0L) && isTRUE(m$sdr$pdHess)
}

build_estimand <- function(d, cl) {
  pd <- pollin_data(d, "y", "trt", unit = cl$unit, denominator = cl$denominator)
  ds <- pollin_design(pd, cl$arms, cl$layout)
  pollin_estimand(ds, cl$channel, cl$scale, cl$basis)
}

# Task: one data set, both engines, both glmmTMB intervals ----------------------

run_task <- function(task) {
  id <- task$cell
  cl <- cells[[id]]
  s <- task$seed
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed_data(cl$k, s))
  sim <- cl$simulate()
  hl <- contrast_of(cl$arms, cl$channel)
  truth <- truth_of(sim, hl, cl$linkinv, cl$scale, cl$basis)
  est <- build_estimand(sim$data, cl)
  t0 <- Sys.time()
  gfit <- run_engine(function() {
    fit <- pollin_fit(est, engine = "glmmTMB", family = cl$family)
    if (!glmm_converged(fit)) stop("glmmTMB: not converged")
    fit
  })
  g <- if (is.null(gfit$value)) gfit else run_engine(function() {
    pollin_effects(gfit$value, interval = "delta", conf_level = conf_level)$table
  })
  gf <- if (is.null(gfit$value)) gfit else run_engine(function() {
    pollin_effects(gfit$value, interval = "delta_fixed", conf_level = conf_level)$table
  })
  t1 <- Sys.time()
  b <- run_engine(function() {
    fit <- pollin_fit(est, engine = "flexyBayes", family = cl$family,
                      aggregate = FALSE)
    pollin_effects(fit, interval = "posterior", n_draws = n_draws,
                   seed = seed_draw(cl$k, s), conf_level = conf_level)$table
  })
  t2 <- Sys.time()
  pick <- function(r, col) if (is.null(r$value)) NA_real_ else r$value[[col]][1L]
  data.frame(
    cell = id, seed = s,
    n_rows = nrow(sim$data),
    truth = truth,
    g_estimate = pick(g, "estimate"), g_se = pick(g, "std.error"),
    g_low = pick(g, "conf.low"), g_high = pick(g, "conf.high"),
    f_estimate = pick(gf, "estimate"), f_se = pick(gf, "std.error"),
    f_low = pick(gf, "conf.low"), f_high = pick(gf, "conf.high"),
    b_estimate = pick(b, "estimate"), b_se = pick(b, "std.error"),
    b_low = pick(b, "conf.low"), b_high = pick(b, "conf.high"),
    g_error = g$error, g_n_warn = g$n_warn + gfit$n_warn,
    g_first_warn = if (is.na(gfit$first_warn)) g$first_warn else gfit$first_warn,
    g_fallback = isTRUE(g$fallback),
    f_error = gf$error,
    b_error = b$error, b_n_warn = b$n_warn, b_first_warn = b$first_warn,
    g_secs = as.numeric(difftime(t1, t0, units = "secs")),
    b_secs = as.numeric(difftime(t2, t1, units = "secs")),
    stringsAsFactors = FALSE
  )
}

tasks <- unlist(lapply(names(cells), function(id) {
  lapply(seq_len(cells[[id]]$n_seeds), function(s) list(cell = id, seed = s))
}), recursive = FALSE)

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
rows <- do.call(rbind, res)

# Summaries (rule Section 2) ---------------------------------------------------

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

cover <- function(truth, lo, hi) truth >= lo & truth <= hi

summary_tab <- do.call(rbind, lapply(names(cells), function(id) {
  cl <- cells[[id]]
  z <- rows[rows$cell == id, ]
  if (nrow(z) == 0L) {
    return(NULL)
  }
  fin <- function(p) {
    is.finite(z[[paste0(p, "_estimate")]]) & is.finite(z[[paste0(p, "_se")]]) &
      is.finite(z[[paste0(p, "_low")]]) & is.finite(z[[paste0(p, "_high")]])
  }
  g_ok <- is.na(z$g_error) & fin("g")
  f_ok <- is.na(z$f_error) & fin("f")
  b_ok <- is.na(z$b_error) & fin("b")
  ok <- g_ok & f_ok & b_ok
  zo <- z[ok, ]
  n_ok <- nrow(zo)
  n_g_failed <- sum(!g_ok)
  n_b_failed <- sum(!b_ok)
  lim <- fail_limit(nrow(z))
  fits_ok <- cl$slope || (n_g_failed <= lim && n_b_failed <= lim)
  ratio <- (zo$b_high - zo$b_low) / (zo$g_high - zo$g_low)
  width_ratio <- stats::median(ratio)
  width_met <- isTRUE(width_ratio >= 0.8 && width_ratio <= 1.25)
  g_cov_ev <- cover(zo$truth, zo$g_low, zo$g_high)
  g_cov <- mean(g_cov_ev)
  g_mcse <- sqrt(g_cov * (1 - g_cov) / n_ok)
  ci <- wilson(sum(g_cov_ev), n_ok)
  cov_verdict <- band_verdict(g_cov, g_mcse, 0.93, 0.97)
  gate <- if (!cl$slope) "G3: coverage + width" else if (id == "A1W") "G4: coverage" else "G4: width"
  parts <- if (!cl$slope) {
    c(cov_verdict, if (width_met) "PASS" else "FAIL")
  } else if (id == "A1W") {
    cov_verdict
  } else {
    if (width_met) "PASS" else "FAIL"
  }
  verdict <- if (!fits_ok) {
    "FAIL (fits)"
  } else if (all(parts == "PASS")) {
    "PASS"
  } else if (any(parts == "FAIL")) {
    "FAIL"
  } else {
    "BORDERLINE"
  }
  f_cov <- mean(cover(zo$truth, zo$f_low, zo$f_high))
  b_cov <- mean(cover(zo$truth, zo$b_low, zo$b_high))
  zg <- (zo$g_estimate - zo$truth) / zo$g_se
  zf <- (zo$f_estimate - zo$truth) / zo$f_se
  zb <- (zo$b_estimate - zo$truth) / zo$b_se
  data.frame(
    cell = id,
    family = cl$family$family,
    slope = cl$slope,
    gate = gate,
    n_seeds = nrow(z),
    n_ok = n_ok,
    n_failed_glmmTMB = n_g_failed,
    n_failed_flexyBayes = n_b_failed,
    n_fallback = sum(z$g_fallback),
    n_warned_glmmTMB = sum(z$g_n_warn > 0L),
    truth_mean = mean(z$truth),
    coverage_joint = g_cov,
    coverage_joint_mcse = g_mcse,
    coverage_joint_lo = ci[1L],
    coverage_joint_hi = ci[2L],
    coverage_verdict = cov_verdict,
    median_width_ratio = width_ratio,
    width_ratio_q25 = stats::quantile(ratio, 0.25, names = FALSE),
    width_ratio_q75 = stats::quantile(ratio, 0.75, names = FALSE),
    verdict = verdict,
    misses_below_joint = sum(zo$truth < zo$g_low),
    misses_above_joint = sum(zo$truth > zo$g_high),
    sd_z_joint = stats::sd(zg),
    mean_z_joint = mean(zg),
    coverage_fixed = f_cov,
    sd_z_fixed = stats::sd(zf),
    median_se_ratio_fixed_over_joint = stats::median(zo$f_se / zo$g_se),
    coverage_flexyBayes = b_cov,
    sd_z_flexyBayes = stats::sd(zb),
    max_abs_estimate_gap_joint_fixed = max(abs(zo$g_estimate - zo$f_estimate)),
    n_near_coverage = sum(pmin(abs(zo$truth - zo$g_low), abs(zo$truth - zo$g_high)) < 0.15 * zo$g_se),
    median_secs_glmmTMB = stats::median(z$g_secs),
    median_secs_flexyBayes = stats::median(z$b_secs),
    stringsAsFactors = FALSE
  )
}))

all_met <- all(summary_tab$verdict == "PASS")

t_end <- Sys.time()
wall <- difftime(t_end, t_start, units = "mins")

# Write ------------------------------------------------------------------------

rows <- rows[order(rows$cell, rows$seed), ]
rownames(rows) <- NULL
utils::write.csv(summary_tab, file.path(out_dir, "planted_truth_S3.2_results.csv"), row.names = FALSE)
gz <- gzfile(file.path(out_dir, "planted_truth_S3.2_seeds.csv.gz"), "w")
utils::write.csv(rows, gz, row.names = FALSE)
close(gz)

sink(file.path(out_dir, "planted_truth_S3.2_session.txt"))
cat("Planted-truth suite S3.2 (rule: validation/recovery_rule_S3.2.md)\n")
cat("Started :", format(t_start, "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Finished:", format(t_end, "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Wall time (min):", round(as.numeric(wall), 2), "\n")
cat("Seed scale:", scale, "; cores:", cores, "; tasks:", length(tasks), "\n")
cat("Posterior draws per flexyBayes call:", n_draws, "\n")
cat("All gates met:", all_met, "\n")
cat("pollinatoR git HEAD:", system("git rev-parse HEAD", intern = TRUE), "\n")
cat("INLA num.threads:", INLA::inla.getOption("num.threads"),
    "; flexyBayes.inla_best_of:", getOption("flexyBayes.inla_best_of"), "\n")
cat("\nPackage versions:\n")
for (p in c("pollinatoR", "flexyBayes", "INLA", "glmmTMB", "TMB", "marginaleffects", "Matrix")) {
  cat(sprintf("  %-16s %s\n", p, as.character(utils::packageVersion(p))))
}
for (x in c("g_error", "f_error", "b_error")) {
  v <- rows[[x]]
  cat("\nErrors in", x, "(message, rows):\n")
  if (any(!is.na(v))) print(table(v[!is.na(v)])) else cat("  none\n")
}
for (x in c("g_first_warn", "b_first_warn")) {
  v <- rows[[x]]
  cat("\nWarnings in", x, "(first per row, rows; marginaleffects advisory excluded):\n")
  if (any(!is.na(v))) print(table(v[!is.na(v)])) else cat("  none\n")
}
cat("\n")
print(utils::sessionInfo())
sink()

print(summary_tab[c("cell", "gate", "n_ok", "n_failed_glmmTMB", "n_fallback",
                    "coverage_joint", "coverage_joint_mcse", "median_width_ratio",
                    "verdict", "coverage_fixed", "coverage_flexyBayes")], digits = 3)
message("All gates met: ", all_met)
message("Wall time: ", round(as.numeric(wall), 2), " min")
