# Planted-truth recovery suite, Stage 3.1 (criteria C2 engine agreement and C3
# Bayesian calibration on the binomial-denominator cell B1 and the controls A2
# and A5).
#
# Implements validation/recovery_rule_S3.1.md exactly. Run from the package root:
#
#   Rscript --vanilla validation/planted_truth_S3.1.R [n_seeds_b1] [n_seeds_ctrl] [cores] [out_dir]
#
# Defaults: 200 seeds for B1, 100 for each of A2 and A5, 8 cores,
# out_dir = "validation". A smoke run uses smaller seed counts and a scratch
# out_dir; only the full run writes the results named in the rule.

t_start <- Sys.time()
args <- commandArgs(trailingOnly = TRUE)
n_seeds_b1 <- if (length(args) >= 1L) as.integer(args[1L]) else 200L
n_seeds_ctrl <- if (length(args) >= 2L) as.integer(args[2L]) else 100L
cores <- if (length(args) >= 3L) as.integer(args[3L]) else 8L
out_dir <- if (length(args) >= 4L) args[4L] else "validation"
stopifnot(n_seeds_b1 >= 0L, n_seeds_ctrl >= 0L, cores >= 1L, dir.exists(out_dir))

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

# Data designs: k is the seed block; n_seeds the design's seed count
data_designs <- list(
  B1 = list(
    k = 9L,
    n_seeds = n_seeds_b1,
    simulate = function() {
      sim_rcbd_trials(8L, 6L, two, stats::qlogis(c(Cage = 0.55, Open = 0.75)),
                      0.40, 20:40)
    },
    unit = "block", layout = "rcbd", arms = arms_two, denominator = "n"
  ),
  A2 = list(
    k = 2L,
    n_seeds = n_seeds_ctrl,
    simulate = function() {
      sim_rcbd(8L, 6L, two, log(c(Cage = 20, Open = 30)), 0.40,
               function(mu) stats::rpois(length(mu), mu), exp)
    },
    unit = "block", layout = "rcbd", arms = arms_two, denominator = NULL
  ),
  A5 = list(
    k = 5L,
    n_seeds = n_seeds_ctrl,
    simulate = function() {
      sim_rcbd(8L, 4L, three, log(c(Cage = 20, Open = 30, Hand = 36)), 0.40,
               function(mu) rnb(mu, 8), exp)
    },
    unit = "block", layout = "rcbd", arms = arms_three, denominator = NULL
  )
)

# Cells (rule Section 3)
cells <- list(
  B1 = list(data = "B1", family = stats::binomial(), linkinv = stats::plogis,
            channel = "dependence", scale = "ratio", basis = "loss"),
  A2 = list(data = "A2", family = stats::poisson(), linkinv = exp,
            channel = "dependence", scale = "ratio", basis = "loss"),
  A5 = list(data = "A5", family = glmmTMB::nbinom2(), linkinv = exp,
            channel = "deficit", scale = "ratio", basis = "gain")
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

# Truth: arm means at the realised random effects, averaged over the rows
# with equal weight (for B1, over plants, whatever their ovule count)
truth_of <- function(sim, hl, linkinv, scale, basis) {
  h <- linkinv(sim$eta[, hl[["high"]]])
  l <- linkinv(sim$eta[, hl[["low"]]])
  effect_of(mean(h), mean(l), scale, basis)
}

# Seeds (rule Section 5): disjoint from S1, S2, S3 and the S3.1 probe
seed_data <- function(k, s) 41000000L + 100000L * k + s
seed_draw <- function(k, s) 61000000L + 100000L * k + s

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

build_estimand <- function(d, dd, channel, scale, basis) {
  pd <- pollin_data(d, "y", "trt", unit = dd$unit, denominator = dd$denominator)
  ds <- pollin_design(pd, dd$arms, dd$layout)
  pollin_estimand(ds, channel, scale, basis)
}

# Task: one data set, its cell, both engines ------------------------------------

run_task <- function(task) {
  dd <- data_designs[[task$data]]
  s <- task$seed
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed_data(dd$k, s))
  sim <- dd$simulate()
  id <- task$data
  cl <- cells[[id]]
  hl <- contrast_of(dd$arms, cl$channel)
  truth <- truth_of(sim, hl, cl$linkinv, cl$scale, cl$basis)
  est <- build_estimand(sim$data, dd, cl$channel, cl$scale, cl$basis)
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
                   seed = seed_draw(dd$k, s), conf_level = conf_level)$table
  })
  t2 <- Sys.time()
  pick <- function(r, col) if (is.null(r$value)) NA_real_ else r$value[[col]][1L]
  data.frame(
    cell = id, seed = s,
    n_rows = nrow(sim$data),
    n_trials = if (is.null(dd$denominator)) NA_real_ else sum(sim$data$n),
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
}

# Tasks: B1 first (the longest), then the controls --------------------------------

tasks <- unlist(lapply(names(data_designs), function(k) {
  lapply(seq_len(data_designs[[k]]$n_seeds), function(s) list(data = k, seed = s))
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

# Summaries (rule Section 4) ---------------------------------------------------

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

summary_tab <- do.call(rbind, lapply(names(cells), function(id) {
  z <- rows[rows$cell == id, ]
  if (nrow(z) == 0L) {
    return(NULL)
  }
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
  width_met <- isTRUE(width_ratio >= 0.8 && width_ratio <= 1.25)
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
  g_cov <- mean(g_cov_ev)
  c3_cov <- band_verdict(b_cov, b_mcse, 0.93, 0.97)
  c3 <- if (!fits_ok) "FAIL" else c3_cov
  zb <- (zo$b_estimate - zo$truth) / zo$b_se
  data.frame(
    cell = id,
    family = cells[[id]]$family$family,
    channel = cells[[id]]$channel,
    scale = cells[[id]]$scale,
    basis = cells[[id]]$basis,
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
    width_ratio_q25 = stats::quantile(
      (zo$b_high - zo$b_low) / (zo$g_high - zo$g_low), 0.25, names = FALSE
    ),
    width_ratio_q75 = stats::quantile(
      (zo$b_high - zo$b_low) / (zo$g_high - zo$g_low), 0.75, names = FALSE
    ),
    c2_verdict = if (c2_pass) "PASS" else paste0("FAIL (", paste(c2_missed, collapse = ", "), ")"),
    coverage_flexyBayes = b_cov,
    coverage_mcse = b_mcse,
    coverage_lo = ci[1L],
    coverage_hi = ci[2L],
    c3_coverage_alone = c3_cov,
    c3_verdict = c3,
    misses_below_flexyBayes = sum(zo$truth < zo$b_low),
    misses_above_flexyBayes = sum(zo$truth > zo$b_high),
    sd_z_flexyBayes = stats::sd(zb),
    mean_z_flexyBayes = mean(zb),
    bias_flexyBayes = mean(zo$b_estimate - zo$truth),
    n_near_coverage = sum(pmin(abs(zo$truth - zo$b_low), abs(zo$truth - zo$b_high)) < 0.15 * zo$b_se),
    n_near_agreement = sum(abs(abs(zo$b_estimate - zo$g_estimate) - 0.25 * zo$g_se) < 0.15 * zo$b_se),
    coverage_glmmTMB = g_cov,
    coverage_glmmTMB_mcse = sqrt(g_cov * (1 - g_cov) / n_ok),
    bias_glmmTMB = mean(zo$g_estimate - zo$truth),
    median_secs_glmmTMB = stats::median(z$g_secs),
    median_secs_flexyBayes = stats::median(z$b_secs),
    median_trials = if (all(is.na(z$n_trials))) NA_real_ else stats::median(z$n_trials),
    stringsAsFactors = FALSE
  )
}))

c2_met <- all(summary_tab$c2_verdict == "PASS")
c3_met <- all(summary_tab$c3_verdict == "PASS")

t_end <- Sys.time()
wall <- difftime(t_end, t_start, units = "mins")

# Write ------------------------------------------------------------------------

rows <- rows[order(rows$cell, rows$seed), ]
rownames(rows) <- NULL
utils::write.csv(summary_tab, file.path(out_dir, "planted_truth_S3.1_results.csv"), row.names = FALSE)
gz <- gzfile(file.path(out_dir, "planted_truth_S3.1_seeds.csv.gz"), "w")
utils::write.csv(rows, gz, row.names = FALSE)
close(gz)

sink(file.path(out_dir, "planted_truth_S3.1_session.txt"))
cat("Planted-truth suite S3.1 (rule: validation/recovery_rule_S3.1.md)\n")
cat("Started :", format(t_start, "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Finished:", format(t_end, "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Wall time (min):", round(as.numeric(wall), 2), "\n")
cat(
  "Seeds B1:", n_seeds_b1, "; A2 and A5:", n_seeds_ctrl,
  "; cores:", cores, "; tasks:", length(tasks), "\n"
)
cat("Posterior draws per flexyBayes call:", n_draws, "\n")
cat("C2 met:", c2_met, "; C3 met:", c3_met, "\n")
cat("pollinatoR git HEAD:", system("git rev-parse HEAD", intern = TRUE), "\n")
cat("INLA num.threads:", INLA::inla.getOption("num.threads"),
    "; flexyBayes.inla_best_of:", getOption("flexyBayes.inla_best_of"), "\n")
cat("\nPackage versions:\n")
for (p in c("pollinatoR", "flexyBayes", "INLA", "glmmTMB", "TMB", "marginaleffects", "Matrix")) {
  cat(sprintf("  %-16s %s\n", p, as.character(utils::packageVersion(p))))
}
for (x in c("g_error", "b_error")) {
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

print(summary_tab[c("cell", "n_ok", "n_failed_glmmTMB", "n_failed_flexyBayes",
                    "agreement_share", "median_width_ratio", "c2_verdict",
                    "coverage_flexyBayes", "coverage_mcse", "c3_verdict",
                    "coverage_glmmTMB")], digits = 3)
message("C2 met: ", c2_met, "; C3 met: ", c3_met)
message("Wall time: ", round(as.numeric(wall), 2), " min")
