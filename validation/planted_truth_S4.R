# Planted-truth recovery suite, Stage 4 (identification, 0.5.0): gates G-M1,
# G-F1, G-F2a/b, G-F3, G-X1, G-X1b, G-X2a/b.
#
# Implements validation/recovery_rule_S4.md exactly. Run from the package root:
#
#   Rscript --vanilla validation/planted_truth_S4.R [scale] [cores] [out_dir]
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
zq <- stats::qnorm(0.975)
me_advisory <- "only takes into account the uncertainty in fixed-effect"

n_of <- function(n) max(1L, as.integer(round(scale * n)))
rnb <- function(mu) stats::rnbinom(length(mu), mu = mu, size = 10)
rgam <- function(mu) stats::rgamma(length(mu), shape = 10, rate = 10 / mu)

# Planted values (rule Section 1) -----------------------------------------------

m1_R <- exp(0.15) *
  ((1 - 0.012 * 40 / 15) / (1 - 0.012 * 40 * exp(0.25) / 15))^15
truth_m1_total <- 1 - 1 / m1_R
truth_m1_direct <- 1 - exp(-0.15)
truth_mod <- 0.10
truth_k <- 1.10
truth_net <- 0.25

# Simulators ---------------------------------------------------------------------

sim_m1 <- function(n_field) {
  fields <- sprintf("F%02d", seq_len(n_field))
  g <- expand.grid(
    plant = seq_len(5L),
    trt = c("Cage", "Open"),
    field = fields,
    stringsAsFactors = FALSE
  )
  u <- stats::rnorm(n_field, 0, 0.3)
  names(u) <- fields
  open <- g$trt == "Open"
  m <- ifelse(open, 40 * exp(0.25), 40)
  g$veg <- stats::rgamma(nrow(g), shape = 15, rate = 15 / m)
  g$y <- rnb(exp(3 + u[g$field] + 0.15 * open + 0.012 * g$veg))
  g
}

# The standardised between-unit component, formed as pollin_moderation() does
zb_of <- function(bees, site) .z(pollin_within_between(bees, site)$between)

# F1-F3 share one site layout; `kind` picks the planted structure
sim_field <- function(kind) {
  n_site <- 44L
  sites <- sprintf("s%02d", seq_len(n_site))
  g <- expand.grid(
    plant = seq_len(6L),
    trt = c("Cage", "Open"),
    site = sites,
    stringsAsFactors = FALSE
  )
  si <- match(g$site, sites)
  open <- g$trt == "Open"
  x <- stats::rnorm(n_site)
  b_site <- if (kind == "F3") {
    stats::rpois(n_site, exp(2))
  } else {
    stats::rpois(n_site, exp(2 + 0.5 * x))
  }
  g$bees <- b_site[si]
  if (kind != "F3") {
    g$bees <- g$bees + stats::rpois(nrow(g), 2)
  }
  u <- stats::rnorm(n_site, 0, 0.3)
  if (kind == "F1") {
    tau <- 0.3 + stats::rnorm(n_site, 0, 0.15)
    eta <- 3 + u[si] + 0.4 * x[si] + tau[si] * open
    b_z <- NA_real_
  } else {
    zb <- zb_of(g$bees, g$site)
    zb_site <- zb[match(sites, g$site)]
    tau <- 0.3 + truth_mod * zb_site + stats::rnorm(n_site, 0, 0.10)
    if (kind == "F2") {
      tau <- tau + 0.20 * x
      b_z <- unname(stats::coef(stats::lm(x ~ zb_site))[2L])
      eta <- 3 + u[si] + 0.3 * x[si] + tau[si] * open
    } else {
      b_z <- NA_real_
      eta <- 3 + u[si] + tau[si] * open
    }
  }
  g$z <- x[si]
  g$y <- rnb(exp(eta))
  list(data = g, b_z = b_z)
}

sim_x <- function(k, sham) {
  arms <- if (sham) c("Cage", "Sham", "Open") else c("Cage", "Open")
  blocks <- sprintf("b%02d", seq_len(10L))
  g <- expand.grid(
    plant = seq_len(8L),
    trt = arms,
    block = blocks,
    stringsAsFactors = FALSE
  )
  u <- stats::rnorm(10L, 0, 0.3)
  names(u) <- blocks
  mult <- c(Open = 1, Sham = k, Cage = 0.75 * k)[g$trt]
  g$y <- rgam(mult * exp(1 + u[g$block]))
  g
}

# Running one step, recording errors and warnings --------------------------------

run_step <- function(expr_fun) {
  warns <- character(0)
  res <- tryCatch(
    withCallingHandlers(
      suppressMessages(expr_fun()),
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
  list(
    value = if (is.na(err)) res else NULL,
    error = err,
    n_warn = length(warns),
    first_warn = if (length(warns)) warns[1L] else NA_character_
  )
}

glmm_ok <- function(m) {
  isTRUE(m$fit$convergence == 0L) && isTRUE(m$sdr$pdHess)
}

# One fitted interval as a row fragment, prefixed by `p`
fit_effect <- function(est, p, engine = "glmmTMB", family, seed = NULL) {
  t0 <- Sys.time()
  r <- run_step(function() {
    if (engine == "glmmTMB") {
      fit <- pollin_fit(est, engine = "glmmTMB", family = family)
      if (!glmm_ok(fit$fits[[1L]])) stop("glmmTMB: not converged")
      pollin_effects(fit, interval = "delta", conf_level = conf_level)$table
    } else {
      fit <- pollin_fit(est, engine = "flexyBayes", family = family,
                        aggregate = FALSE)
      pollin_effects(fit, interval = "posterior", n_draws = n_draws,
                     seed = seed, conf_level = conf_level)$table
    }
  })
  pick <- function(col) if (is.null(r$value)) NA else r$value[[col]][1L]
  out <- list(
    estimate = pick("estimate"), low = pick("conf.low"), high = pick("conf.high"),
    label = if (is.null(r$value)) NA_character_ else as.character(pick("effect")),
    error = r$error, n_warn = r$n_warn, first_warn = r$first_warn,
    secs = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  )
  stats::setNames(out, paste0(p, "_", names(out)))
}

# One moderation term as a row fragment
fit_moderation <- function(ds, p, dag = NULL) {
  t0 <- Sys.time()
  r <- run_step(function() {
    m <- pollin_moderation(ds, outcome = "y", family = glmmTMB::nbinom2(),
                           dag = dag)
    if (!glmm_ok(m$fit)) stop("glmmTMB: not converged")
    m
  })
  m <- r$value
  term <- function(tab, nm, col) {
    if (is.null(tab)) return(NA_real_)
    i <- grep(nm, tab$term, fixed = TRUE)
    if (length(i) != 1L) NA_real_ else tab[[col]][i]
  }
  tab <- if (is.null(m)) NULL else m$table
  nc <- if (is.null(m)) NULL else m$negative_control
  out <- list(
    b_est = term(tab, "bees_between", "estimate"),
    b_se = term(tab, "bees_between", "std.error"),
    b_p = term(tab, "bees_between", "p.value"),
    w_est = term(tab, "bees_within", "estimate"),
    w_se = term(tab, "bees_within", "std.error"),
    w_p = term(tab, "bees_within", "p.value"),
    ncb_p = term(nc, "bees_between", "p.value"),
    ncw_p = term(nc, "bees_within", "p.value"),
    between_only = if (is.null(m)) NA else "bees" %in% m$between_only,
    error = r$error, n_warn = r$n_warn, first_warn = r$first_warn,
    secs = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  )
  stats::setNames(out, paste0(p, "_", names(out)))
}

arms2 <- c(excluded = "Cage", open = "Open")
arms3 <- c(excluded = "Cage", open = "Open", sham = "Sham")
f2_dag <- "dag { z -> bees; z -> yield; bees -> yield; trt -> yield }"

# Cells (rule Section 1) ----------------------------------------------------------

cells <- list(
  M1a = list(k = 1L, n = n_of(1000L)),
  M1b = list(k = 2L, n = n_of(1000L)),
  F1 = list(k = 3L, n = n_of(1000L)),
  F2 = list(k = 4L, n = n_of(1000L)),
  F3 = list(k = 5L, n = n_of(1000L)),
  X1 = list(k = 6L, n = n_of(1000L)),
  X1B = list(k = 6L, n = n_of(200L)),
  X2 = list(k = 7L, n = n_of(1000L))
)

seed_data <- function(k, s) 81000000L + 100000L * k + s
seed_draw <- function(k, s) 91000000L + 100000L * k + s

run_task <- function(task) {
  id <- task$cell
  s <- task$seed
  k <- cells[[id]]$k
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed_data(k, s))
  row <- list(cell = id, seed = s)

  if (id %in% c("M1a", "M1b")) {
    d <- sim_m1(if (id == "M1a") 24L else 44L)
    ds_t <- pollin_design(pollin_data(d, "y", "trt", unit = "field"),
                          arms2, "matched_pairs")
    ds_d <- pollin_design(pollin_data(d, "y", "trt", unit = "field",
                                      adjust = "veg"),
                          arms2, "matched_pairs")
    est_t <- pollin_estimand(ds_t, "dependence", "ratio", "loss")
    est_d <- pollin_estimand(ds_d, "dependence", "ratio", "loss",
                             effect = "direct", mediators = "veg")
    row <- c(row, list(truth_a = truth_m1_total, truth_b = truth_m1_direct),
             fit_effect(est_t, "a", family = glmmTMB::nbinom2()),
             fit_effect(est_d, "b", family = glmmTMB::nbinom2()))
  } else if (id %in% c("F1", "F2", "F3")) {
    sim <- sim_field(id)
    d <- sim$data
    ds0 <- pollin_design(pollin_data(d, "y", "trt", unit = "site",
                                     exposures = "bees"),
                         arms2, "multi_site")
    row <- c(row, list(b_z = sim$b_z), fit_moderation(ds0, "u"))
    if (id == "F2") {
      ds1 <- pollin_design(pollin_data(d, "y", "trt", unit = "site",
                                       exposures = "bees", adjust = "z"),
                           arms2, "multi_site")
      row <- c(row, fit_moderation(ds1, "g", dag = f2_dag_obj))
    }
  } else if (id %in% c("X1", "X1B")) {
    d <- sim_x(truth_k, sham = TRUE)
    ds <- pollin_design(pollin_data(d, "y", "trt", unit = "block"),
                        arms3, "rcbd")
    est_c <- pollin_estimand(ds, "cage", "ratio")
    fam <- stats::Gamma(link = "log")
    if (id == "X1") {
      ds_s <- pollin_design(pollin_data(d, "y", "trt", unit = "block"),
                            arms3, "rcbd", reference = "sham")
      est_n <- pollin_estimand(ds_s, "dependence", "ratio", "loss")
      row <- c(row, list(truth_a = truth_k, truth_b = truth_net),
               fit_effect(est_c, "a", family = fam),
               fit_effect(est_n, "b", family = fam))
    } else {
      row <- c(row, list(truth_a = truth_k),
               fit_effect(est_c, "a", "flexyBayes", fam, seed_draw(k, s)),
               fit_effect(est_c, "g", family = fam))
    }
  } else if (id == "X2") {
    kk <- stats::runif(1L, 0.95, 1.25)
    d <- sim_x(kk, sham = FALSE)
    ds <- pollin_design(pollin_data(d, "y", "trt", unit = "block"),
                        arms2, "rcbd")
    est <- pollin_estimand(ds, "dependence", "ratio", "loss")
    t0 <- Sys.time()
    r <- run_step(function() {
      fit <- pollin_fit(est, engine = "glmmTMB", family = stats::Gamma(link = "log"))
      if (!glmm_ok(fit$fits[[1L]])) stop("glmmTMB: not converged")
      dr <- pollin_draws(fit, n = n_draws, seed = seed_draw(k, s))
      b <- pollin_bound(dr, cage_factor = c(1, 1.2))$table
      # G-X2a: the same call on the planted shares
      dp <- dr
      dp$low[] <- 0.75 * kk
      dp$high[] <- 1
      bp <- pollin_bound(dp, cage_factor = c(1, 1.2))$table
      list(b = b, bp = bp)
    })
    v <- r$value
    row <- c(row, list(
      k_planted = kk,
      inside = kk >= 1 && kk <= 1.2,
      x_observed = if (is.null(v)) NA else v$b$observed,
      x_lower = if (is.null(v)) NA else v$b$lower,
      x_upper = if (is.null(v)) NA else v$b$upper,
      x_low = if (is.null(v)) NA else v$b$conf.low,
      x_high = if (is.null(v)) NA else v$b$conf.high,
      x_alg_lower = if (is.null(v)) NA else v$bp$lower,
      x_alg_upper = if (is.null(v)) NA else v$bp$upper,
      x_error = r$error, x_n_warn = r$n_warn, x_first_warn = r$first_warn,
      x_secs = as.numeric(difftime(Sys.time(), t0, units = "secs"))
    ))
  }
  row
}

# dagitty runs on V8, which macOS kills in a forked child, so the F2 graph is
# built once here. Its sets depend only on the data roles, the same in every seed.
f2_dag_obj <- local({
  set.seed(1L)
  d <- sim_field("F2")$data
  pollin_dag(
    pollin_design(pollin_data(d, "y", "trt", unit = "site",
                              exposures = "bees", adjust = "z"),
                  arms2, "multi_site"),
    f2_dag, "yield"
  )
})
stopifnot(identical(f2_dag_obj$exposure_sets$bees[[1L]], "z"))

tasks <- unlist(lapply(names(cells), function(id) {
  lapply(seq_len(cells[[id]]$n), function(s) list(cell = id, seed = s))
}), recursive = FALSE)

message("Running ", length(tasks), " tasks on ", cores, " cores ...")
res <- parallel::mclapply(tasks, run_task, mc.cores = cores,
                          mc.preschedule = FALSE)
bad <- vapply(res, function(r) is.null(r) || inherits(r, "try-error"), logical(1))
if (any(bad)) {
  stop(sum(bad), " task(s) crashed outside the cell handler; first: ",
       as.character(res[[which(bad)[1L]]]))
}

# Bind rows with differing fields
all_cols <- unique(unlist(lapply(res, names)))
rows <- do.call(rbind, lapply(res, function(r) {
  r[setdiff(all_cols, names(r))] <- NA
  as.data.frame(r[all_cols], stringsAsFactors = FALSE)
}))

# Gates (rule Section 2) ------------------------------------------------------------

band_verdict <- function(x, mcse, lo, hi) {
  if (!is.finite(x) || !is.finite(mcse)) return("FAIL")
  if (x >= lo && x <= hi) return("PASS")
  gap <- if (x < lo) lo - x else x - hi
  if (gap <= 2 * mcse) "BORDERLINE" else "FAIL"
}
prop_mcse <- function(p, n) sqrt(p * (1 - p) / n)
cover <- function(t, lo, hi) t >= lo & t <= hi

gates <- list()
add_gate <- function(gate, cell, quantity, value, mcse, requirement, verdict,
                     n_ok, n_failed, note = "") {
  gates[[length(gates) + 1L]] <<- data.frame(
    gate = gate, cell = cell, quantity = quantity, value = value, mcse = mcse,
    requirement = requirement, verdict = verdict, n_ok = n_ok,
    n_failed = n_failed, note = note, stringsAsFactors = FALSE
  )
}
fits_verdict <- function(n_failed, n) if (n_failed <= floor(0.02 * n)) "PASS" else "FAIL (fits)"

for (id in c("M1a", "M1b", "X1")) {
  z <- rows[rows$cell == id, ]
  if (!nrow(z)) next
  ok <- is.na(z$a_error) & is.na(z$b_error) & is.finite(z$a_low) & is.finite(z$b_low)
  zo <- z[ok, ]
  fv <- fits_verdict(sum(!ok), nrow(z))
  gname <- if (id == "X1") "G-X1" else "G-M1"
  qa <- if (id == "X1") "cage factor coverage" else "total coverage"
  qb <- if (id == "X1") "net dependence coverage" else "direct coverage"
  for (p in c("a", "b")) {
    cv <- mean(cover(zo[[paste0("truth_", p)]], zo[[paste0(p, "_low")]],
                     zo[[paste0(p, "_high")]]))
    se <- prop_mcse(cv, nrow(zo))
    v <- band_verdict(cv, se, 0.93, 0.97)
    if (fv != "PASS") v <- fv
    add_gate(gname, id, if (p == "a") qa else qb, cv, se, "[0.93, 0.97]", v,
             nrow(zo), sum(!ok),
             sprintf("bias %.4f; median width %.4f; below %d, above %d",
                     mean(zo[[paste0(p, "_estimate")]] - zo[[paste0("truth_", p)]]),
                     stats::median(zo[[paste0(p, "_high")]] - zo[[paste0(p, "_low")]]),
                     sum(zo[[paste0("truth_", p)]] < zo[[paste0(p, "_low")]]),
                     sum(zo[[paste0("truth_", p)]] > zo[[paste0(p, "_high")]])))
  }
  if (id != "X1") {
    lab_ok <- all(zo$a_label == "total") && all(zo$b_label == "direct")
    add_gate("G-M1", id, "labels as fitted", as.numeric(lab_ok), NA,
             "all seeds", if (lab_ok) "PASS" else "FAIL", nrow(zo), sum(!ok),
             sprintf("mean total %.4f, mean direct %.4f",
                     mean(zo$a_estimate), mean(zo$b_estimate)))
  }
}

z <- rows[rows$cell == "F1", ]
if (nrow(z)) {
  ok <- is.na(z$u_error) & is.finite(z$u_b_p) & is.finite(z$u_w_p)
  zo <- z[ok, ]
  n <- nrow(zo)
  one_sided <- function(x, se, bound, upper) {
    if (!is.finite(x)) return("FAIL")
    if (if (upper) x <= bound else x >= bound) return("PASS")
    if (abs(x - bound) <= 2 * se) "BORDERLINE" else "FAIL"
  }
  for (t in c("w", "b")) {
    fc <- mean(zo[[paste0("u_", t, "_p")]] < 0.05)
    se <- prop_mcse(max(fc, 0.05), n)
    add_gate("G-F1", "F1", paste("false-call", if (t == "w") "within" else "between"),
             fc, se, "<= 0.07", one_sided(fc, se, 0.07, TRUE), n, sum(!ok),
             "failure share reported, not gated")
  }
  nc <- mean(zo$u_ncb_p < 0.05)
  se <- prop_mcse(nc, n)
  add_gate("G-F1", "F1", "negative control rejects (between)", nc, se, ">= 0.80",
           one_sided(nc, se, 0.80, FALSE), n, sum(!ok),
           sprintf("within negative control rejects %.3f", mean(zo$u_ncw_p < 0.05)))
}

z <- rows[rows$cell == "F2", ]
if (nrow(z)) {
  ok <- is.na(z$g_error) & is.finite(z$g_b_est) & is.finite(z$g_b_se)
  zo <- z[ok, ]
  cv <- mean(cover(truth_mod, zo$g_b_est - zq * zo$g_b_se, zo$g_b_est + zq * zo$g_b_se))
  se <- prop_mcse(cv, nrow(zo))
  add_gate("G-F2a", "F2", "moderation coverage with the graph's set", cv, se,
           "[0.93, 0.97]", band_verdict(cv, se, 0.93, 0.97), nrow(zo), sum(!ok),
           sprintf("bias %.4f", mean(zo$g_b_est - truth_mod)))
  ok <- is.na(z$u_error) & is.finite(z$u_b_est) & is.finite(z$b_z)
  zo <- z[ok, ]
  n <- nrow(zo)
  dev <- zo$u_b_est - truth_mod
  gap <- mean(dev) - mean(0.20 * zo$b_z)
  se_gap <- stats::sd(dev - 0.20 * zo$b_z) / sqrt(n)
  se_dev <- stats::sd(dev) / sqrt(n)
  v1 <- abs(gap) <= 2 * se_gap
  v2 <- abs(mean(dev)) > 2 * se_dev
  add_gate("G-F2b", "F2", "unadjusted bias minus planted bias", gap, se_gap,
           "within 2 MC SE of 0; bias itself > 2 MC SE",
           if (v1 && v2) "PASS" else "FAIL", n, sum(!ok),
           sprintf("mean bias %.4f (MC SE %.4f); mean planted bias %.4f; mean b_z %.4f",
                   mean(dev), se_dev, mean(0.20 * zo$b_z), mean(zo$b_z)))
}

z <- rows[rows$cell == "F3", ]
if (nrow(z)) {
  conv <- is.na(z$u_error)
  zo <- z[conv, ]
  rank_ok <- all(is.finite(zo$u_b_est) & is.finite(zo$u_b_se)) &&
    all(zo$u_between_only %in% TRUE)
  cv <- mean(cover(truth_mod, zo$u_b_est - zq * zo$u_b_se, zo$u_b_est + zq * zo$u_b_se))
  se <- prop_mcse(cv, nrow(zo))
  v <- band_verdict(cv, se, 0.93, 0.97)
  if (!rank_ok) v <- "FAIL (rank)"
  add_gate("G-F3", "F3", "between-only term coverage", cv, se,
           "[0.93, 0.97]; no missing term; bees between only", v, nrow(zo),
           sum(!conv), sprintf("bias %.4f; rank and label check %s",
                               mean(zo$u_b_est - truth_mod), rank_ok))
}

z <- rows[rows$cell == "X1B", ]
if (nrow(z)) {
  ok <- is.na(z$a_error) & is.finite(z$a_low)
  share <- mean(ok)
  okg <- ok & is.na(z$g_error) & is.finite(z$g_low)
  zo <- z[okg, ]
  cv <- mean(cover(truth_k, zo$a_low, zo$a_high))
  wr <- stats::median((zo$a_high - zo$a_low) / (zo$g_high - zo$g_low))
  add_gate("G-X1b", "X1B", "flexyBayes cage rows returned", share,
           prop_mcse(share, nrow(z)), ">= 0.98",
           if (share >= 0.98) "PASS" else "FAIL", sum(ok), sum(!ok),
           sprintf("posterior coverage %.3f; median width ratio fB/glmmTMB %.3f",
                   cv, wr))
}

z <- rows[rows$cell == "X2", ]
if (nrow(z)) {
  ok <- is.na(z$x_error) & is.finite(z$x_low)
  zo <- z[ok, ]
  ins <- zo[zo$inside, ]
  outs <- zo[!zo$inside, ]
  alg <- all(ins$x_alg_lower <= truth_net + 1e-12 & truth_net <= ins$x_alg_upper + 1e-12)
  add_gate("G-X2a", "X2", "algebraic bound contains 0.25 (inside seeds)",
           mean(ins$x_alg_lower <= truth_net + 1e-12 & truth_net <= ins$x_alg_upper + 1e-12),
           NA, "1 (every inside seed)",
           if (alg) "PASS" else "FAIL", nrow(ins), sum(!ok))
  cv <- mean(cover(truth_net, ins$x_low, ins$x_high))
  se <- prop_mcse(cv, nrow(ins))
  v <- if (cv >= 0.95) "PASS" else if (0.95 - cv <= 2 * se) "BORDERLINE" else "FAIL"
  if (fits_verdict(sum(!ok), nrow(z)) != "PASS") v <- "FAIL (fits)"
  add_gate("G-X2b", "X2", "outer interval coverage (inside seeds)", cv, se,
           ">= 0.95", v, nrow(ins), sum(!ok),
           sprintf(paste("outside seeds n %d coverage %.3f; median bounds contain",
                         "0.25: inside %.3f, outside %.3f"),
                   nrow(outs), mean(cover(truth_net, outs$x_low, outs$x_high)),
                   mean(cover(truth_net, ins$x_lower, ins$x_upper)),
                   mean(cover(truth_net, outs$x_lower, outs$x_upper))))
}

gate_tab <- do.call(rbind, gates)
all_met <- all(gate_tab$verdict == "PASS")
t_end <- Sys.time()
wall <- difftime(t_end, t_start, units = "mins")

# Write ---------------------------------------------------------------------------

rows <- rows[order(rows$cell, rows$seed), ]
rownames(rows) <- NULL
utils::write.csv(gate_tab, file.path(out_dir, "planted_truth_S4_results.csv"),
                 row.names = FALSE)
gz <- gzfile(file.path(out_dir, "planted_truth_S4_seeds.csv.gz"), "w")
utils::write.csv(rows, gz, row.names = FALSE)
close(gz)

sink(file.path(out_dir, "planted_truth_S4_session.txt"))
cat("Planted-truth suite S4 (rule: validation/recovery_rule_S4.md)\n")
cat("Started :", format(t_start, "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Finished:", format(t_end, "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Wall time (min):", round(as.numeric(wall), 2), "\n")
cat("Seed scale:", scale, "; cores:", cores, "; tasks:", length(tasks), "\n")
cat("Planted M1 total:", sprintf("%.6f", truth_m1_total),
    "; direct:", sprintf("%.6f", truth_m1_direct), "\n")
cat("All gates met:", all_met, "\n")
cat("pollinatoR git HEAD:", system("git rev-parse HEAD", intern = TRUE), "\n")
cat("\nPackage versions:\n")
for (p in c("pollinatoR", "flexyBayes", "INLA", "glmmTMB", "TMB", "dagitty", "Matrix")) {
  cat(sprintf("  %-16s %s\n", p, as.character(utils::packageVersion(p))))
}
for (x in grep("_error$", names(rows), value = TRUE)) {
  v <- rows[[x]]
  cat("\nErrors in", x, "(message, rows):\n")
  if (any(!is.na(v))) print(table(v[!is.na(v)])) else cat("  none\n")
}
for (x in grep("_first_warn$", names(rows), value = TRUE)) {
  v <- rows[[x]]
  cat("\nWarnings in", x, "(first per row, rows):\n")
  if (any(!is.na(v))) print(table(v[!is.na(v)])) else cat("  none\n")
}
cat("\nMedian seconds per step:\n")
for (x in grep("_secs$", names(rows), value = TRUE)) {
  cat(sprintf("  %-10s %.2f\n", x, stats::median(rows[[x]], na.rm = TRUE)))
}
cat("\n")
print(utils::sessionInfo())
sink()

print(gate_tab[c("gate", "cell", "quantity", "value", "mcse", "verdict", "n_ok",
                 "n_failed")], digits = 3)
message("All gates met: ", all_met)
message("Wall time: ", round(as.numeric(wall), 2), " min")
