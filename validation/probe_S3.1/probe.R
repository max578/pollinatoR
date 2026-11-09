# S3.1 coverage probe: flexyBayes posterior intervals on site-slope designs.
# Synthetic data only. Simulators copied from validation/planted_truth_S3.R (v0.3.0).
# Run from the package root (the probe ran on the v0.3.0 tag):
#   Rscript validation/probe_S3.1/probe.R <setting> <n_seeds> <cores> <out_csv>
args <- commandArgs(trailingOnly = TRUE)
setting <- args[1L]
n_seeds <- as.integer(args[2L])
cores <- as.integer(args[3L])
out_csv <- args[4L]

suppressMessages(pkgload::load_all(".", quiet = TRUE))
options(flexyBayes.inla_best_of = 1L, flexyBayes.silence_default_prior_note = TRUE)
INLA::inla.setOption(num.threads = "1:1")

rnb <- function(mu, theta) stats::rnbinom(length(mu), mu = mu, size = theta)
rbern <- function(p) stats::rbinom(length(p), 1L, p)
sim_multi_site <- function(n_site, n_rep, arms, beta, sd_site, sd_cell, draw, linkinv) {
  sites <- sprintf("s%02d", seq_len(n_site))
  g <- expand.grid(rep = seq_len(n_rep), trt = arms, site = sites, stringsAsFactors = FALSE)
  u_site <- stats::rnorm(n_site, 0, sd_site); names(u_site) <- sites
  u_cell <- matrix(stats::rnorm(n_site * length(arms), 0, sd_cell), n_site,
                   dimnames = list(sites, arms))
  eta <- sapply(arms, function(a) beta[[a]] + u_site[g$site] + u_cell[cbind(g$site, a)])
  own <- eta[cbind(seq_len(nrow(g)), match(g$trt, arms))]
  g$y <- draw(linkinv(own))
  list(data = g, eta = eta)
}
sim_rcbd <- function(n_block, n_rep, arms, beta, sd_block, draw, linkinv) {
  blocks <- sprintf("b%02d", seq_len(n_block))
  g <- expand.grid(rep = seq_len(n_rep), trt = arms, block = blocks, stringsAsFactors = FALSE)
  u_block <- stats::rnorm(n_block, 0, sd_block); names(u_block) <- blocks
  eta <- sapply(arms, function(a) beta[[a]] + u_block[g$block])
  own <- eta[cbind(seq_len(nrow(g)), match(g$trt, arms))]
  g$y <- draw(linkinv(own))
  list(data = g, eta = eta)
}
two <- c("Cage", "Open")
arms_two <- c(excluded = "Cage", open = "Open")

# Settings: design, number of sites, extra pollin_fit() args, INLA injection
cell_prior <- function(spec) {
  # Names of the engine route's random terms for unit "site", treatment "trt"
  eval(bquote(flexyBayes::fb_prior(
    sd(group = "site") ~ .(spec),
    sd(group = "pollin_site_trt") ~ .(spec)
  )))
}
S <- list(
  default     = list(design = "A1", n_site = 12L),
  eb          = list(design = "A1", n_site = 12L, inject = "eb"),
  finegrid    = list(design = "A1", n_site = 12L, inject = "finegrid"),
  hn1         = list(design = "A1", n_site = 12L, dots = list(prior = quote(cell_prior(quote(half_normal(scale = 1)))))),
  pc          = list(design = "A1", n_site = 12L, dots = list(prior = quote(cell_prior(quote(pc(upper = 1, prob = 0.01)))))),
  lognorm1    = list(design = "A1", n_site = 12L, dots = list(prior_vc_sd = 1)),
  sites24     = list(design = "A1", n_site = 24L),
  sites48     = list(design = "A1", n_site = 48L),
  A4_default  = list(design = "A4"),
  A4_finegrid = list(design = "A4", inject = "finegrid"),
  A4_lognorm1 = list(design = "A4", dots = list(prior_vc_sd = 1))
)
st <- S[[setting]]
stopifnot(!is.null(st))

if (!is.null(st$inject)) {
  # Not reachable through flexybayes()'s public arguments: wraps the single
  # INLA call site of flexyBayes, in this probe process only.
  strat <- st$inject
  ns <- asNamespace("flexyBayes")
  orig <- get(".inla_call", envir = ns)
  utils::assignInNamespace(".inla_call", function(args) {
    args$control.inla <- if (strat == "finegrid") list(int.strategy = "grid", dz = 0.25, diff.logdens = 10) else list(int.strategy = strat)
    orig(args)
  }, "flexyBayes")
}
dots <- lapply(st$dots, eval)

simulate <- function() {
  if (st$design == "A1") {
    sim_multi_site(st$n_site, 8L, two, log(c(Cage = 20, Open = 30)), 0.40, 0.20,
                   function(mu) rnb(mu, 8), exp)
  } else {
    sim_rcbd(8L, 24L, two, stats::qlogis(c(Cage = 0.55, Open = 0.75)), 0.40,
             rbern, stats::plogis)
  }
}
unit <- if (st$design == "A1") "site" else "block"
layout <- if (st$design == "A1") "multi_site" else "rcbd"
fam <- if (st$design == "A1") glmmTMB::nbinom2() else stats::binomial()
linkinv <- if (st$design == "A1") exp else stats::plogis

# Probe seeds: disjoint from the S3 seeds (3,000,000 + ...); data seed
# 900000 + s (A1 designs) / 950000 + s (A4); draw seed 990000 + s.
run_seed <- function(s) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed((if (st$design == "A1") 900000L else 950000L) + s)
  sim <- simulate()
  h <- linkinv(sim$eta[, "Open"]); l <- linkinv(sim$eta[, "Cage"])
  truth <- 1 - mean(l) / mean(h)
  t0 <- Sys.time()
  out <- tryCatch({
    pd <- pollin_data(sim$data, "y", "trt", unit = unit)
    ds <- pollin_design(pd, arms_two, layout)
    est <- pollin_estimand(ds, "dependence", "ratio", "loss")
    fit <- suppressWarnings(do.call(pollin_fit, c(list(est, engine = "flexyBayes",
                                                     family = fam, aggregate = FALSE), dots)))
    m <- fit$fits[[1L]]$inla
    tab <- pollin_effects(fit, interval = "posterior", n_draws = 1000L,
                          seed = 990000L + s, conf_level = 0.95)$table
    hp <- m$summary.hyperpar
    prec_cell <- grep("pollin_site_trt", rownames(hp))
    list(est = tab$estimate[1L], low = tab$conf.low[1L], high = tab$conf.high[1L],
         sd = tab$std.error[1L], nconfig = m$misc$configs$nconfig,
         strat = as.character(m$.args$control.inla$int.strategy)[1L],
         sd_cell_med = if (length(prec_cell)) 1 / sqrt(hp[prec_cell, "0.5quant"]) else NA_real_,
         err = NA_character_)
  }, error = function(e) list(est = NA, low = NA, high = NA, sd = NA, nconfig = NA,
                               strat = NA, sd_cell_med = NA, err = conditionMessage(e)))
  data.frame(setting = setting, seed = s, truth = truth, out,
             secs = as.numeric(difftime(Sys.time(), t0, units = "secs")),
             stringsAsFactors = FALSE)
}

t0 <- Sys.time()
res <- parallel::mclapply(seq_len(n_seeds), run_seed, mc.cores = cores, mc.preschedule = FALSE)
bad <- vapply(res, inherits, logical(1), what = "try-error")
if (any(bad)) stop(sum(bad), " crashed: ", as.character(res[[which(bad)[1L]]]))
res <- do.call(rbind, res)
utils::write.csv(res, out_csv, row.names = FALSE)
ok <- is.na(res$err)
cov <- mean(res$truth[ok] >= res$low[ok] & res$truth[ok] <= res$high[ok])
message(sprintf("%s: n_ok=%d/%d cov=%.3f (MCSE %.3f) medwidth=%.4f strat=%s nconfig=%s wall=%.1f min",
                setting, sum(ok), nrow(res), cov, sqrt(cov * (1 - cov) / sum(ok)),
                median(res$high[ok] - res$low[ok]), paste(unique(res$strat), collapse = ","),
                paste(range(res$nconfig, na.rm = TRUE), collapse = "-"),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
if (any(!ok)) print(table(res$err[!ok]))
