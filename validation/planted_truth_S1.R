# Planted-truth recovery suite, Stage 1 (criterion C2).
#
# Implements validation/recovery_rule_S1.md exactly. Run from the package root:
#
#   Rscript --vanilla validation/planted_truth_S1.R [n_seeds] [cores] [out_dir]
#
# Defaults: 1000 seeds, 8 cores, out_dir = "validation". A smoke run uses a
# smaller n_seeds and a scratch out_dir; only the full run writes the results
# named in the rule.

t_start <- Sys.time()
args <- commandArgs(trailingOnly = TRUE)
n_seeds <- if (length(args) >= 1L) as.integer(args[1L]) else 1000L
cores <- if (length(args) >= 2L) as.integer(args[2L]) else 8L
out_dir <- if (length(args) >= 3L) args[3L] else "validation"
stopifnot(n_seeds >= 1L, cores >= 1L, dir.exists(out_dir))

suppressMessages(devtools::load_all(".", quiet = TRUE))

# Designs (rule Section 2) -----------------------------------------------------

designs <- list(
  D1 = list(
    id = 1L,
    name = "three-arm trial",
    layout = "rcbd",
    unit = "block",
    arms = c(excluded = "Cage", open = "Open", supplemented = "Hand"),
    mean = c(Cage = 20, Open = 30, Hand = 36),
    theta = 8,
    families = c("deficit", "dependence"),
    simulate = function(d) {
      g <- expand.grid(
        plant = 1:4,
        trt = c("Cage", "Open", "Hand"),
        block = sprintf("b%02d", 1:8),
        stringsAsFactors = FALSE
      )
      u_block <- stats::rnorm(8L, 0, 0.40)
      names(u_block) <- sprintf("b%02d", 1:8)
      list(data = g, u = unname(u_block[g$block]))
    }
  ),
  D2 = list(
    id = 2L,
    name = "split-plot trial",
    layout = "split_plot",
    unit = c("block", "wplot"),
    arms = c(excluded = "Cage", open = "Open", supplemented = "Hand"),
    mean = c(Cage = 50, Open = 70, Hand = 77),
    theta = 15,
    families = c("deficit", "dependence"),
    simulate = function(d) {
      g <- expand.grid(
        trt = c("Cage", "Open", "Hand"),
        wplot = sprintf("w%d", 1:4),
        block = sprintf("b%d", 1:6),
        stringsAsFactors = FALSE
      )
      u_block <- stats::rnorm(6L, 0, 0.30)
      names(u_block) <- sprintf("b%d", 1:6)
      wp_key <- paste(rep(sprintf("b%d", 1:6), each = 4L), sprintf("w%d", 1:4))
      u_wp <- stats::rnorm(24L, 0, 0.25)
      names(u_wp) <- wp_key
      u <- u_block[g$block] + u_wp[paste(g$block, g$wplot)]
      list(data = g, u = unname(u))
    }
  ),
  D3 = list(
    id = 3L,
    name = "matched-pair field study",
    layout = "matched_pairs",
    unit = c("site", "pair"),
    arms = c(excluded = "Cage", open = "Open"),
    mean = c(Cage = 150, Open = 187.5),
    theta = 5,
    families = "dependence",
    simulate = function(d) {
      g <- expand.grid(
        trt = c("Cage", "Open"),
        pair = sprintf("p%d", 1:5),
        site = sprintf("s%02d", 1:20),
        stringsAsFactors = FALSE
      )
      u_site <- stats::rnorm(20L, 0, 0.50)
      names(u_site) <- sprintf("s%02d", 1:20)
      pr_key <- paste(rep(sprintf("s%02d", 1:20), each = 5L), sprintf("p%d", 1:5))
      u_pair <- stats::rnorm(100L, 0, 0.30)
      names(u_pair) <- pr_key
      u <- u_site[g$site] + u_pair[paste(g$site, g$pair)]
      list(data = g, u = unname(u))
    }
  )
)

# Scale x basis combinations (rule Section 3) ---------------------------------
scales <- data.frame(
  scale = c("difference", "ratio", "ratio", "rii", "lnrr"),
  basis = c("loss", "loss", "gain", "loss", "loss"),
  stringsAsFactors = FALSE
)

# Simulation and truth (rule Sections 1, 4) ------------------------------------

simulate_design <- function(des, s) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(100000L * des$id + s)
  sim <- des$simulate(des)
  d <- sim$data
  mu <- des$mean[d$trt] * exp(sim$u)
  d$y <- stats::rnbinom(nrow(d), mu = unname(mu), size = des$theta)
  list(data = d, u = sim$u)
}

truth_value <- function(des, u, family, scale, basis) {
  levels_hl <- if (family == "deficit") {
    c(des$arms[["supplemented"]], des$arms[["open"]])
  } else {
    c(des$arms[["open"]], des$arms[["excluded"]])
  }
  H <- mean(des$mean[[levels_hl[1L]]] * exp(u))
  L <- mean(des$mean[[levels_hl[2L]]] * exp(u))
  switch(
    scale,
    difference = H - L,
    ratio = if (basis == "loss") 1 - L / H else H / L - 1,
    rii = (H - L) / (H + L),
    lnrr = log(H / L)
  )
}

# One cell on one simulated data set (public API only) -------------------------

me_advisory <- "only takes into account the uncertainty in fixed-effect"

run_cell <- function(des, sim, family, scale, basis) {
  warns <- character(0)
  res <- tryCatch(
    withCallingHandlers(
      {
        pd <- pollin_data(sim$data, "y", "trt", unit = des$unit)
        ds <- pollin_design(pd, des$arms, des$layout)
        est <- pollin_estimand(ds, family, scale, basis, population = "trial")
        fit <- pollin_fit(est, engine = "glmmTMB", family = glmmTMB::nbinom2())
        pollin_effects(fit, interval = "delta", conf_level = 0.95)$table
      },
      warning = function(w) {
        msg <- conditionMessage(w)
        if (!grepl(me_advisory, msg, fixed = TRUE)) {
          warns <<- c(warns, msg)
        }
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) conditionMessage(e)
  )
  out <- data.frame(
    estimate = NA_real_,
    std.error = NA_real_,
    conf.low = NA_real_,
    conf.high = NA_real_,
    error = NA_character_,
    n_warn = length(warns),
    first_warn = if (length(warns)) warns[1L] else NA_character_,
    stringsAsFactors = FALSE
  )
  if (is.character(res)) {
    out$error <- res
  } else {
    out[c("estimate", "std.error", "conf.low", "conf.high")] <-
      res[1L, c("estimate", "std.error", "conf.low", "conf.high")]
  }
  out
}

# Tasks: one per design x seed --------------------------------------------------

run_task <- function(task) {
  des <- designs[[task$design]]
  sim <- simulate_design(des, task$seed)
  rows <- list()
  for (fam in des$families) {
    for (j in seq_len(nrow(scales))) {
      sc <- scales$scale[j]
      bs <- scales$basis[j]
      r <- run_cell(des, sim, fam, sc, bs)
      rows[[length(rows) + 1L]] <- cbind(
        data.frame(
          design = task$design,
          seed = task$seed,
          family = fam,
          scale = sc,
          basis = bs,
          truth = truth_value(des, sim$u, fam, sc, bs),
          stringsAsFactors = FALSE
        ),
        r
      )
    }
  }
  do.call(rbind, rows)
}

tasks <- unlist(
  lapply(names(designs), function(k) {
    lapply(seq_len(n_seeds), function(s) list(design = k, seed = s))
  }),
  recursive = FALSE
)

message(
  "Running ", length(tasks), " design x seed tasks on ", cores, " cores ..."
)
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
seeds <- do.call(rbind, res)
seeds <- seeds[order(seeds$design, seeds$family, seeds$scale, seeds$basis, seeds$seed), ]
rownames(seeds) <- NULL

# Summary per cell (rule Section 6) --------------------------------------------

wilson <- function(x, n, z = stats::qnorm(0.975)) {
  p <- x / n
  centre <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  c(centre - half, centre + half)
}

cells <- unique(seeds[c("design", "family", "scale", "basis")])
summ <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
  key <- cells[i, ]
  z <- seeds[
    seeds$design == key$design & seeds$family == key$family &
      seeds$scale == key$scale & seeds$basis == key$basis,
  ]
  ok <- is.na(z$error) &
    is.finite(z$estimate) & is.finite(z$std.error) &
    is.finite(z$conf.low) & is.finite(z$conf.high)
  zo <- z[ok, ]
  n_ok <- nrow(zo)
  covered <- zo$truth >= zo$conf.low & zo$truth <= zo$conf.high
  cov <- mean(covered)
  ci <- wilson(sum(covered), n_ok)
  bias <- mean(zo$estimate - zo$truth)
  emp_se <- stats::sd(zo$estimate)
  bias_se <- bias / emp_se
  n_failed <- nrow(z) - n_ok
  pass_cov <- cov >= 0.93 && cov <= 0.97
  pass_bias <- abs(bias_se) < 0.1
  pass_fit <- n_failed <= 20L
  missed <- c(
    if (!pass_cov) "coverage",
    if (!pass_bias) "bias",
    if (!pass_fit) "fits"
  )
  data.frame(
    design = key$design,
    design_name = designs[[key$design]]$name,
    distribution = "nbinom2",
    family = key$family,
    scale = key$scale,
    basis = key$basis,
    n_seeds = nrow(z),
    n_ok = n_ok,
    n_failed = n_failed,
    n_warned = sum(z$n_warn > 0L),
    truth_mean = mean(z$truth),
    coverage = cov,
    coverage_lo = ci[1L],
    coverage_hi = ci[2L],
    bias = bias,
    emp_se = emp_se,
    bias_over_se = bias_se,
    mean_reported_se = mean(zo$std.error),
    mean_width = mean(zo$conf.high - zo$conf.low),
    verdict = if (length(missed)) {
      paste0("FAIL (", paste(missed, collapse = ", "), ")")
    } else {
      "PASS"
    },
    stringsAsFactors = FALSE
  )
}))
summ$scale <- factor(summ$scale, c("difference", "ratio", "rii", "lnrr"))
summ <- summ[order(summ$design, summ$family, summ$scale, summ$basis), ]
summ$scale <- as.character(summ$scale)
rownames(summ) <- NULL

t_end <- Sys.time()
wall <- difftime(t_end, t_start, units = "mins")

# Write -----------------------------------------------------------------------

utils::write.csv(
  summ,
  file.path(out_dir, "planted_truth_S1_results.csv"),
  row.names = FALSE
)
gz <- gzfile(file.path(out_dir, "planted_truth_S1_seeds.csv.gz"), "w")
utils::write.csv(seeds, gz, row.names = FALSE)
close(gz)

warn_tab <- table(seeds$first_warn[!is.na(seeds$first_warn)])
err_tab <- table(seeds$error[!is.na(seeds$error)])
sink(file.path(out_dir, "planted_truth_S1_session.txt"))
cat("Planted-truth suite S1 (rule: validation/recovery_rule_S1.md)\n")
cat("Started :", format(t_start, "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Finished:", format(t_end, "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Wall time (min):", round(as.numeric(wall), 2), "\n")
cat("Seeds per cell:", n_seeds, "; cores:", cores, "; tasks:", length(tasks), "\n")
cat("Fits run:", nrow(seeds), "\n")
cat("pollinatoR git HEAD:", system("git rev-parse HEAD", intern = TRUE), "\n")
cat("\nPackage versions:\n")
for (p in c("pollinatoR", "glmmTMB", "TMB", "marginaleffects", "Matrix")) {
  cat(sprintf("  %-16s %s\n", p, as.character(utils::packageVersion(p))))
}
cat("\nErrors (first message, count):\n")
if (length(err_tab)) print(err_tab) else cat("  none\n")
cat("\nWarnings (first per seed-cell, count; marginaleffects advisory excluded):\n")
if (length(warn_tab)) print(warn_tab) else cat("  none\n")
cat("\n")
print(utils::sessionInfo())
sink()

print(summ[c(
  "design", "family", "scale", "basis", "n_ok", "n_failed",
  "coverage", "bias_over_se", "verdict"
)], digits = 3)
message("Wall time: ", round(as.numeric(wall), 2), " min")
