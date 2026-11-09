# Planted-truth recovery suite, Stage 2 (criterion C2: value intervals).
#
# Implements validation/recovery_rule_S2.md exactly. Run from the package root:
#
#   Rscript --vanilla validation/planted_truth_S2.R [n_seeds] [n_seeds_boot] [cores] [out_dir]
#
# Defaults: 1000 seeds, 200 bootstrap seeds, 8 cores, out_dir = "validation". A
# smoke run uses smaller seed counts and a scratch out_dir; only the full run
# writes the results named in the rule.

t_start <- Sys.time()
args <- commandArgs(trailingOnly = TRUE)
n_seeds <- if (length(args) >= 1L) as.integer(args[1L]) else 1000L
n_seeds_boot <- if (length(args) >= 2L) as.integer(args[2L]) else 200L
cores <- if (length(args) >= 3L) as.integer(args[3L]) else 8L
out_dir <- if (length(args) >= 4L) args[4L] else "validation"
stopifnot(n_seeds >= 1L, n_seeds_boot >= 0L, cores >= 1L, dir.exists(out_dir))

suppressMessages(devtools::load_all(".", quiet = TRUE))

n_draws <- 2000L
n_boot <- 199L
market_point <- c(price = 650, quantity = 12000)
market_sdlog <- c(price = 0.10, quantity = 0.05)

# Designs (rule Section 2) -----------------------------------------------------
# `mean_of(arm, data)` gives exp(beta[arm]) for every row of the data.

designs <- list(
  D1 = list(
    id = 1L,
    name = "three-arm trial",
    layout = "rcbd",
    unit = "block",
    arms = c(excluded = "Cage", open = "Open", supplemented = "Hand"),
    mean = c(Cage = 20, Open = 30, Hand = 36),
    theta = 8,
    channels = c("dependence", "deficit"),
    simulate = function() {
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
    channels = c("dependence", "deficit"),
    simulate = function() {
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
    channels = "dependence",
    simulate = function() {
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
  ),
  D4 = list(
    id = 4L,
    name = "two-stratum population",
    layout = "rcbd",
    unit = "block",
    arms = c(excluded = "Cage", open = "Open"),
    mean_band = rbind(
      near = c(Cage = 20, Open = 32),
      far = c(Cage = 20, Open = 24)
    ),
    theta = 8,
    channels = "dependence",
    population = pollin_population("band", c(near = 1, far = 3)),
    simulate = function() {
      blocks <- c(sprintf("n%02d", 1:10), sprintf("f%02d", 1:5))
      g <- expand.grid(
        plant = 1:4,
        trt = c("Cage", "Open"),
        block = blocks,
        stringsAsFactors = FALSE
      )
      g$band <- ifelse(substr(g$block, 1L, 1L) == "n", "near", "far")
      u_block <- stats::rnorm(15L, 0, 0.40)
      names(u_block) <- blocks
      list(data = g, u = unname(u_block[g$block]))
    }
  )
)

mean_of <- function(des, arm, d) {
  if (is.null(des$mean_band)) {
    rep(des$mean[[arm]], nrow(d))
  } else {
    unname(des$mean_band[d$band, arm])
  }
}

# Simulation and truth (rule Sections 1, 3, 5) ---------------------------------

simulate_design <- function(des, s) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(2000000L + 100000L * des$id + s)
  sim <- des$simulate()
  d <- sim$data
  mu <- mean_of(des, "Cage", d)
  for (a in unique(d$trt)) {
    mu[d$trt == a] <- mean_of(des, a, d)[d$trt == a]
  }
  mu <- mu * exp(sim$u)
  d$y <- stats::rnbinom(nrow(d), mu = mu, size = des$theta)
  mk <- list(
    price_true = stats::rlnorm(1L, log(market_point[["price"]]), market_sdlog[["price"]]),
    quantity_true = stats::rlnorm(
      1L,
      log(market_point[["quantity"]]),
      market_sdlog[["quantity"]]
    ),
    price_draws = stats::rlnorm(
      n_draws,
      log(market_point[["price"]]),
      market_sdlog[["price"]]
    ),
    quantity_draws = stats::rlnorm(
      n_draws,
      log(market_point[["quantity"]]),
      market_sdlog[["quantity"]]
    )
  )
  list(data = d, u = sim$u, market = mk)
}

contrast_of <- function(des, channel) {
  if (channel == "deficit") {
    c(high = des$arms[["supplemented"]], low = des$arms[["open"]])
  } else {
    c(high = des$arms[["open"]], low = des$arms[["excluded"]])
  }
}

truth_ratio <- function(des, sim, channel) {
  hl <- contrast_of(des, channel)
  d <- sim$data
  h <- mean_of(des, hl[["high"]], d) * exp(sim$u)
  l <- mean_of(des, hl[["low"]], d) * exp(sim$u)
  if (is.null(des$population)) {
    H <- mean(h)
    L <- mean(l)
  } else {
    w <- des$population$weights
    H <- sum(w * tapply(h, d$band, mean)[names(w)])
    L <- sum(w * tapply(l, d$band, mean)[names(w)])
  }
  if (channel == "dependence") 1 - L / H else H / L - 1
}

markets <- function(sim) {
  list(
    off = pollin_market(
      market_point[["quantity"]],
      market_point[["price"]],
      "t",
      "AUD",
      2024,
      "Planted-truth suite S2, market off"
    ),
    on = pollin_market(
      market_point[["quantity"]],
      market_point[["price"]],
      "t",
      "AUD",
      2024,
      "Planted-truth suite S2, market on",
      quantity_draws = sim$market$quantity_draws,
      price_draws = sim$market$price_draws
    )
  )
}

# One channel on one simulated data set (public API only) ----------------------

me_advisory <- "only takes into account the uncertainty in fixed-effect"

empty_row <- function() {
  data.frame(
    median = NA_real_,
    conf.low = NA_real_,
    conf.high = NA_real_,
    sd_draws = NA_real_,
    ratio_median = NA_real_,
    ratio_low = NA_real_,
    ratio_high = NA_real_,
    error = NA_character_,
    n_warn = 0L,
    first_warn = NA_character_,
    stringsAsFactors = FALSE
  )
}

# Returns one row per market in `which_market`
run_channel <- function(des, sim, channel, method, which_market, draws_seed) {
  warns <- character(0)
  res <- tryCatch(
    withCallingHandlers(
      {
        pd <- pollin_data(sim$data, "y", "trt", unit = des$unit)
        ds <- pollin_design(pd, des$arms, des$layout)
        basis <- if (channel == "dependence") "loss" else "gain"
        pop <- if (is.null(des$population)) "trial" else des$population
        est <- pollin_estimand(ds, channel, "ratio", basis, population = pop)
        fit <- pollin_fit(
          est,
          engine = "glmmTMB",
          family = glmmTMB::nbinom2(),
          by = if (is.null(des$population)) NULL else des$population$by
        )
        n <- if (method == "simulation") n_draws else n_boot
        dr <- pollin_draws(fit, n = n, method = method, seed = draws_seed)
        level <- if (is.null(des$population)) "overall" else "population"
        rs <- summary(dr)
        rs <- rs[rs$level == level, ]
        mk <- markets(sim)
        lapply(which_market, function(m) {
          market <- mk[[m]]
          if (m == "on" && dr$n != n_draws) {
            stop("Market draws need ", n_draws, " estimand draws.")
          }
          v <- pollin_value(dr, market, channel = channel, level = level)
          s <- summary(v)
          list(
            median = s$median,
            conf.low = s$conf.low,
            conf.high = s$conf.high,
            sd_draws = stats::sd(v$draws),
            ratio_median = rs$median,
            ratio_low = rs$conf.low,
            ratio_high = rs$conf.high
          )
        })
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
  out <- do.call(rbind, lapply(seq_along(which_market), function(j) {
    r <- empty_row()
    r$n_warn <- length(warns)
    r$first_warn <- if (length(warns)) warns[1L] else NA_character_
    if (is.character(res)) {
      r$error <- res
    } else {
      r[names(res[[j]])] <- res[[j]]
    }
    r
  }))
  cbind(data.frame(market = which_market, stringsAsFactors = FALSE), out)
}

# Tasks (rule Sections 4, 5) ---------------------------------------------------

run_task <- function(task) {
  des <- designs[[task$design]]
  s <- task$seed
  sim <- simulate_design(des, s)
  rows <- list()
  if (task$method == "simulation") {
    which_market <- if (is.null(des$population)) c("off", "on") else "off"
    for (ch in des$channels) {
      r <- run_channel(
        des,
        sim,
        ch,
        "simulation",
        which_market,
        draws_seed = 1000000L * des$id + s
      )
      rows[[length(rows) + 1L]] <- cbind(
        data.frame(design = task$design, channel = ch, method = "simulation"),
        r,
        data.frame(ratio_truth = truth_ratio(des, sim, ch))
      )
    }
  } else {
    r <- run_channel(
      des,
      sim,
      "dependence",
      "bootstrap",
      "off",
      draws_seed = 5000000L + s
    )
    rows[[1L]] <- cbind(
      data.frame(design = task$design, channel = "dependence", method = "bootstrap"),
      r,
      data.frame(ratio_truth = truth_ratio(des, sim, "dependence"))
    )
  }
  out <- do.call(rbind, rows)
  pq_true <- ifelse(
    out$market == "on",
    sim$market$price_true * sim$market$quantity_true,
    market_point[["price"]] * market_point[["quantity"]]
  )
  out$truth <- out$ratio_truth * pq_true
  out$seed <- s
  out
}

tasks <- c(
  unlist(
    lapply(names(designs), function(k) {
      lapply(seq_len(n_seeds), function(s) {
        list(design = k, seed = s, method = "simulation")
      })
    }),
    recursive = FALSE
  ),
  lapply(seq_len(n_seeds_boot), function(s) {
    list(design = "D1", seed = s, method = "bootstrap")
  })
)
# Bootstrap tasks are the slowest; start them first
tasks <- tasks[order(vapply(tasks, function(t) t$method != "bootstrap", logical(1)))]

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
seeds <- do.call(rbind, res)

# Cell labels (rule Section 4) -------------------------------------------------

cell_key <- paste(seeds$design, seeds$channel, seeds$market, seeds$method)
cell_ids <- c(
  "D1 dependence off simulation" = "V01",
  "D1 dependence on simulation" = "V02",
  "D1 deficit off simulation" = "V03",
  "D1 deficit on simulation" = "V04",
  "D2 dependence off simulation" = "V05",
  "D2 dependence on simulation" = "V06",
  "D2 deficit off simulation" = "V07",
  "D2 deficit on simulation" = "V08",
  "D3 dependence off simulation" = "V09",
  "D3 dependence on simulation" = "V10",
  "D4 dependence off simulation" = "V11",
  "D1 dependence off bootstrap" = "B01"
)
seeds$cell <- unname(cell_ids[cell_key])
stopifnot(!anyNA(seeds$cell))
seeds <- seeds[order(seeds$cell, seeds$seed), c(
  "cell", "design", "channel", "market", "method", "seed", "truth",
  "ratio_truth", "median", "conf.low", "conf.high", "sd_draws",
  "ratio_median", "ratio_low", "ratio_high", "error", "n_warn", "first_warn"
)]
rownames(seeds) <- NULL

# Summary per cell (rule Section 7) --------------------------------------------

wilson <- function(x, n, z = stats::qnorm(0.975)) {
  p <- x / n
  centre <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  c(centre - half, centre + half)
}

summ <- do.call(rbind, lapply(sort(unique(seeds$cell)), function(id) {
  z <- seeds[seeds$cell == id, ]
  ok <- is.na(z$error) &
    is.finite(z$median) & is.finite(z$conf.low) &
    is.finite(z$conf.high) & is.finite(z$sd_draws)
  zo <- z[ok, ]
  n_ok <- nrow(zo)
  covered <- zo$truth >= zo$conf.low & zo$truth <= zo$conf.high
  cov <- mean(covered)
  ci <- wilson(sum(covered), n_ok)
  bias <- mean(zo$median - zo$truth)
  se <- mean(zo$sd_draws)
  bias_se <- bias / se
  n_failed <- nrow(z) - n_ok
  fail_limit <- if (id == "B01") 4L else 20L
  pass_cov <- cov >= 0.93 && cov <= 0.97
  pass_bias <- abs(bias_se) < 0.1
  pass_fit <- n_failed <= fail_limit
  missed <- c(
    if (!pass_cov) "coverage",
    if (!pass_bias) "bias",
    if (!pass_fit) "fits"
  )
  ratio_cov <- mean(zo$ratio_truth >= zo$ratio_low & zo$ratio_truth <= zo$ratio_high)
  data.frame(
    cell = id,
    design = z$design[1L],
    design_name = designs[[z$design[1L]]]$name,
    channel = z$channel[1L],
    market = z$market[1L],
    method = z$method[1L],
    n_seeds = nrow(z),
    n_ok = n_ok,
    n_failed = n_failed,
    n_warned = sum(z$n_warn > 0L),
    truth_mean = mean(z$truth),
    ratio_truth_mean = mean(z$ratio_truth),
    coverage = cov,
    coverage_mcse = sqrt(cov * (1 - cov) / n_ok),
    coverage_lo = ci[1L],
    coverage_hi = ci[2L],
    bias = bias,
    se = se,
    bias_over_se = bias_se,
    mc_sd_median = stats::sd(zo$median),
    median_width = stats::median(zo$conf.high - zo$conf.low),
    ratio_coverage = ratio_cov,
    ratio_bias = mean(zo$ratio_median - zo$ratio_truth),
    in_c2 = id != "B01",
    verdict = if (length(missed)) {
      paste0("FAIL (", paste(missed, collapse = ", "), ")")
    } else {
      "PASS"
    },
    stringsAsFactors = FALSE
  )
}))
rownames(summ) <- NULL
c2_met <- all(summ$verdict[summ$in_c2] == "PASS")

t_end <- Sys.time()
wall <- difftime(t_end, t_start, units = "mins")

# Write -----------------------------------------------------------------------

utils::write.csv(
  summ,
  file.path(out_dir, "planted_truth_S2_results.csv"),
  row.names = FALSE
)
gz <- gzfile(file.path(out_dir, "planted_truth_S2_seeds.csv.gz"), "w")
utils::write.csv(seeds, gz, row.names = FALSE)
close(gz)

warn_tab <- table(seeds$first_warn[!is.na(seeds$first_warn)])
err_tab <- table(seeds$error[!is.na(seeds$error)])
sink(file.path(out_dir, "planted_truth_S2_session.txt"))
cat("Planted-truth suite S2 (rule: validation/recovery_rule_S2.md)\n")
cat("Started :", format(t_start, "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Finished:", format(t_end, "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Wall time (min):", round(as.numeric(wall), 2), "\n")
cat(
  "Seeds per cell:", n_seeds, "; bootstrap seeds:", n_seeds_boot,
  "; cores:", cores, "; tasks:", length(tasks), "\n"
)
cat("Estimand draws:", n_draws, "; bootstrap refits:", n_boot, "\n")
cat("C2 met (V01-V11 all PASS):", c2_met, "\n")
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
  "cell", "design", "channel", "market", "method", "n_ok", "n_failed",
  "coverage", "coverage_mcse", "bias_over_se", "median_width",
  "ratio_coverage", "verdict"
)], digits = 3)
message("C2 met: ", c2_met)
message("Wall time: ", round(as.numeric(wall), 2), " min")
