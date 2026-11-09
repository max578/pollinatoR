# Effects on the estimand's scale (layer L4) -----------------------------------

#' Estimate the effect on the estimand's scale
#'
#' Reads the estimand's contrast from each fitted model by G-computation: every
#' row of the model data is predicted under the `high` arm and under the `low`
#' arm, the predictions are averaged, and the effect is formed from the two
#' averaged predictions. A ratio-type scale is a ratio of averaged predictions,
#' never an average of unit ratios. With `rr = high / low`:
#'
#' | Scale | Basis | Estimate |
#' |---|---|---|
#' | `"difference"` | either | `high - low` |
#' | `"ratio"` | `"loss"` | `D = 1 - 1 / rr` |
#' | `"ratio"` | `"gain"` | `rr - 1` |
#' | `"rii"` | either | `(rr - 1) / (rr + 1)`, equal to `D / (2 - D)` |
#' | `"lnrr"` | either | `log(rr)` |
#' | `"ratio"` | none (cage family) | `rr`, the cage factor `sham / open` |
#'
#' With `interval = "delta"`, the interval for a ratio-type scale is a delta
#' interval for `log(rr)`, transformed to the scale, so the loss share stays
#' below 1 and the RII inside (-1, 1); its p-value tests `rr = 1`. The delta
#' interval is taken over the fixed effects and the random effects jointly,
#' with their covariance from the inverse of the joint precision of the
#' Laplace approximation: the estimand averages over the blocks and slopes of
#' this trial, so the random effects enter at their realised values, with
#' their uncertainty and their correlation with the fixed effects. A
#' zero-inflated fit, or one whose joint precision is singular (a variance at
#' zero), is an error that points to the bootstrap. `interval = "delta_fixed"`
#' is the delta interval of pollinatoR 0.3.1 and earlier, from
#' `marginaleffects`: it varies the fixed effects and holds the random effects
#' at their modes, so it is wider than the joint interval wherever the
#' estimand depends on both, as on a design with a random treatment slope or
#' under a logit link. With `interval = "bootstrap"`, responses are simulated
#' from the fitted model (new random effects included), the model is refitted
#' to each simulated set,
#' and the interval is the percentile interval of the refitted estimates;
#' `std.error` is their standard deviation and `p.value` is `NA`.
#'
#' A fit with `engine = "flexyBayes"` has the `"posterior"` interval, its
#' default and only one: the estimate is the median of the posterior draws of
#' the estimand from [pollin_draws()], the interval their equal-tailed
#' interval, `std.error` their standard deviation, and `p.value` is `NA`.
#' With the same `seed` and `n_draws`, the row equals the summary of
#' `pollin_draws(fit, n = n_draws, seed = seed)`.
#'
#' When the data declare a denominator ([pollin_data()]), each response is a
#' count of successes out of that many trials, fitted as binomial, and the
#' arm means are expected proportions (successes / trials), averaged over the
#' rows with equal weight whatever their number of trials; the contrast is
#' between those averaged proportions.
#'
#' The table reports the trial contrast, overall or per level of `by`. A
#' population declared with [pollin_population()] is not used here; the
#' population estimand is read with [pollin_draws()].
#'
#' @param fit A [pollin_fit()] object.
#' @param by The column the fit was made with, `pollin_fit(..., by = )`.
#'   Defaults to the fit's `by`, giving one contrast per level from the fixed
#'   treatment-by-level effects, each averaged over the rows at that level.
#'   `NULL` gives one contrast averaged over all rows. A column the fit was not
#'   made with is an error: refit with `pollin_fit(..., by = )`.
#' @param interval `"delta"` (the default for a glmmTMB fit),
#'   `"delta_fixed"`, `"bootstrap"`, or `"posterior"` (for a flexyBayes fit).
#' @param n_boot Number of parametric bootstrap refits, used when
#'   `interval = "bootstrap"`.
#' @param seed Optional integer seed for the bootstrap simulation or the
#'   posterior draws. Each response is drawn from this seed, and the
#'   session's random number state is restored afterwards.
#' @param conf_level Confidence level of the interval.
#' @param n_draws Number of posterior draws, used when
#'   `interval = "posterior"`.
#'
#' @return An object of class `pollin_effects`: a list holding `table`, with
#'   one row per response (and per `by` level) and the columns `response`,
#'   `contrast`, any `by` column, `estimate`, `std.error`, `conf.low`,
#'   `conf.high`, `p.value`, `family` (the estimand family), `scale`, `basis`,
#'   `effect` (`"total"` or `"direct"`), `distribution` (the model family) and
#'   `interval`; and `interval`, `n_boot`,
#'   `seed`, `conf_level`, `by`, the `estimand` and `n_draws` (posterior
#'   intervals only).
#' @seealso [pollin_estimand()] for the scales and the basis.
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- data.frame(
#'   trt = rep(rep(c("Cage", "Open"), each = 10L), 8L),
#'   Site = rep(paste0("s", seq_len(8L)), each = 20L)
#' )
#' d$pods <- rnbinom(nrow(d), mu = exp(3 + 0.3 * (d$trt == "Open")), size = 10)
#' pd <- pollin_data(d, "pods", "trt", unit = "Site")
#' ds <- pollin_design(pd, c(excluded = "Cage", open = "Open"), "rcbd")
#' est <- pollin_estimand(ds, "dependence", "ratio", "loss")
#' fit <- pollin_fit(est)
#' pollin_effects(fit)
#' pollin_effects(fit, interval = "bootstrap", n_boot = 50L, seed = 1L)
#' }
pollin_effects <- function(
  fit,
  by = fit$by,
  interval = c("delta", "delta_fixed", "bootstrap", "posterior"),
  n_boot = 999L,
  seed = NULL,
  conf_level = 0.95,
  n_draws = 4000L
) {
  .check_class(fit, "pollin_fit", arg = "fit")
  if (missing(interval)) {
    interval <- if (fit$engine == "flexyBayes") "posterior" else "delta"
  }
  .check_choice(
    interval,
    "interval",
    c("delta", "delta_fixed", "bootstrap", "posterior")
  )
  .effects_check_interval(interval, fit$engine)
  .check_conf_level(conf_level)
  if (interval == "bootstrap") {
    .check_n_boot(n_boot)
    .check_seed(seed)
  }
  if (interval == "posterior") {
    .draws_check_n(n_draws)
    .check_seed(seed)
  }
  .check_effects_by(by, fit$by)
  for (d in fit$data) {
    .check_by(by, d)
  }

  # One block of rows per response ---------------------------------------------
  est <- fit$estimand
  rows <- lapply(names(fit$fits), function(r) {
    tab <- if (interval == "delta") {
      .effects_joint(fit, r, est, by, conf_level)
    } else if (interval == "delta_fixed") {
      .muffle_fixed_only(
        .effects_delta(fit$fits[[r]], fit$treatment, est, by, conf_level)
      )
    } else if (interval == "posterior") {
      .effects_posterior(fit, r, by, n_draws, seed, conf_level)
    } else {
      .effects_boot(fit, r, by, n_boot, seed, conf_level)
    }
    cbind(
      data.frame(
        response = rep(r, nrow(tab)),
        contrast = rep(
          paste(est$contrast[["high"]], "vs", est$contrast[["low"]]),
          nrow(tab)
        ),
        stringsAsFactors = FALSE
      ),
      tab,
      data.frame(
        family = est$family,
        scale = est$scale,
        basis = est$basis,
        effect = est$effect,
        distribution = fit$families[[r]]$family,
        interval = interval,
        stringsAsFactors = FALSE
      )
    )
  })
  tab <- do.call(rbind, rows)
  rownames(tab) <- NULL

  structure(
    list(
      table = tab,
      interval = interval,
      n_boot = if (interval == "bootstrap") as.integer(n_boot),
      seed = seed,
      conf_level = conf_level,
      by = by,
      estimand = est,
      n_draws = if (interval == "posterior") as.integer(n_draws)
    ),
    class = "pollin_effects"
  )
}

#' @export
print.pollin_effects <- function(x, ...) {
  e <- x$estimand
  high <- e$contrast[["high"]]
  low <- e$contrast[["low"]]
  cat("<pollin_effects>\n")
  cat("  family   :", .family_label(e$family), "\n")
  cat("  scale    :", e$scale, "\n")
  cat("  basis    :", .basis_label(e$basis), "\n")
  cat("  effect   :", .effect_label(e$effect, e$mediators), "\n")
  cat("  estimate :", .scale_definition(e$scale, e$basis, high, low), "\n")
  cat(
    "  interval :",
    sprintf("%g%%", 100 * x$conf_level),
    switch(
      x$interval,
      delta = "delta, fixed and random effects jointly",
      delta_fixed = "delta, fixed effects only",
      bootstrap = paste0("parametric bootstrap (", x$n_boot, " refits)"),
      posterior = paste0("posterior, equal-tailed (", x$n_draws, " draws)")
    ),
    "\n"
  )
  if (x$interval == "delta_fixed") {
    cat("  note     : delta intervals cover fixed-effect uncertainty only\n")
  }
  print(x$table, row.names = FALSE)
  invisible(x)
}

# Scales ----------------------------------------------------------------------

#' The estimate's definition in words, for printing
#' @noRd
#' @keywords internal
.scale_definition <- function(scale, basis, high, low) {
  r <- sprintf("%s/%s", high, low)
  switch(
    scale,
    difference = sprintf("%s - %s (difference)", high, low),
    ratio = if (is.na(basis)) {
      sprintf("%s (cage factor)", r)
    } else if (basis == "loss") {
      sprintf("1 - %s/%s (loss share)", low, high)
    } else {
      sprintf("%s - 1 (gain)", r)
    },
    rii = sprintf("(%s - %s)/(%s + %s) (RII)", high, low, high, low),
    lnrr = sprintf("log(%s) (lnRR)", r)
  )
}

#' Map the ratio `high / low` to a ratio-type scale
#' @noRd
#' @keywords internal
.from_ratio <- function(r, scale, basis) {
  switch(
    scale,
    ratio = if (is.na(basis)) {
      r
    } else if (basis == "loss") {
      1 - 1 / r
    } else {
      r - 1
    },
    rii = (r - 1) / (r + 1),
    lnrr = log(r)
  )
}

#' Derivative of `.from_ratio()` with respect to `log(r)`, for delta standard
#' errors on the scale
#' @noRd
#' @keywords internal
.from_ratio_slope <- function(r, scale, basis) {
  switch(
    scale,
    ratio = if (is.na(basis)) {
      r
    } else if (basis == "loss") {
      1 / r
    } else {
      r
    },
    rii = 2 * r / (r + 1)^2,
    lnrr = rep(1, length(r))
  )
}

# Delta intervals -------------------------------------------------------------

#' Delta-method effects for one fitted model, via `marginaleffects`
#' @noRd
#' @keywords internal
.effects_delta <- function(model, treatment, estimand, by, conf_level) {
  vars <- list(c(estimand$contrast[["low"]], estimand$contrast[["high"]]))
  names(vars) <- treatment
  ratio <- estimand$scale != "difference"
  cmp <- marginaleffects::avg_comparisons(
    model,
    variables = vars,
    comparison = if (ratio) "lnratioavg" else "difference",
    by = if (is.null(by)) TRUE else by,
    conf_level = conf_level
  )
  cmp <- as.data.frame(cmp)
  out <- data.frame(
    estimate = cmp$estimate,
    std.error = cmp$std.error,
    conf.low = cmp$conf.low,
    conf.high = cmp$conf.high,
    p.value = cmp$p.value
  )
  if (ratio) {
    s <- estimand$scale
    b <- estimand$basis
    r <- exp(cmp$estimate)
    out$estimate <- .from_ratio(r, s, b)
    out$std.error <- cmp$std.error * .from_ratio_slope(r, s, b)
    out$conf.low <- .from_ratio(exp(cmp$conf.low), s, b)
    out$conf.high <- .from_ratio(exp(cmp$conf.high), s, b)
  }
  if (!is.null(by)) {
    lv <- data.frame(as.character(cmp[[by]]), stringsAsFactors = FALSE)
    names(lv) <- by
    out <- cbind(lv, out)
  }
  out
}

#' Evaluate `expr`, muffling only the `marginaleffects` warning that its
#' intervals for mixed models cover the fixed effects only. That statement is
#' printed with every delta result instead.
#' @noRd
#' @keywords internal
.muffle_fixed_only <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      msg <- conditionMessage(w)
      if (
        grepl("marginaleffects", msg, fixed = TRUE) &&
          grepl("uncertainty in fixed-effect parameters", msg, fixed = TRUE)
      ) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

# Parametric bootstrap --------------------------------------------------------

#' Averaged predictions under the `high` and `low` arms, overall or per level
#' of `by`
#' @noRd
#' @keywords internal
.arm_means <- function(model, d, treatment, high, low, by) {
  lev <- levels(d[[treatment]])
  nd <- d
  nd[[treatment]] <- factor(rep(high, nrow(d)), levels = lev)
  p_high <- stats::predict(model, newdata = nd, type = "response")
  nd[[treatment]] <- factor(rep(low, nrow(d)), levels = lev)
  p_low <- stats::predict(model, newdata = nd, type = "response")
  g <- if (is.null(by)) rep("all", nrow(d)) else as.character(d[[by]])
  list(
    high = tapply(p_high, g, mean),
    low = tapply(p_low, g, mean)
  )
}

#' The effect on the estimand's scale from averaged arm predictions
#' @noRd
#' @keywords internal
.effect_from_means <- function(m, scale, basis) {
  if (scale == "difference") {
    return(m$high - m$low)
  }
  .from_ratio(m$high / m$low, scale, basis)
}

#' Parametric-bootstrap arm means for one response of a `pollin_fit`: the
#' point arm means, and the arm means of each refit as matrices with one row
#' per refit and one column per `by` level (failed refits are rows of `NA`).
#' Warns when refits failed or warned, and stops when fewer than two
#' succeeded.
#' @noRd
#' @keywords internal
.boot_arm_means <- function(fit, response, by, n_boot, seed) {
  e <- fit$estimand
  high <- e$contrast[["high"]]
  low <- e$contrast[["low"]]
  model <- fit$fits[[response]]
  d <- fit$data[[response]]
  trt <- fit$treatment
  point <- .arm_means(model, d, trt, high, low, by)
  k <- length(point$high)

  # Simulate, refit and recompute --------------------------------------------
  sims <- stats::simulate(model, nsim = n_boot, seed = seed)
  lhs <- .response_lhs(response, fit$denominator)
  n_warned <- 0L
  first_warning <- NULL
  reps <- vapply(
    seq_len(n_boot),
    function(b) {
      d_b <- d
      y_b <- sims[[b]]
      # A binomial fit with trials simulates successes and failures as columns.
      d_b[[response]] <- if (is.matrix(y_b)) y_b[, 1L] else y_b
      warned <- FALSE
      m_b <- tryCatch(
        withCallingHandlers(
          do.call(
            .fit_glmm,
            c(
              list(lhs, fit$rhs, d_b, fit$families[[response]]),
              fit$args
            )
          ),
          warning = function(w) {
            warned <<- TRUE
            if (is.null(first_warning)) {
              first_warning <<- conditionMessage(w)
            }
            invokeRestart("muffleWarning")
          }
        ),
        error = function(err) NULL
      )
      if (warned) {
        n_warned <<- n_warned + 1L
      }
      if (is.null(m_b)) {
        return(rep(NA_real_, 2L * k))
      }
      m <- .arm_means(m_b, d_b, trt, high, low, by)
      c(m$high, m$low)
    },
    numeric(2L * k)
  )
  reps <- matrix(reps, nrow = 2L * k)

  # Report refits that failed or warned ----------------------------------------
  n_failed <- sum(is.na(reps[1L, ]))
  if (n_failed > 0L || n_warned > 0L) {
    warning(
      "Bootstrap for '",
      response,
      "': ",
      .n_noun(n_failed, "refit"),
      " failed and ",
      n_warned,
      " gave a warning, of ",
      n_boot,
      if (!is.null(first_warning)) paste0("; first warning: ", first_warning),
      ".",
      call. = FALSE
    )
  }
  if (n_boot - n_failed < 2L) {
    stop(
      "Bootstrap for '",
      response,
      "': fewer than two refits succeeded.",
      call. = FALSE
    )
  }

  draws <- function(rows) {
    out <- t(reps[rows, , drop = FALSE])
    colnames(out) <- names(point$high)
    out
  }
  list(
    point = point,
    high = draws(seq_len(k)),
    low = draws(k + seq_len(k))
  )
}

#' Parametric-bootstrap effects for one response of a `pollin_fit`
#' @noRd
#' @keywords internal
.effects_boot <- function(fit, response, by, n_boot, seed, conf_level) {
  e <- fit$estimand
  bt <- .boot_arm_means(fit, response, by, n_boot, seed)
  point <- .effect_from_means(bt$point, e$scale, e$basis)
  reps <- t(.effect_from_means(bt, e$scale, e$basis))

  a <- (1 - conf_level) / 2
  ci <- apply(reps, 1L, stats::quantile, probs = c(a, 1 - a), na.rm = TRUE)
  out <- data.frame(
    estimate = unname(point),
    std.error = apply(reps, 1L, stats::sd, na.rm = TRUE),
    conf.low = unname(ci[1L, ]),
    conf.high = unname(ci[2L, ]),
    p.value = NA_real_
  )
  if (!is.null(by)) {
    lv <- data.frame(names(point), stringsAsFactors = FALSE)
    names(lv) <- by
    out <- cbind(lv, out)
  }
  out
}

# Posterior intervals ---------------------------------------------------------

#' Posterior effects for one response of a flexyBayes fit: the median,
#' standard deviation and equal-tailed interval of the estimand's draws
#' @noRd
#' @keywords internal
.effects_posterior <- function(fit, response, by, n_draws, seed, conf_level) {
  e <- fit$estimand
  m <- .flexy_arm_draws(fit, response, n_draws, seed, by)
  reps <- .effect_from_means(m, e$scale, e$basis)
  a <- (1 - conf_level) / 2
  q <- apply(
    reps,
    2L,
    stats::quantile,
    probs = c(0.5, a, 1 - a),
    names = FALSE
  )
  out <- data.frame(
    estimate = q[1L, ],
    std.error = apply(reps, 2L, stats::sd),
    conf.low = q[2L, ],
    conf.high = q[3L, ],
    p.value = NA_real_
  )
  if (!is.null(by)) {
    lv <- data.frame(colnames(reps), stringsAsFactors = FALSE)
    names(lv) <- by
    out <- cbind(lv, out)
  }
  out
}

# Internal checks -------------------------------------------------------------

#' Stop unless the interval suits the engine: posterior for a flexyBayes fit,
#' delta or bootstrap for a glmmTMB fit
#' @noRd
#' @keywords internal
.effects_check_interval <- function(interval, engine) {
  if (engine == "flexyBayes" && interval != "posterior") {
    stop(
      "A flexyBayes fit has posterior intervals; `interval = \"",
      interval,
      "\"` applies to glmmTMB fits.",
      call. = FALSE
    )
  }
  if (engine != "flexyBayes" && interval == "posterior") {
    stop(
      "`interval = \"posterior\"` needs a fit with `engine = \"flexyBayes\"`.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Stop unless `by` is NULL or one column of the model data
#' @noRd
#' @keywords internal
.check_by <- function(by, data) {
  if (is.null(by)) {
    return(invisible(TRUE))
  }
  if (!is.character(by) || length(by) != 1L || !by %in% names(data)) {
    stop("`by` must name one column of the data.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Stop when `by` asks for levels the fit was not made with
#' @noRd
#' @keywords internal
.check_effects_by <- function(by, fit_by) {
  if (is.null(by) || identical(by, fit_by)) {
    return(invisible(TRUE))
  }
  stop(
    "The fit was not made with `by = \"",
    paste(by, collapse = ", "),
    "\"`; refit with `pollin_fit(..., by = \"",
    paste(by, collapse = ", "),
    "\")`.",
    call. = FALSE
  )
}

#' Stop unless `conf_level` is one number strictly between 0 and 1
#' @noRd
#' @keywords internal
.check_conf_level <- function(conf_level) {
  if (
    !is.numeric(conf_level) ||
      length(conf_level) != 1L ||
      !is.finite(conf_level) ||
      conf_level <= 0 ||
      conf_level >= 1
  ) {
    stop("`conf_level` must be a number between 0 and 1.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Stop unless `n_boot` is one whole number of at least 2
#' @noRd
#' @keywords internal
.check_n_boot <- function(n_boot) {
  if (
    !is.numeric(n_boot) ||
      length(n_boot) != 1L ||
      !is.finite(n_boot) ||
      n_boot < 2 ||
      n_boot != round(n_boot)
  ) {
    stop("`n_boot` must be a whole number of at least 2.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Stop unless `seed` is NULL or one finite number
#' @noRd
#' @keywords internal
.check_seed <- function(seed) {
  if (
    !is.null(seed) &&
      (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed))
  ) {
    stop("`seed` must be NULL or one number.", call. = FALSE)
  }
  invisible(TRUE)
}
