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
#'
#' With `interval = "delta"`, the interval for a ratio-type scale is a delta
#' interval for `log(rr)`, transformed to the scale, so the loss share stays
#' below 1 and the RII inside (-1, 1); its p-value tests `rr = 1`. The delta
#' interval carries the uncertainty in the fixed effects only. With
#' `interval = "bootstrap"`, responses are simulated from the fitted model
#' (new random effects included), the model is refitted to each simulated set,
#' and the interval is the percentile interval of the refitted estimates;
#' `std.error` is their standard deviation and `p.value` is `NA`.
#'
#' @param fit A [pollin_fit()] object.
#' @param by Optional name of a column of the model data. Gives one unpooled
#'   contrast per level, each averaged over the rows at that level.
#' @param interval `"delta"` (the default) or `"bootstrap"`.
#' @param n_boot Number of parametric bootstrap refits, used when
#'   `interval = "bootstrap"`.
#' @param seed Optional integer seed for the bootstrap simulation. Each
#'   response is simulated from this seed, and the session's random number
#'   state is restored afterwards.
#' @param conf_level Confidence level of the interval.
#'
#' @return An object of class `pollin_effects`: a list holding `table`, with
#'   one row per response (and per `by` level) and the columns `response`,
#'   `contrast`, any `by` column, `estimate`, `std.error`, `conf.low`,
#'   `conf.high`, `p.value`, `family` (the estimand family), `scale`, `basis`,
#'   `distribution` (the model family) and `interval`; and `interval`, `n_boot`,
#'   `seed`, `conf_level`, `by` and the `estimand`.
#' @seealso [pollin_estimand()] for the scales and the basis.
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- data.frame(
#'   trt = rep(rep(c("Cage", "Open"), each = 10L), 8L),
#'   Site = rep(paste0("s", seq_len(8L)), each = 20L)
#' )
#' d$pods <- rpois(nrow(d), exp(3 + 0.3 * (d$trt == "Open")))
#' pd <- pollin_data(d, "pods", "trt", unit = "Site")
#' ds <- pollin_design(pd, c(excluded = "Cage", open = "Open"), "rcbd")
#' est <- pollin_estimand(ds, "dependence", "ratio", "loss")
#' fit <- pollin_fit(est)
#' pollin_effects(fit)
#' pollin_effects(fit, interval = "bootstrap", n_boot = 50L, seed = 1L)
#' }
pollin_effects <- function(
  fit,
  by = NULL,
  interval = c("delta", "bootstrap"),
  n_boot = 999L,
  seed = NULL,
  conf_level = 0.95
) {
  .check_class(fit, "pollin_fit", arg = "fit")
  if (missing(interval)) {
    interval <- "delta"
  }
  .check_choice(interval, "interval", c("delta", "bootstrap"))
  .check_conf_level(conf_level)
  if (interval == "bootstrap") {
    .check_n_boot(n_boot)
    .check_seed(seed)
  }
  for (d in fit$data) {
    .check_by(by, d)
  }

  # One block of rows per response ---------------------------------------------
  est <- fit$estimand
  rows <- lapply(names(fit$fits), function(r) {
    tab <- if (interval == "delta") {
      .effects_delta(fit$fits[[r]], fit$treatment, est, by, conf_level)
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
      estimand = est
    ),
    class = "pollin_effects"
  )
}

#' @export
print.pollin_effects <- function(x, ...) {
  e <- x$estimand
  high <- e$contrast[["high"]]
  low <- e$contrast[["low"]]
  label <- if (e$family == "deficit") "deficit (E1)" else "dependence (E2)"
  cat("<pollin_effects>\n")
  cat("  family   :", label, "\n")
  cat("  scale    :", e$scale, "\n")
  cat("  basis    :", e$basis, "\n")
  cat("  estimate :", .scale_definition(e$scale, e$basis, high, low), "\n")
  cat(
    "  interval :",
    sprintf("%g%%", 100 * x$conf_level),
    if (x$interval == "delta") {
      "delta"
    } else {
      paste0("parametric bootstrap (", x$n_boot, " refits)")
    },
    "\n"
  )
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
    ratio = if (basis == "loss") {
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
    ratio = if (basis == "loss") 1 - 1 / r else r - 1,
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
    ratio = if (basis == "loss") 1 / r else r,
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
    out <- cbind(cmp[by], out)
  }
  out
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

#' Parametric-bootstrap effects for one response of a `pollin_fit`
#' @noRd
#' @keywords internal
.effects_boot <- function(fit, response, by, n_boot, seed, conf_level) {
  e <- fit$estimand
  high <- e$contrast[["high"]]
  low <- e$contrast[["low"]]
  model <- fit$fits[[response]]
  d <- fit$data[[response]]
  trt <- fit$treatment
  point <- .effect_from_means(
    .arm_means(model, d, trt, high, low, by),
    e$scale,
    e$basis
  )

  # Simulate, refit and recompute --------------------------------------------
  sims <- stats::simulate(model, nsim = n_boot, seed = seed)
  n_warned <- 0L
  first_warning <- NULL
  reps <- vapply(
    seq_len(n_boot),
    function(b) {
      d_b <- d
      d_b[[response]] <- sims[[b]]
      warned <- FALSE
      m_b <- tryCatch(
        withCallingHandlers(
          do.call(
            .fit_glmm,
            c(
              list(response, fit$rhs, d_b, fit$families[[response]]),
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
        return(rep(NA_real_, length(point)))
      }
      .effect_from_means(
        .arm_means(m_b, d_b, trt, high, low, by),
        e$scale,
        e$basis
      )
    },
    numeric(length(point))
  )
  reps <- matrix(reps, nrow = length(point))

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

# Internal checks -------------------------------------------------------------

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
