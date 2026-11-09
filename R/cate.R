# Conditional effects with partial pooling (layer L4) ---------------------------

#' Estimate the effect at each level of a grouping, with partial pooling
#'
#' Reads the estimand's contrast at each level of a column that carries a
#' random treatment slope: the unit column of the design's slope
#' (`pollin_design(slope = )`, a slope on the outermost unit), or the column
#' declared with `pollin_fit(pooled = )`. Each level's effect is formed by
#' G-computation over the rows at that level, every row keeping its random
#' effects, so a level with little data is pulled towards the overall
#' effect.
#'
#' On the glmmTMB engine the estimate uses the fixed effects and the
#' conditional modes of the random effects. The interval is approximate: the
#' fixed effects are drawn from their sampling distribution and, independently
#' of them, each coefficient of the slope term from its conditional
#' variance, the other random terms held at their modes; it is the percentile
#' interval of the effects recomputed per draw, and `std.error` their standard
#' deviation. On the flexyBayes engine the effects come from the posterior
#' draws of [pollin_draws()]: the estimate is their median and the interval
#' their equal-tailed interval.
#'
#' The table holds pooled contrasts only. A fit that pools a column does not
#' also fit its levels as fixed effects, so the unpooled per-level contrasts
#' come from a second fit, `pollin_fit(estimand, by = )` without the slope on
#' that column, read with [pollin_effects()]; the two side by side show what
#' the pooling moved. The print names which table is which.
#'
#' @param fit A [pollin_fit()] object with a random treatment slope.
#' @param by The column whose levels are reported. `NULL` (the default) takes
#'   the fit's `pooled` column, or else the unit column of the design's slope.
#' @param conf_level Level of the interval.
#' @param n Number of draws for the interval.
#' @param seed Optional seed for the draws. The session's random number state
#'   is restored afterwards.
#'
#' @return An object of class `pollin_cate`: a list holding `table`, with one
#'   row per response and level and the columns `response`, `contrast`,
#'   `level`, `estimate`, `std.error`, `conf.low`, `conf.high`, `pooled`
#'   (`TRUE`), `n` (rows of the model data at the level), `family`, `scale`,
#'   `basis`, `distribution` and `interval` (`"conditional"` on glmmTMB,
#'   `"posterior"` on flexyBayes); and `by`, `engine`, `conf_level`, `n`,
#'   `seed` and the `estimand`.
#' @seealso [pollin_effects()] for the trial contrast and the fixed per-level
#'   contrasts.
#' @export
#' @examples
#' \donttest{
#' set.seed(3)
#' d <- data.frame(
#'   trt = rep(rep(c("Cage", "Open"), each = 8L), 12L),
#'   Site = rep(sprintf("s%02d", seq_len(12L)), each = 16L)
#' )
#' slope <- rnorm(12L, 0, 0.3)[as.integer(factor(d$Site))]
#' d$pods <- rnbinom(
#'   nrow(d),
#'   mu = exp(3 + (0.3 + slope) * (d$trt == "Open")),
#'   size = 20
#' )
#' pd <- pollin_data(d, "pods", "trt", unit = "Site")
#' ds <- pollin_design(pd, c(excluded = "Cage", open = "Open"), "multi_site")
#' est <- pollin_estimand(ds, "dependence", "ratio", "loss")
#' pollin_cate(pollin_fit(est), n = 1000L, seed = 1L)
#' if (requireNamespace("flexyBayes", quietly = TRUE) &&
#'   requireNamespace("INLA", quietly = TRUE)) {
#'   fit_b <- pollin_fit(est, engine = "flexyBayes")
#'   pollin_cate(fit_b, n = 1000L, seed = 1L)
#' }
#' }
pollin_cate <- function(
  fit,
  by = NULL,
  conf_level = 0.95,
  n = 4000L,
  seed = NULL
) {
  .check_class(fit, "pollin_fit", arg = "fit")
  by <- .cate_by(fit, by)
  .check_conf_level(conf_level)
  .draws_check_n(n)
  .check_seed(seed)

  # Pooled contrasts, one block of rows per response ---------------------------
  est <- fit$estimand
  rows <- lapply(names(fit$fits), function(r) {
    tab <- if (fit$engine == "flexyBayes") {
      .cate_posterior(fit, r, by, n, seed, conf_level)
    } else {
      .cate_conditional(fit, r, by, n, seed, conf_level)
    }
    counts <- table(as.character(fit$data[[r]][[by]]))
    cbind(
      data.frame(
        response = r,
        contrast = paste(est$contrast[["high"]], "vs", est$contrast[["low"]]),
        stringsAsFactors = FALSE
      ),
      tab,
      data.frame(
        pooled = TRUE,
        n = as.integer(counts[tab$level]),
        family = est$family,
        scale = est$scale,
        basis = est$basis,
        distribution = fit$families[[r]]$family,
        interval = if (fit$engine == "flexyBayes") {
          "posterior"
        } else {
          "conditional"
        },
        stringsAsFactors = FALSE
      )
    )
  })
  tab <- do.call(rbind, rows)
  rownames(tab) <- NULL

  structure(
    list(
      table = tab,
      by = by,
      engine = fit$engine,
      conf_level = conf_level,
      n = as.integer(n),
      seed = seed,
      estimand = est
    ),
    class = "pollin_cate"
  )
}

#' @export
print.pollin_cate <- function(x, ...) {
  e <- x$estimand
  cat("<pollin_cate>\n")
  cat("  by       :", x$by, "\n")
  cat("  scale    :", e$scale, "\n")
  cat("  basis    :", e$basis, "\n")
  cat(
    "  estimate :",
    .scale_definition(
      e$scale,
      e$basis,
      e$contrast[["high"]],
      e$contrast[["low"]]
    ),
    "\n"
  )
  cat(
    "  interval :",
    sprintf("%g%%", 100 * x$conf_level),
    if (x$engine == "flexyBayes") {
      paste0("posterior, equal-tailed (", x$n, " draws)")
    } else {
      paste0(
        "approximate, fixed effects and conditional variance of the slope (",
        x$n,
        " draws)"
      )
    },
    "\n"
  )
  cat(
    "Pooled: partial pooling through the random treatment slope on",
    x$by,
    "\n"
  )
  print(x$table, row.names = FALSE)
  cat(
    "Unpooled: not carried by this fit, which pools ",
    x$by,
    "; the fixed per-level contrasts come from `pollin_fit(by = \"",
    x$by,
    "\")` without the slope, read with pollin_effects().\n",
    sep = ""
  )
  invisible(x)
}

# Internal ---------------------------------------------------------------------

#' The column to report: `by` when it carries a random treatment slope, or
#' by default the fit's pooled column, else the design's slope unit
#' @noRd
#' @keywords internal
.cate_by <- function(fit, by) {
  slope <- fit$estimand$design$slope
  if (is.null(by)) {
    by <- if (is.null(fit$pooled)) slope else fit$pooled
    if (is.null(by)) {
      stop(
        "The fit has no random treatment slope; declare one with ",
        "`pollin_design(slope = )` or `pollin_fit(pooled = )`.",
        call. = FALSE
      )
    }
  }
  if (!is.character(by) || length(by) != 1L) {
    stop("`by` must be NULL or one column name.", call. = FALSE)
  }
  if (identical(by, fit$pooled)) {
    return(by)
  }
  if (identical(by, slope)) {
    if (!identical(fit$slope_key, slope)) {
      stop(
        "The design's slope is on the nested unit ",
        fit$slope_key,
        "; pollin_cate() reads a slope on one column.",
        call. = FALSE
      )
    }
    return(by)
  }
  stop(
    "`by = \"",
    by,
    "\"` carries no random treatment slope in this fit; declare it with ",
    "`pollin_fit(pooled = \"",
    by,
    "\")`.",
    call. = FALSE
  )
}

#' Per-level effects on the glmmTMB engine: the estimate from the fixed
#' effects and conditional modes, the interval from draws of the fixed
#' effects and of the slope coefficients from their conditional variances
#' @noRd
#' @keywords internal
.cate_conditional <- function(fit, response, by, n, seed, conf_level) {
  e <- fit$estimand
  model <- fit$fits[[response]]
  d <- fit$data[[response]]
  trt <- fit$treatment
  point <- .arm_means(
    model,
    d,
    trt,
    e$contrast[["high"]],
    e$contrast[["low"]],
    by
  )
  estimate <- .effect_from_means(point, e$scale, e$basis)

  # Draws of the fixed effects and of the slope term ---------------------------
  beta <- glmmTMB::fixef(model)$cond
  v <- stats::vcov(model, full = FALSE)$cond
  r <- tryCatch(chol(v), error = function(err) NULL)
  re <- glmmTMB::ranef(model, condVar = TRUE)$cond[[by]]
  cv <- attr(re, "condVar")
  sds <- if (is.null(cv)) NULL else sqrt(apply(cv, 3L, diag))
  if (is.null(r) || is.null(sds) || any(!is.finite(sds))) {
    stop(
      "The fit of '",
      response,
      "' gives no covariance for its fixed effects or no conditional ",
      "variance for the slope on '",
      by,
      "'; the slope variance may be at zero.",
      call. = FALSE
    )
  }
  sds <- matrix(sds, nrow = ncol(re))
  p <- length(beta)
  k <- ncol(re)
  g <- nrow(re)
  z <- matrix(.draws_rnorm(n * (p + k * g), seed), nrow = n)
  b <- z[, seq_len(p), drop = FALSE] %*% r + rep(beta, each = n)
  shift <- z[, p + seq_len(k * g), drop = FALSE] *
    rep(as.vector(sds), each = n)

  # Arm means per draw ---------------------------------------------------------
  grp <- match(as.character(d[[by]]), rownames(re))
  lev <- as.character(d[[by]])
  arm <- function(level) {
    a <- .draws_arm_linear(model, d, trt, level, beta)
    nd <- d
    nd[[trt]] <- factor(rep(level, nrow(d)), levels = levels(d[[trt]]))
    zr <- stats::model.matrix(
      stats::reformulate(trt),
      nd
    )[, colnames(re), drop = FALSE]
    eta <- a$x %*% t(b) + a$re
    for (j in seq_len(k)) {
      cols <- (seq_len(g) - 1L) * k + j
      eta <- eta + zr[, j] * t(shift[, cols, drop = FALSE])[grp, , drop = FALSE]
    }
    mu <- stats::family(model)$linkinv(eta)
    t(rowsum(mu, lev) / as.vector(table(lev)))
  }
  m <- list(high = arm(e$contrast[["high"]]), low = arm(e$contrast[["low"]]))
  .cate_table(estimate, .effect_from_means(m, e$scale, e$basis), conf_level)
}

#' Per-level effects on the flexyBayes engine, from the posterior draws
#' @noRd
#' @keywords internal
.cate_posterior <- function(fit, response, by, n, seed, conf_level) {
  e <- fit$estimand
  m <- .flexy_arm_draws(fit, response, n, seed, by)
  reps <- .effect_from_means(m, e$scale, e$basis)
  .cate_table(apply(reps, 2L, stats::median), reps, conf_level)
}

#' The per-level table from point estimates and a matrix of draws, one column
#' per level
#' @noRd
#' @keywords internal
.cate_table <- function(estimate, reps, conf_level) {
  a <- (1 - conf_level) / 2
  q <- apply(reps, 2L, stats::quantile, probs = c(a, 1 - a), names = FALSE)
  data.frame(
    level = colnames(reps),
    estimate = unname(as.vector(estimate[colnames(reps)])),
    std.error = unname(apply(reps, 2L, stats::sd)),
    conf.low = q[1L, ],
    conf.high = q[2L, ],
    stringsAsFactors = FALSE
  )
}
