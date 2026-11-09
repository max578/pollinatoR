# Draws of the estimand (layer L4b) -------------------------------------------

#' Draw the estimand from its sampling distribution
#'
#' Returns draws of the two G-computation arm means of a fitted estimand (the
#' predictions under the `high` and the `low` arm, averaged over the rows of
#' the model data, as in [pollin_effects()]) and the estimand formed from each
#' pair of draws on the estimand's scale and basis. A quantity computed from
#' the estimand is then computed per draw and summarised last.
#'
#' Three methods are available:
#'
#' * `"simulation"` (the default for a glmmTMB fit) draws the fixed effects from their
#'   estimated sampling distribution, a multivariate normal centred on the
#'   estimates with their covariance matrix, and recomputes the arm means for
#'   each draw. Every row keeps its estimated random effects, the prediction
#'   rule of [pollin_effects()], so the draws carry the uncertainty in the
#'   fixed effects only, like the delta interval.
#' * `"bootstrap"` keeps the arm means of the parametric bootstrap refits
#'   that [pollin_effects()] runs with `interval = "bootstrap"`: responses are
#'   simulated from the fitted model (new random effects included) and the
#'   model is refitted to each simulated set. With the same `seed` and
#'   `n = n_boot`, the percentile interval of the draws is the bootstrap
#'   interval of [pollin_effects()]. Refits that fail are dropped, with a
#'   warning, so the number of draws can be smaller than `n`.
#' * `"posterior"`, the default and the only method for a fit with
#'   `engine = "flexyBayes"`, draws the fixed effects, the random effects and
#'   the variance parameters jointly from the posterior
#'   (`flexyBayes::fb_inla_posterior_sample()`) and recomputes the arm means
#'   for each draw with the prediction rule of [pollin_effects()]: every row
#'   keeps its own random effects, here the drawn ones, and under each arm a
#'   row takes its unit's effect for that arm.
#'
#' When the estimand refers to a population declared with
#' [pollin_population()], a further column `"population"` holds the
#' population arm means, the area-weighted sums of the stratum arm means, and
#' the estimand formed from them: a ratio of population means, not a weighted
#' mean of the stratum estimands.
#'
#' @param x A [pollin_fit()] object.
#' @param n Number of draws.
#' @param method `"simulation"` (the default for a glmmTMB fit),
#'   `"bootstrap"`, or `"posterior"` (for a flexyBayes fit).
#' @param response The response to draw. Required when the fit holds more than
#'   one response.
#' @param seed Optional seed for the draws. The session's random number state
#'   is restored afterwards. Posterior draws with a seed set both R's
#'   generator and INLA's, and are reproducible.
#'
#' @return An object of class `pollin_draws`: a list holding
#'   * `high`, `low`: numeric matrices of draws of the arm means, one row per
#'     draw and one column per level of the fit's `by` column, or a single
#'     column `"overall"` when the fit has no `by`;
#'   * `estimate`: the matrix of the estimand on its scale and basis, formed
#'     from `high` and `low` draw by draw;
#'   * `estimand`, `response`, `by`, `method`, `engine`, `n` (the number of
#'     draws), `seed`, and `population` (the [pollin_population()] object,
#'     or `NULL` for the trial population).
#' @seealso [pollin_effects()] for estimates and intervals;
#'   [pollin_population()] for a population of strata.
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
#' dr <- pollin_draws(pollin_fit(est), n = 2000L, seed = 1L)
#' dr
#' summary(dr)
#' }
pollin_draws <- function(
  x,
  n = 4000L,
  method = c("simulation", "bootstrap", "posterior"),
  response = NULL,
  seed = NULL
) {
  .check_class(x, "pollin_fit")
  if (missing(method)) {
    method <- if (x$engine == "flexyBayes") "posterior" else "simulation"
  }
  .check_choice(method, "method", c("simulation", "bootstrap", "posterior"))
  .draws_check_method(method, x$engine)
  .draws_check_n(n)
  .check_seed(seed)
  response <- .draws_response(response, names(x$fits))
  est <- x$estimand
  population <- .draws_population(est, x$by)

  m <- if (method == "posterior") {
    .flexy_arm_draws(x, response, n, seed, x$by)
  } else if (method == "simulation") {
    .draws_simulation(x, response, n, seed)
  } else {
    bt <- .boot_arm_means(x, response, x$by, n, seed)
    ok <- !is.na(bt$high[, 1L])
    list(
      high = bt$high[ok, , drop = FALSE],
      low = bt$low[ok, , drop = FALSE]
    )
  }
  if (is.null(x$by)) {
    colnames(m$high) <- colnames(m$low) <- "overall"
  }
  .new_pollin_draws(
    high = m$high,
    low = m$low,
    estimand = est,
    response = response,
    by = x$by,
    method = method,
    engine = x$engine,
    n = nrow(m$high),
    seed = seed,
    population = population
  )
}

#' @export
print.pollin_draws <- function(x, ...) {
  e <- x$estimand
  label <- if (e$family == "deficit") "deficit (E1)" else "dependence (E2)"
  cat("<pollin_draws>\n")
  cat("  family   :", label, "\n")
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
  cat("  response :", x$response, "\n")
  cat("  method   :", .draws_method_label(x$method), "\n")
  cat("  engine   :", x$engine, "\n")
  cat("  draws    :", x$n, "\n")
  cat("  columns  :", paste(colnames(x$estimate), collapse = ", "), "\n")
  if (x$method == "simulation") {
    cat("  note     : simulation draws cover fixed-effect uncertainty only\n")
  }
  invisible(x)
}

#' Summarise draws of the estimand
#'
#' The median and percentile interval of the estimand draws, per column.
#' `as.matrix()` returns the matrix of estimand draws.
#'
#' @param object,x A [pollin_draws()] object.
#' @param conf_level Level of the percentile interval.
#' @param ... Not used.
#'
#' @return A data frame with one row per column of the draws and the columns
#'   `response`, `level` (the column name: a `by` level, `"overall"` or
#'   `"population"`), `median`, `conf.low` and `conf.high` (the percentile
#'   interval), `family`, `scale`, `basis`, `method` and `n`.
#'   `as.matrix()` returns the `estimate` matrix of the draws.
#' @export
summary.pollin_draws <- function(object, conf_level = 0.95, ...) {
  .check_conf_level(conf_level)
  a <- (1 - conf_level) / 2
  q <- apply(
    object$estimate,
    2L,
    stats::quantile,
    probs = c(0.5, a, 1 - a),
    names = FALSE
  )
  e <- object$estimand
  data.frame(
    response = object$response,
    level = colnames(object$estimate),
    median = q[1L, ],
    conf.low = q[2L, ],
    conf.high = q[3L, ],
    family = e$family,
    scale = e$scale,
    basis = e$basis,
    method = object$method,
    n = object$n,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' @rdname summary.pollin_draws
#' @export
as.matrix.pollin_draws <- function(x, ...) {
  x$estimate
}

# Constructor -----------------------------------------------------------------

#' Build a `pollin_draws` object from draws of the arm means. Every engine
#' (simulation, bootstrap, or a Bayesian posterior) returns its draws through
#' this constructor: `high` and `low` are numeric matrices with one row per
#' draw and named columns (the `by` levels, or `"overall"`). The estimand is
#' formed per draw on its scale and basis, and a `"population"` column is
#' added when `population` is a `pollin_population`.
#' @noRd
#' @keywords internal
.new_pollin_draws <- function(
  high,
  low,
  estimand,
  response,
  by,
  method,
  engine,
  n,
  seed,
  population
) {
  .check_class(estimand, "pollin_estimand", arg = "estimand")
  .draws_check_matrix(high, "high")
  .draws_check_matrix(low, "low")
  if (
    !identical(dim(high), dim(low)) || !identical(colnames(high), colnames(low))
  ) {
    stop(
      "`high` and `low` must have the same dimensions and column names.",
      call. = FALSE
    )
  }
  if (!identical(as.integer(n), nrow(high))) {
    stop("`n` must equal the number of rows of `high`.", call. = FALSE)
  }
  if (!is.character(response) || length(response) != 1L) {
    stop("`response` must be one string.", call. = FALSE)
  }
  if (!is.character(method) || length(method) != 1L) {
    stop("`method` must be one string.", call. = FALSE)
  }
  if (!is.character(engine) || length(engine) != 1L) {
    stop("`engine` must be one string.", call. = FALSE)
  }
  if (!is.null(population)) {
    m <- .population_means(high, low, population)
    high <- m$high
    low <- m$low
  }
  structure(
    list(
      high = high,
      low = low,
      estimate = .effect_from_means(
        list(high = high, low = low),
        estimand$scale,
        estimand$basis
      ),
      estimand = estimand,
      response = response,
      by = by,
      method = method,
      engine = engine,
      n = as.integer(n),
      seed = seed,
      population = population
    ),
    class = "pollin_draws"
  )
}

# Simulation ------------------------------------------------------------------

#' Arm-mean draws from a multivariate normal draw of the fixed effects. The
#' random-effect part of each row's linear predictor, `predict(type = "link")`
#' less the fixed part, is held at its estimate, so at the estimated fixed
#' effects the arm means equal those of `.arm_means()`.
#' @noRd
#' @keywords internal
.draws_simulation <- function(fit, response, n, seed) {
  e <- fit$estimand
  model <- fit$fits[[response]]
  d <- fit$data[[response]]
  beta <- glmmTMB::fixef(model)$cond
  v <- stats::vcov(model, full = FALSE)$cond
  r <- tryCatch(chol(v), error = function(err) NULL)
  if (is.null(r)) {
    stop(
      "The covariance matrix of the fixed effects of '",
      response,
      "' is not positive definite; use `method = \"bootstrap\"`.",
      call. = FALSE
    )
  }
  z <- matrix(.draws_rnorm(n * length(beta), seed), nrow = n)
  b <- z %*% r + rep(beta, each = n)

  arm <- function(level) {
    p <- .draws_arm_linear(model, d, fit$treatment, level, beta)
    mu <- stats::family(model)$linkinv(p$x %*% t(b) + p$re)
    g <- if (is.null(fit$by)) rep("all", nrow(d)) else as.character(d[[fit$by]])
    t(rowsum(mu, g) / as.vector(table(g)))
  }
  list(high = arm(e$contrast[["high"]]), low = arm(e$contrast[["low"]]))
}

#' Fixed-effect design matrix of the model data with every row set to one
#' treatment level, and the random-effect offset of each row's linear
#' predictor
#' @noRd
#' @keywords internal
.draws_arm_linear <- function(model, d, treatment, level, beta) {
  nd <- d
  nd[[treatment]] <- factor(
    rep(level, nrow(d)),
    levels = levels(d[[treatment]])
  )
  tt <- stats::delete.response(stats::terms(model))
  x <- stats::model.matrix(tt, nd)
  missing_cols <- setdiff(names(beta), colnames(x))
  if (length(missing_cols)) {
    stop(
      "Fixed effects not found in the model matrix: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  x <- x[, names(beta), drop = FALSE]
  eta <- stats::predict(model, newdata = nd, type = "link")
  list(x = x, re = as.vector(eta - x %*% beta))
}

#' `n` standard normal deviates, from `seed` when given, restoring the
#' session's random number state afterwards
#' @noRd
#' @keywords internal
.draws_rnorm <- function(n, seed) {
  if (!is.null(seed)) {
    if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      stats::runif(1L)
    }
    old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = globalenv()))
    set.seed(seed)
  }
  stats::rnorm(n)
}

# Internal ---------------------------------------------------------------------

#' The method in words, for printing
#' @noRd
#' @keywords internal
.draws_method_label <- function(method) {
  switch(
    method,
    simulation = "simulation of the fixed effects",
    bootstrap = "parametric bootstrap refits",
    posterior = "joint posterior draws (INLA)",
    method
  )
}

#' The population of the estimand, or NULL for the trial; stops when the fit
#' was not made with the population's `by` column
#' @noRd
#' @keywords internal
.draws_population <- function(estimand, fit_by) {
  pop <- estimand$population
  if (!inherits(pop, "pollin_population")) {
    return(NULL)
  }
  .population_check_fit_by(pop, fit_by)
  pop
}

#' Stop unless the method suits the engine: posterior draws for a
#' flexyBayes fit, simulation or bootstrap for a glmmTMB fit
#' @noRd
#' @keywords internal
.draws_check_method <- function(method, engine) {
  if (engine == "flexyBayes" && method != "posterior") {
    stop(
      "A flexyBayes fit is drawn from its posterior; `method = \"",
      method,
      "\"` applies to glmmTMB fits.",
      call. = FALSE
    )
  }
  if (engine != "flexyBayes" && method == "posterior") {
    stop(
      "`method = \"posterior\"` needs a fit with `engine = \"flexyBayes\"`.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Stop unless `n` is one whole number of at least 2
#' @noRd
#' @keywords internal
.draws_check_n <- function(n) {
  if (
    !is.numeric(n) ||
      length(n) != 1L ||
      !is.finite(n) ||
      n < 2 ||
      n != round(n)
  ) {
    stop("`n` must be a whole number of at least 2.", call. = FALSE)
  }
  invisible(TRUE)
}

#' The response to draw: the only one, or the one named
#' @noRd
#' @keywords internal
.draws_response <- function(response, responses) {
  if (is.null(response)) {
    if (length(responses) > 1L) {
      stop(
        "The fit holds several responses (",
        paste(responses, collapse = ", "),
        "); name one with `response =`.",
        call. = FALSE
      )
    }
    return(responses)
  }
  if (
    !is.character(response) ||
      length(response) != 1L ||
      !response %in% responses
  ) {
    stop(
      "`response` must be one of the fitted responses: ",
      paste(responses, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  response
}

#' Stop unless `x` is a numeric matrix with at least one row and named columns
#' @noRd
#' @keywords internal
.draws_check_matrix <- function(x, arg) {
  if (
    !is.matrix(x) ||
      !is.numeric(x) ||
      nrow(x) < 1L ||
      is.null(colnames(x)) ||
      anyDuplicated(colnames(x))
  ) {
    stop(
      "`",
      arg,
      "` must be a numeric matrix with one row per draw and uniquely named ",
      "columns.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
