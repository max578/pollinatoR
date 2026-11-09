# Model fit (layer L3) ---------------------------------------------------------

#' Fit the outcome models for an estimand
#'
#' Fits one generalised linear mixed model per response in the estimand's data,
#' with the pollination treatment as the fixed effect and the random terms of
#' the design: `response ~ treatment + <design random terms>`. With `by`, the
#' fixed part is `treatment * by`, so the treatment effect has its own fixed
#' estimate at each level of `by`. The treatment is releveled so the
#' estimand's `low` arm is the reference level. Every treatment level in the
#' data is kept in the fit; the contrast is read from the fitted model by
#' [pollin_effects()].
#'
#' For an estimand with `effect = "direct"`, the fixed part adds the
#' estimand's adjustment columns (its mediators, and with a graph the rest of
#' its direct-effect set) as main effects, `response ~ treatment + <adjustment>
#' + <random terms>`. The arm means then keep each row's observed adjustment
#' values, so the contrast is the controlled direct effect averaged over the
#' observed mediator values. With no treatment-by-mediator interaction in the
#' model, that is the effect of the treatment at any fixed mediator value.
#'
#' Rows with a missing value in the response, the treatment, a unit column,
#' the denominator, the `by` column or a direct effect's adjustment column are
#' dropped for that response's model.
#'
#' When the data declare a denominator (`pollin_data(denominator = )`), each
#' response is a count of successes out of that many trials and is fitted as
#' binomial with a logit link, `cbind(response, denominator - response)` on
#' glmmTMB. The arm means read by [pollin_effects()] and [pollin_draws()]
#' are then expected proportions, averaged over the rows with equal weight.
#'
#' @param estimand A [pollin_estimand()] object.
#' @param engine The fitting engine: `"glmmTMB"` (the default) or
#'   `"flexyBayes"`, a Bayesian fit on INLA through flexyBayes
#'   (see Details). On designs with a random treatment slope and few units
#'   (12 sites tested), its 95% posterior intervals cover about 91% to 93%,
#'   as an exact posterior does on those designs.
#' @param family `NULL` (the default) chooses a family per response from its
#'   values (count, proportion in (0, 1), positive continuous, or continuous
#'   with zeros), or `binomial()` when the data declare a denominator. A
#'   single `glmmTMB` family object applies to every response; a named list
#'   gives one per response, and responses it does not name are chosen as
#'   for `NULL`. With a denominator, a family other than binomial is an
#'   error.
#' @param by Optional name of a column of the data (for example variety or
#'   region). The treatment effect is fitted separately at each of its levels,
#'   as an interaction, and [pollin_effects()] then reports one contrast per
#'   level. The column is treated as a factor. An estimand whose population
#'   is a [pollin_population()] must be fitted with `by` set to the
#'   population's `by` column.
#' @param pooled Optional name of a column of the data that is not a unit
#'   column (for example variety). Its levels get a random intercept and a
#'   random treatment slope, so the treatment effect at each level is pooled
#'   towards the overall effect; [pollin_cate()] reads the per-level effects.
#'   It must differ from `by`. A slope on a unit column is declared with
#'   `pollin_design(slope = )` instead. Rows with a missing value in the
#'   column are dropped.
#' @param ... Passed to [glmmTMB::glmmTMB()], or to `flexyBayes::flexybayes()`
#'   on the flexyBayes engine (for example `prior` or `seed`).
#'
#' @details
#' With `engine = "flexyBayes"` each response is fitted by
#' `flexyBayes::flexybayes(backend = "inla")`, which needs the flexyBayes and
#' INLA packages. The fixed part is that of the glmmTMB engine. Each unit
#' level is a random intercept; a nested level enters as one grouping column
#' `pollin_<outer>_<inner>` added to the model data. flexyBayes fits a
#' correlated random slope on brms only, so on INLA the treatment slope of
#' the design is a random intercept per unit plus an independent effect per
#' unit-by-arm cell (column `pollin_<unit>_<treatment>`): equal variances
#' across arms and a non-negative correlation between them, the exchangeable
#' form of the slope. A `pooled` column enters the same way. The families
#' map as `nbinom2()` to `"negative_binomial"`, `poisson()` to `"poisson"`,
#' `binomial()` to `"binomial"`, `Gamma(link = "log")` to `"gamma"` and
#' `beta_family()` to `"beta"`; any other family is refused. With a
#' denominator, flexyBayes takes a binomial response as one 0/1 row per
#' trial, so each row is expanded to as many rows as its denominator, its
#' successes coded 1 and the rest 0; the arm means are still averaged over
#' the original rows. An expansion beyond two million rows is an error, and
#' such data are fitted with the glmmTMB engine. Priors are
#' passed through `...` as `flexybayes()` takes them; the seed of the
#' posterior draws is set in [pollin_draws()]. The backend is INLA only, since the posterior draws of a brms fit are not in
#' the flexyBayes extension API version 1. A fit whose convergence verdict
#' (`flexyBayes::fb_convergence()`) is `FALSE` is an error, and one with no
#' verdict a warning. Estimates and intervals of a flexyBayes fit come from
#' its posterior draws; see [pollin_draws()].
#'
#' @return An object of class `pollin_fit`: a list holding the `estimand`, the
#'   `engine`, the fitted models `fits` (named by response), the model `data`
#'   for each response, the right-hand-side terms `rhs`, the `families` (as
#'   glmmTMB family objects on both engines), the `by` column (or `NULL`),
#'   the `denominator` column (or `NULL`), and the further fitting arguments
#'   `args`.
#' @seealso [pollin_effects()] for estimates and intervals on the estimand's
#'   scale.
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
#' pollin_fit(est)
#' }
pollin_fit <- function(
  estimand,
  engine = "glmmTMB",
  family = NULL,
  by = NULL,
  pooled = NULL,
  ...
) {
  .check_class(estimand, "pollin_estimand", arg = "estimand")
  .check_engine(engine)
  if (inherits(estimand$population, "pollin_population")) {
    .population_check_fit_by(estimand$population, by)
  }
  inputs <- .pooled_inputs(.estimand_inputs(estimand), pooled, by)
  .check_family_arg(family, inputs$responses)
  if (engine == "flexyBayes") {
    return(.flexy_new_fit(estimand, inputs, family, ..., by = by))
  }
  .new_fit(estimand, inputs, inputs$responses, family, by = by, ...)
}

#' @export
print.pollin_fit <- function(x, ...) {
  e <- x$estimand
  cat("<pollin_fit>\n")
  cat("  engine    :", x$engine, "\n")
  cat("  family    :", e$family, "\n")
  cat("  contrast  :", e$contrast[["high"]], "vs", e$contrast[["low"]], "\n")
  cat("  model     : <response> ~", paste(x$rhs, collapse = " + "), "\n")
  if (!is.null(x$by)) {
    cat("  by        :", x$by, "\n")
  }
  if (!is.null(x$denominator)) {
    cat("  trials    :", x$denominator, "(binomial: successes out of trials)\n")
  }
  if (!is.null(x$pooled)) {
    cat("  pooled    :", x$pooled, "(random treatment slope)\n")
  }
  for (r in names(x$fits)) {
    cat(sprintf(
      "  %-10s: %s, %d rows\n",
      r,
      x$families[[r]]$family,
      nrow(x$data[[r]])
    ))
  }
  invisible(x)
}

# Fit inputs ------------------------------------------------------------------

#' Fit inputs shared by every model-based verb: the data with the treatment
#' releveled on `low`, the responses, the random terms, the slope term, the
#' grouping for the within-between split and the denominator
#' @noRd
#' @keywords internal
.fit_inputs <- function(
  data,
  responses,
  treatment,
  low,
  unit,
  random,
  slope_key,
  wb_group,
  exposures,
  adjust,
  denominator = NULL
) {
  data[[treatment]] <- stats::relevel(factor(data[[treatment]]), ref = low)
  list(
    data = data,
    responses = responses,
    treatment = treatment,
    low = low,
    unit = unit,
    random = random,
    slope_key = slope_key,
    wb_group = wb_group,
    exposures = exposures,
    adjust = adjust,
    denominator = denominator
  )
}

#' Fit inputs from a `pollin_design`, with `low` as the reference level
#' @noRd
#' @keywords internal
.design_inputs <- function(design, low) {
  pd <- design$data
  r <- pd$roles
  unit <- r$unit
  depth <- if (is.null(design$slope)) 1L else match(design$slope, unit)
  key <- unit[seq_len(depth)]
  .fit_inputs(
    data = pd$data,
    responses = r$response,
    treatment = r$treatment,
    low = low,
    unit = unit,
    random = design$random,
    slope_key = if (is.null(design$slope)) NULL else paste(key, collapse = ":"),
    wb_group = interaction(pd$data[key], drop = TRUE),
    exposures = r$exposures,
    adjust = r$adjust,
    denominator = r$denominator
  )
}

#' Fit inputs with a random treatment slope on the levels of `pooled`: rows
#' missing it dropped, the column a factor, and the slope term added to the
#' random terms
#' @noRd
#' @keywords internal
.pooled_inputs <- function(inputs, pooled, by) {
  if (is.null(pooled)) {
    return(inputs)
  }
  d <- inputs$data
  if (
    !is.character(pooled) ||
      length(pooled) != 1L ||
      !pooled %in% names(d) ||
      pooled %in% c(inputs$treatment, inputs$responses)
  ) {
    stop(
      "`pooled` must name one column of the data other than the treatment ",
      "and the responses.",
      call. = FALSE
    )
  }
  if (pooled %in% inputs$unit) {
    stop(
      "`pooled` names the unit column '",
      pooled,
      "'; declare a slope on a unit with `pollin_design(slope = )`.",
      call. = FALSE
    )
  }
  if (identical(pooled, by)) {
    stop(
      "`pooled` and `by` name the same column; a level's effect is either ",
      "fixed (`by`) or pooled (`pooled`).",
      call. = FALSE
    )
  }
  d <- d[!is.na(d[[pooled]]), , drop = FALSE]
  d[[pooled]] <- factor(d[[pooled]])
  inputs$data <- d
  inputs$random <- c(
    inputs$random,
    sprintf("(%s | %s)", inputs$treatment, pooled)
  )
  inputs$pooled <- pooled
  inputs
}

#' Fit inputs from a `pollin_estimand`: the design's inputs with the estimand's
#' `low` arm as the reference level and, for a direct effect, the adjustment
#' columns as `covariates`
#' @noRd
#' @keywords internal
.estimand_inputs <- function(estimand) {
  inputs <- .design_inputs(estimand$design, estimand$contrast[["low"]])
  if (identical(estimand$effect, "direct")) {
    inputs$covariates <- estimand$adjustment
  }
  inputs
}

# Fitting ---------------------------------------------------------------------

#' Fit one `glmmTMB` model; the formula is built in this frame and its
#' environment holds `d`. `response` is a column name or a left-hand side
#' from `.response_lhs()`.
#' @noRd
#' @keywords internal
.fit_glmm <- function(response, rhs, d, family, ...) {
  form <- stats::reformulate(rhs, response = response)
  glmmTMB::glmmTMB(form, data = d, family = family, ...)
}

#' Fit `treatment + random` for each response, or `treatment * by + random`
#' when `by` is given. A response whose model fails is dropped with a warning;
#' no fitted response is an error.
#' @noRd
#' @keywords internal
.fit_responses <- function(inputs, responses, family, ..., by = NULL) {
  fixed <- if (is.null(by)) {
    inputs$treatment
  } else {
    paste(inputs$treatment, "*", by)
  }
  rhs <- c(fixed, inputs$covariates, inputs$random)
  fits <- list()
  data <- list()
  families <- list()
  den <- inputs$denominator
  for (r in responses) {
    d <- .response_rows(inputs, r, by)
    fam <- .family_for(family, r, d[[r]], den)
    fit <- tryCatch(
      .fit_glmm(.response_lhs(r, den), rhs, d, fam, ...),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      warning(
        "Outcome '",
        r,
        "' failed to fit: ",
        conditionMessage(fit),
        call. = FALSE
      )
      next
    }
    fits[[r]] <- fit
    data[[r]] <- d
    families[[r]] <- fam
  }
  if (length(fits) == 0L) {
    stop("No outcome could be fitted.", call. = FALSE)
  }
  list(fits = fits, data = data, families = families, rhs = rhs)
}

#' Build a `pollin_fit` for the given responses
#' @noRd
#' @keywords internal
.new_fit <- function(estimand, inputs, responses, family, ..., by = NULL) {
  .check_fit_by(by, inputs$data, inputs$treatment)
  f <- .fit_responses(inputs, responses, family, ..., by = by)
  structure(
    list(
      estimand = estimand,
      engine = "glmmTMB",
      fits = f$fits,
      data = f$data,
      rhs = f$rhs,
      families = f$families,
      treatment = inputs$treatment,
      by = by,
      slope_key = inputs$slope_key,
      pooled = inputs$pooled,
      denominator = inputs$denominator,
      args = list(...)
    ),
    class = "pollin_fit"
  )
}

#' The rows of the inputs' data fitted for one response: complete in the
#' response, the treatment, the units, the denominator and `by`, with `by` a
#' factor, and in a direct effect's covariates; with a denominator, the counts
#' are checked against it
#' @noRd
#' @keywords internal
.response_rows <- function(inputs, response, by) {
  d <- inputs$data
  den <- inputs$denominator
  cols <- c(
    response,
    inputs$treatment,
    inputs$unit,
    den,
    by,
    inputs$covariates
  )
  keep <- stats::complete.cases(d[cols])
  if (!all(keep)) {
    d <- d[keep, , drop = FALSE]
  }
  if (!is.null(by)) {
    d[[by]] <- factor(d[[by]])
  }
  .check_trials(d[[response]], d[[den]], response, den)
  d
}

#' Left-hand side of the model formula: the response column, or
#' `cbind(successes, trials - successes)` when there is a denominator
#' @noRd
#' @keywords internal
.response_lhs <- function(response, denominator) {
  if (is.null(denominator)) {
    return(response)
  }
  y <- as.name(response)
  deparse1(call("cbind", y, call("-", as.name(denominator), y)))
}

#' The family for one response: from a named list, a single family, or the
#' heuristic in `.pick_family()`; binomial when there is a denominator
#' @noRd
#' @keywords internal
.family_for <- function(family, response, y, denominator = NULL) {
  if (!is.null(denominator)) {
    return(.trials_family(family, response, denominator))
  }
  if (inherits(family, "family")) {
    return(family)
  }
  if (!is.null(family[[response]])) {
    return(family[[response]])
  }
  .pick_family(y)
}

#' Standard deviation of the random treatment slope on the `key` term: the
#' `coef` column, or the second column when `coef` is NULL. NA when the term
#' carries no slope.
#' @noRd
#' @keywords internal
.slope_sd <- function(fit, key, coef = NULL) {
  if (is.null(key)) {
    return(NA_real_)
  }
  vc <- tryCatch(glmmTMB::VarCorr(fit)$cond[[key]], error = function(e) NULL)
  if (is.null(vc) || nrow(vc) < 2L) {
    return(NA_real_)
  }
  sds <- sqrt(diag(vc))
  unname(if (is.null(coef)) sds[2L] else sds[coef])
}

#' The family of a response counted out of a denominator: binomial(logit)
#' unless `family` gives a binomial family for it; any other family is an
#' error
#' @noRd
#' @keywords internal
.trials_family <- function(family, response, denominator) {
  fam <- if (inherits(family, "family")) family else family[[response]]
  if (is.null(fam)) {
    return(stats::binomial())
  }
  if (!identical(fam$family, "binomial")) {
    stop(
      "`family` gives '",
      fam$family,
      "' for '",
      response,
      "', but the data declare the denominator '",
      denominator,
      "', so '",
      response,
      "' is a count out of trials and is fitted as binomial. Leave `family` ",
      "unset or give `binomial()`; to model the count itself, drop ",
      "`denominator` from `pollin_data()`.",
      call. = FALSE
    )
  }
  fam
}

# Internal checks -------------------------------------------------------------

#' Stop unless the successes `y` and trials `n` of the fitted rows are whole
#' numbers with `0 <= y <= n` and `n > 0`
#' @noRd
#' @keywords internal
.check_trials <- function(y, n, response, denominator) {
  if (is.null(denominator)) {
    return(invisible(TRUE))
  }
  whole <- function(x) all(abs(x - round(x)) < 1e-8)
  if (!whole(y)) {
    stop(
      "'",
      response,
      "' is counted out of the denominator '",
      denominator,
      "' and must be whole numbers; it has fractional values.",
      call. = FALSE
    )
  }
  if (!whole(n)) {
    stop(
      "The denominator '",
      denominator,
      "' must be whole numbers of trials; it has fractional values.",
      call. = FALSE
    )
  }
  zero <- sum(n == 0)
  if (zero > 0L) {
    stop(
      "'",
      response,
      "' has ",
      .n_noun(zero, "row"),
      " with the denominator '",
      denominator,
      "' equal to 0, which give no proportion; drop them or set the ",
      "denominator to NA.",
      call. = FALSE
    )
  }
  out <- sum(y < 0 | y > n)
  if (out > 0L) {
    stop(
      "'",
      response,
      "' has ",
      .n_noun(out, "row"),
      " outside 0 to the denominator '",
      denominator,
      "'.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Stop unless `by` is NULL or one column of the data other than the
#' treatment
#' @noRd
#' @keywords internal
.check_fit_by <- function(by, data, treatment) {
  if (is.null(by)) {
    return(invisible(TRUE))
  }
  if (
    !is.character(by) ||
      length(by) != 1L ||
      !by %in% names(data) ||
      identical(by, treatment)
  ) {
    stop(
      "`by` must name one column of the data other than the treatment.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Stop unless `engine` is "glmmTMB" or "flexyBayes"
#' @noRd
#' @keywords internal
.check_engine <- function(engine) {
  .check_choice(engine, "engine", c("glmmTMB", "flexyBayes"))
}

#' Check `family`: NULL, one family object, or a list of family objects named
#' by response
#' @noRd
#' @keywords internal
.check_family_arg <- function(family, responses) {
  if (
    is.null(family) ||
      inherits(family, "family") ||
      (is.list(family) && length(family) == 0L)
  ) {
    return(invisible(TRUE))
  }
  ok <- is.list(family) &&
    !is.null(names(family)) &&
    all(nzchar(names(family))) &&
    all(vapply(family, inherits, logical(1), what = "family"))
  if (!ok) {
    stop(
      "`family` must be NULL, a family object, or a list of family objects ",
      "named by response.",
      call. = FALSE
    )
  }
  unknown <- setdiff(names(family), responses)
  if (length(unknown)) {
    stop(
      "`family` names responses not in the data: ",
      paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}
