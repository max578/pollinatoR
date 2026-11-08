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
#' Rows with a missing value in the response, the treatment, a unit column or
#' the `by` column are dropped for that response's model.
#'
#' @param estimand A [pollin_estimand()] object.
#' @param engine The fitting engine. Only `"glmmTMB"` is available;
#'   `"flexyBayes"` arrives in version 0.2.0.
#' @param family `NULL` (the default) chooses a family per response from its
#'   values (count, proportion in (0, 1), positive continuous, or continuous
#'   with zeros). A single `glmmTMB` family object applies to every response; a
#'   named list gives one per response, and responses it does not name are
#'   chosen from their values.
#' @param by Optional name of a column of the data (for example variety or
#'   region). The treatment effect is fitted separately at each of its levels,
#'   as an interaction, and [pollin_effects()] then reports one contrast per
#'   level. The column is treated as a factor.
#' @param ... Passed to [glmmTMB::glmmTMB()].
#'
#' @return An object of class `pollin_fit`: a list holding the `estimand`, the
#'   `engine`, the fitted models `fits` (named by response), the model `data`
#'   for each response, the right-hand-side terms `rhs`, the `families`, the
#'   `by` column (or `NULL`), and the further fitting arguments `args`.
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
  ...
) {
  .check_class(estimand, "pollin_estimand", arg = "estimand")
  .check_engine(engine)
  inputs <- .estimand_inputs(estimand)
  .check_family_arg(family, inputs$responses)
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
#' releveled on `low`, the responses, the random terms, the slope term and the
#' grouping for the within-between split
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
  adjust
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
    adjust = adjust
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
    adjust = r$adjust
  )
}

#' Fit inputs from a `pollin_estimand`: the design's inputs with the estimand's
#' `low` arm as the reference level
#' @noRd
#' @keywords internal
.estimand_inputs <- function(estimand) {
  .design_inputs(estimand$design, estimand$contrast[["low"]])
}

#' Fit inputs from a deprecated `pollin_spec`, without the `pollin_data` checks
#' @noRd
#' @keywords internal
.spec_inputs <- function(spec) {
  .fit_inputs(
    data = spec$data,
    responses = spec$outcomes,
    treatment = spec$treatment,
    low = spec$control,
    unit = spec$random,
    random = .random_terms(spec),
    slope_key = spec$latent_unit,
    wb_group = spec$data[[spec$latent_unit]],
    exposures = spec$exposures,
    adjust = spec$adjust
  )
}

# Fitting ---------------------------------------------------------------------

#' Fit one `glmmTMB` model; the formula is built in this frame and its
#' environment holds `d`
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
  rhs <- c(fixed, inputs$random)
  fits <- list()
  data <- list()
  families <- list()
  for (r in responses) {
    d <- inputs$data
    keep <- stats::complete.cases(d[c(r, inputs$treatment, inputs$unit, by)])
    if (!all(keep)) {
      d <- d[keep, , drop = FALSE]
    }
    if (!is.null(by)) {
      d[[by]] <- factor(d[[by]])
    }
    fam <- .family_for(family, r, d[[r]])
    fit <- tryCatch(
      .fit_glmm(r, rhs, d, fam, ...),
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
      args = list(...)
    ),
    class = "pollin_fit"
  )
}

#' The family for one response: from a named list, a single family, or the
#' heuristic in `.pick_family()`
#' @noRd
#' @keywords internal
.family_for <- function(family, response, y) {
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

# Internal checks -------------------------------------------------------------

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

#' Stop unless `engine` is "glmmTMB"; "flexyBayes" points to 0.2.0
#' @noRd
#' @keywords internal
.check_engine <- function(engine) {
  if (identical(engine, "flexyBayes")) {
    stop(
      "The `flexyBayes` engine arrives in pollinatoR 0.2.0; ",
      "use `engine = \"glmmTMB\"`.",
      call. = FALSE
    )
  }
  if (!identical(engine, "glmmTMB")) {
    stop("`engine` must be \"glmmTMB\".", call. = FALSE)
  }
  invisible(TRUE)
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
