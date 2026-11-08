# E2 on the difference scale: the pollinator contribution ---------------------

#' Estimate the pollinator contribution per outcome
#'
#' Fits, for each outcome, a generalised linear mixed model with the
#' pollination treatment as the fixed effect of interest and the random
#' structure of the design ([pollin_fit()]), then reports the contribution as a
#' model-based average treatment effect, the reference-minus-exclusion
#' difference in expected response, with propagated uncertainty, via
#' G-computation ([pollin_effects()]). A random treatment slope, where the
#' design carries one, absorbs unobserved between-unit heterogeneity in the
#' contribution (for example an unrecorded genotype). Family is chosen per
#' outcome from its values unless supplied.
#'
#' The treatment effect is causal by design here -- the exclusion contrast is an
#' experimental manipulation -- so no covariate adjustment is applied; see
#' [pollin_moderation()] for the observational "which insects" question.
#'
#' A [pollin_spec()] is still accepted and is fitted as before; `pollin_spec()`
#' is deprecated and is removed at version 0.2.0.
#'
#' @param x A [pollin_estimand()] object with `family = "dependence"` and
#'   `scale = "difference"`, or a deprecated [pollin_spec()] object.
#' @param outcomes Outcomes to model. `NULL` (the default) models every
#'   response.
#' @param by Optional name of a column for conditional contrasts, one per
#'   level (for example a region).
#' @param families Optional named list of `glmmTMB` family objects, one per
#'   outcome; outcomes it does not name are chosen from their values.
#' @param ... Passed to [glmmTMB::glmmTMB()].
#'
#' @return An object of class `pollin_contribution`: a list with `table` (a tidy
#'   data frame of contrasts with FDR-adjusted p-values and the standard
#'   deviation of the random treatment slope) and `fits` (the fitted models,
#'   for diagnostics). From an estimand, `table` carries the columns of
#'   [pollin_effects()], and the list also holds the `estimand`.
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' site <- rep(paste0("s", seq_len(12L)), each = 20L)
#' trt <- rep(rep(c("Cage", "Open"), each = 10L), 12L)
#' y <- rpois(length(trt), exp(1.7 + 0.35 * (trt == "Open")))
#' d <- data.frame(y = y, trt = trt, Site = site)
#' pd <- pollin_data(d, "y", "trt", unit = "Site")
#' ds <- pollin_design(pd, c(excluded = "Cage", open = "Open"), "multi_site")
#' est <- pollin_estimand(ds, "dependence", "difference", "loss")
#' pollin_contribution(est)
#' }
pollin_contribution <- function(
  x,
  outcomes = NULL,
  by = NULL,
  families = list(),
  ...
) {
  if (inherits(x, "pollin_spec")) {
    return(.contribution_spec(x, outcomes, by, families, ...))
  }
  .check_class(x, "pollin_estimand")
  .check_contribution_estimand(x)
  inputs <- .estimand_inputs(x)
  if (is.null(outcomes)) {
    outcomes <- inputs$responses
  }
  .check_outcomes(outcomes, inputs$responses)
  .check_family_arg(families, inputs$responses)

  fit <- .new_fit(x, inputs, outcomes, families, ..., by = by)
  tab <- pollin_effects(fit, by = by)$table
  high <- x$contrast[["high"]]
  slope_sd <- vapply(
    fit$fits,
    .slope_sd,
    numeric(1),
    key = inputs$slope_key,
    coef = paste0(inputs$treatment, high)
  )
  tab$slope_sd <- unname(slope_sd[tab$response])
  tab$p.fdr <- pollin_fdr(tab$p.value)
  structure(
    list(table = tab, fits = fit$fits, estimand = x),
    class = "pollin_contribution"
  )
}

#' @export
print.pollin_contribution <- function(x, ...) {
  cat("<pollin_contribution>  (open - exclusion, model-based ATE)\n")
  print(x$table, row.names = FALSE)
  invisible(x)
}

# Deprecated pollin_spec path -------------------------------------------------

#' Contribution from a `pollin_spec`: the 0.0.0.9000 table, column for column
#' @noRd
#' @keywords internal
.contribution_spec <- function(spec, outcomes, by, families, ...) {
  if (is.null(outcomes)) {
    outcomes <- spec$outcomes
  }
  inputs <- .spec_inputs(spec)
  f <- .fit_responses(inputs, outcomes, families, ...)

  # Contrasts per fitted outcome ------------------------------------------------
  rows <- list()
  fits <- list()
  for (oc in names(f$fits)) {
    fit <- f$fits[[oc]]
    cmp <- tryCatch(
      if (is.null(by)) {
        marginaleffects::avg_comparisons(fit, variables = spec$treatment)
      } else {
        marginaleffects::avg_comparisons(
          fit,
          variables = spec$treatment,
          by = by
        )
      },
      error = function(e) e
    )
    if (inherits(cmp, "error")) {
      warning(
        "Outcome '",
        oc,
        "' contrast failed: ",
        conditionMessage(cmp),
        call. = FALSE
      )
      next
    }
    cmp$outcome <- oc
    cmp$family <- f$families[[oc]]$family
    cmp$slope_sd <- .slope_sd(fit, inputs$slope_key)
    rows[[oc]] <- as.data.frame(cmp)
    fits[[oc]] <- fit
  }
  if (length(rows) == 0L) {
    stop("No outcome could be fitted.", call. = FALSE)
  }

  # Tidy table ----------------------------------------------------------------
  keep <- c(
    "outcome",
    "term",
    "contrast",
    intersect(
      c(if (!is.null(by)) by),
      names(rows[[1]])
    ),
    "estimate",
    "std.error",
    "p.value",
    "conf.low",
    "conf.high",
    "family",
    "slope_sd"
  )
  tab <- do.call(
    rbind,
    lapply(rows, function(r) {
      r[, intersect(keep, names(r)), drop = FALSE]
    })
  )
  rownames(tab) <- NULL
  if ("p.value" %in% names(tab)) {
    tab$p.fdr <- pollin_fdr(tab$p.value)
  }
  structure(
    list(table = tab, fits = fits, spec = spec),
    class = "pollin_contribution"
  )
}

# Internal checks -------------------------------------------------------------

#' Stop unless the estimand is dependence on the difference scale
#' @noRd
#' @keywords internal
.check_contribution_estimand <- function(estimand) {
  if (estimand$family != "dependence" || estimand$scale != "difference") {
    stop(
      "`pollin_contribution()` reports dependence on the difference scale; ",
      "this estimand is ",
      estimand$family,
      " on the ",
      estimand$scale,
      " scale. Use `pollin_effects(pollin_fit(estimand))`.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Stop unless every outcome is a declared response
#' @noRd
#' @keywords internal
.check_outcomes <- function(outcomes, responses) {
  if (!is.character(outcomes) || length(outcomes) == 0L) {
    stop("`outcomes` must be a character vector of responses.", call. = FALSE)
  }
  unknown <- setdiff(outcomes, responses)
  if (length(unknown)) {
    stop(
      "Outcomes not declared as responses: ",
      paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}
