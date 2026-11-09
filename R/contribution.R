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
#' @param x A [pollin_estimand()] object with `family = "dependence"` and
#'   `scale = "difference"`.
#' @param outcomes Outcomes to model. `NULL` (the default) models every
#'   response.
#' @param by Optional name of a column for conditional contrasts, one per
#'   level (for example a region).
#' @param families Optional named list of `glmmTMB` family objects, one per
#'   outcome; outcomes it does not name are chosen from their values.
#' @inheritParams pollin_fit
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
#' d <- data.frame(
#'   trt = rep(rep(c("Cage", "Open"), each = 10L), 8L),
#'   Site = rep(paste0("s", seq_len(8L)), each = 20L)
#' )
#' d$pods <- rnbinom(nrow(d), mu = exp(3 + 0.3 * (d$trt == "Open")), size = 10)
#' pd <- pollin_data(d, "pods", "trt", unit = "Site")
#' ds <- pollin_design(pd, c(excluded = "Cage", open = "Open"), "rcbd")
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
