# E1: the pollinator contribution (open minus exclusion) ----------------------

#' Estimate the pollinator contribution per outcome
#'
#' Fits, for each outcome, a generalised linear mixed model with the exclusion
#' treatment as the fixed effect of interest, the nested random structure from
#' the specification, and a random treatment slope on the latent unit so that
#' unobserved between-unit heterogeneity in the contribution (for example an
#' unrecorded genotype) is absorbed rather than ignored. The contribution is
#' reported as a model-based average treatment effect (the open-minus-exclusion
#' contrast) with propagated uncertainty, via G-computation. Family is chosen
#' per outcome by a simple heuristic unless supplied.
#'
#' The treatment effect is causal by design here -- the exclusion contrast is an
#' experimental manipulation -- so no covariate adjustment is applied; see
#' [pollin_moderation()] for the observational "which insects" question.
#'
#' @param spec A [pollin_spec()] object.
#' @param outcomes Outcomes to model (defaults to all in `spec`).
#' @param by Optional grouping for conditional contrasts (e.g. `spec$group`).
#' @param families Optional named list of `glmmTMB` family objects, one per
#'   outcome; missing entries are chosen heuristically.
#' @param ... Passed to [glmmTMB::glmmTMB()].
#'
#' @return An object of class `pollin_contribution`: a list with `table` (a tidy
#'   data frame of contrasts with FDR-adjusted p-values and the latent-unit
#'   treatment-slope SD) and `fits` (the fitted models, for diagnostics).
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' site <- rep(paste0("s", seq_len(12L)), each = 20L)
#' trt <- rep(rep(c("Cage", "Open"), each = 10L), 12L)
#' y <- rpois(length(trt), exp(1.7 + 0.35 * (trt == "Open")))
#' d <- data.frame(y = y, trt = trt, Site = site)
#' sp <- pollin_spec(d, "y", "trt", "Cage", "Site")
#' pollin_contribution(sp)
#' }
pollin_contribution <- function(
  spec,
  outcomes = spec$outcomes,
  by = NULL,
  families = list(),
  ...
) {
  stopifnot(inherits(spec, "pollin_spec"))
  rows <- list()
  fits <- list()
  for (oc in outcomes) {
    d <- spec$data
    d[[spec$treatment]] <- stats::relevel(
      factor(d[[spec$treatment]]),
      ref = spec$control
    )
    fam <- if (!is.null(families[[oc]])) {
      families[[oc]]
    } else {
      .pick_family(d[[oc]])
    }
    form <- stats::reformulate(
      c(spec$treatment, .random_terms(spec)),
      response = oc
    )
    fit <- tryCatch(
      glmmTMB::glmmTMB(form, data = d, family = fam, ...),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      warning(
        "Outcome '",
        oc,
        "' failed to fit: ",
        conditionMessage(fit),
        call. = FALSE
      )
      next
    }
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
    vc <- tryCatch(
      glmmTMB::VarCorr(fit)$cond[[spec$latent_unit]],
      error = function(e) NULL
    )
    slope_sd <- if (!is.null(vc) && nrow(vc) >= 2L) {
      sqrt(diag(vc))[2]
    } else {
      NA_real_
    }
    cmp$outcome <- oc
    cmp$family <- fam$family
    cmp$slope_sd <- unname(slope_sd)
    rows[[oc]] <- as.data.frame(cmp)
    fits[[oc]] <- fit
  }
  if (length(rows) == 0L) {
    stop("No outcome could be fitted.", call. = FALSE)
  }
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

#' @export
print.pollin_contribution <- function(x, ...) {
  cat("<pollin_contribution>  (open - exclusion, model-based ATE)\n")
  print(x$table, row.names = FALSE)
  invisible(x)
}
