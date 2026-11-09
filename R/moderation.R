# E2: which exposures drive the contribution (within-between robust) ----------

# Standardise with NA handling.
.z <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) {
    return(x - m)
  }
  (x - m) / s
}

#' Estimate exposure moderation of the contribution
#'
#' Asks which exposures (e.g. honey-bee versus other-insect abundance) drive
#' the open-minus-exclusion contribution, as treatment-by-exposure
#' interactions. Each exposure enters through a within-between (Mundlak)
#' decomposition: a between-unit part (the unit mean) and a within-unit part
#' (the deviation from it). The **within-unit** treatment-by-exposure
#' coefficient is robust to any unit-level confounder -- observed or not,
#' including an unrecorded genotype -- because it is identified from variation
#' that occurs within the unit, where
#' the confounder is constant. The between-unit part is reported but flagged as
#' confounder-vulnerable. Adjustment covariates enter as fixed effects.
#'
#' From a [pollin_design()], the unit of the split is the design's slope unit
#' (with the unit levels outside it), or the outermost unit when the design
#' carries no slope; the treatment reference level is the excluded arm, or the
#' open arm when the design has no excluded arm.
#'
#' @param x A [pollin_design()] object whose data declare `exposures`.
#' @param outcome Single outcome to model.
#' @param exposures Exposures to use. `NULL` (the default) uses every declared
#'   exposure.
#' @param family Optional `glmmTMB` family (chosen from the outcome's values if
#'   missing).
#' @inheritParams pollin_fit
#'
#' @return An object of class `pollin_moderation` with a tidy `table` of the
#'   treatment-by-exposure terms (each tagged `within` or `between`),
#'   FDR-adjusted across the within terms, and the fitted model in `fit`.
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' site <- rep(paste0("s", seq_len(14L)), each = 20L)
#' trt <- rep(rep(c("Cage", "Open"), each = 10L), 14L)
#' bees <- rep(rpois(14L, 8), each = 20L) + rpois(280L, 1)
#' y <- rpois(length(trt), exp(1.6 + 0.3 * (trt == "Open") +
#'                             0.02 * bees * (trt == "Open")))
#' d <- data.frame(y = y, trt = trt, Site = site, bees = bees)
#' pd <- pollin_data(d, "y", "trt", unit = "Site", exposures = "bees")
#' ds <- pollin_design(pd, c(excluded = "Cage", open = "Open"), "multi_site")
#' pollin_moderation(ds, outcome = "y")
#' }
pollin_moderation <- function(
  x,
  outcome,
  exposures = NULL,
  family = NULL,
  ...
) {
  .check_class(x, "pollin_design")
  .check_outcomes(outcome, x$data$roles$response)
  low <- x$arms[[intersect(c("excluded", "open"), names(x$arms))[1L]]]
  inputs <- .design_inputs(x, low)
  if (is.null(exposures)) {
    exposures <- inputs$exposures
  }
  if (is.null(exposures) || length(exposures) == 0L) {
    stop("No exposures supplied (declare them in the data).", call. = FALSE)
  }

  # Within-between split of each exposure --------------------------------------
  d <- inputs$data
  ex_terms <- character(0)
  for (ex in exposures) {
    wb <- pollin_within_between(d[[ex]], inputs$wb_group)
    wn <- paste0(ex, "_within")
    bn <- paste0(ex, "_between")
    d[[wn]] <- .z(wb$within)
    d[[bn]] <- .z(wb$between)
    ex_terms <- c(ex_terms, wn, bn)
  }

  # Fit and tidy the interaction terms -----------------------------------------
  fam <- if (!is.null(family)) family else .pick_family(d[[outcome]])
  fixed <- sprintf(
    "%s * (%s)",
    inputs$treatment,
    paste(ex_terms, collapse = " + ")
  )
  rhs <- c(fixed, inputs$adjust, inputs$random)
  fit <- .fit_glmm(outcome, rhs, d, fam, ...)
  co <- summary(fit)$coefficients$cond
  is_int <- grepl(":", rownames(co))
  tab <- data.frame(
    term = rownames(co)[is_int],
    estimate = co[is_int, 1L],
    std.error = co[is_int, 2L],
    p.value = co[is_int, 4L],
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  tab$component <- ifelse(
    grepl("_within", tab$term),
    "within (robust)",
    ifelse(grepl("_between", tab$term), "between (vulnerable)", "other")
  )
  # FDR across the robust (within) interaction terms only
  tab$p.fdr <- NA_real_
  w <- tab$component == "within (robust)"
  if (any(w)) {
    tab$p.fdr[w] <- pollin_fdr(tab$p.value[w])
  }
  structure(
    list(table = tab, fit = fit, outcome = outcome, design = x),
    class = "pollin_moderation"
  )
}

#' @export
print.pollin_moderation <- function(x, ...) {
  cat("<pollin_moderation>  outcome:", x$outcome, " (treatment x exposure)\n")
  print(x$table, row.names = FALSE)
  invisible(x)
}
