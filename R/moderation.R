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
#' (the deviation from it), each standardised. The **within-unit**
#' treatment-by-exposure coefficient is robust to any unit-level confounder --
#' observed or not, including an unrecorded genotype -- because it is
#' identified from variation that occurs within the unit, where the confounder
#' is constant. The between-unit part is reported but flagged as
#' confounder-vulnerable. An exposure that does not vary within units (a
#' site-level count) has no within part: it is fitted between units only, and
#' a message says so.
#'
#' From a [pollin_design()], the unit of the split is the design's slope unit
#' (with the unit levels outside it), or the outermost unit when the design
#' carries no slope; the treatment reference level is the excluded arm, or the
#' open arm when the design has no excluded arm.
#'
#' **Adjustment.** Without a graph, the declared `adjust` columns enter as
#' main effects. With `dag`, the adjustment set of each exposure (its first
#' minimal set for the exposure's effect on the outcome, from [pollin_dag()])
#' enters as main effects and treatment interactions, since a confounder that
#' changes the contrast is removed only by its interaction. A graph cannot
#' draw effect modification: the set is borrowed from the exposure's effect on
#' the outcome, and the assumptions say so.
#'
#' **Negative control.** Insects cannot reach the excluded plants, so an
#' exposure's slope in the excluded arm is not pollination: it is confounding,
#' or a site effect correlated with the cage. With the excluded arm as the
#' reference level that slope is the exposure's main effect, reported in
#' `negative_control`. It detects a site variable that shifts both arms alike
#' on the link scale; it does not detect one that changes the contrast itself.
#'
#' A fit without a positive-definite Hessian, or that did not converge, has
#' unreliable standard errors and p-values; it returns with a warning naming
#' the outcome.
#'
#' @param x A [pollin_design()] object whose data declare `exposures`.
#' @param outcome Single outcome to model.
#' @param exposures Exposures to use. `NULL` (the default) uses every declared
#'   exposure.
#' @param family Optional `glmmTMB` family (chosen from the outcome's values if
#'   missing).
#' @param dag Optional [pollin_dag()] object built from `x`; its adjustment
#'   sets replace the declared `adjust` columns.
#' @param ... Passed to [glmmTMB::glmmTMB()].
#' @inheritParams pollin_fit
#'
#' @return An object of class `pollin_moderation`: a list holding the tidy
#'   `table` of the treatment-by-exposure terms (each tagged `within`,
#'   `between` or, for a graph's adjustment columns, `adjustment`),
#'   FDR-adjusted across the within terms; the `negative_control` table (the
#'   exposure terms' slopes in the excluded arm, on the link scale per
#'   standard deviation; `NULL` without an excluded arm); `between_only` (the
#'   exposures with no within-unit variation); `adjustment` (the columns
#'   adjusted for); the `assumptions` (a data frame of `assumption` and
#'   `statement`); the fitted model in `fit`; the `outcome` and the `design`.
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
  dag = NULL,
  ...
) {
  .check_class(x, "pollin_design")
  .check_outcomes(outcome, x$data$roles$response)
  ref <- intersect(c("excluded", "open"), names(x$arms))[1L]
  low <- x$arms[[ref]]
  inputs <- .design_inputs(x, low)
  if (is.null(exposures)) {
    exposures <- inputs$exposures
  }
  if (is.null(exposures) || length(exposures) == 0L) {
    stop("No exposures supplied (declare them in the data).", call. = FALSE)
  }
  adjust <- .moderation_adjustment(x, dag, exposures)

  # Within-between split of each exposure --------------------------------------
  d <- inputs$data
  ex_terms <- character(0)
  between_only <- character(0)
  for (ex in exposures) {
    wb <- pollin_within_between(d[[ex]], inputs$wb_group)
    wn <- paste0(ex, "_within")
    bn <- paste0(ex, "_between")
    tol <- 1e-8 * max(1, abs(d[[ex]]), na.rm = TRUE)
    if (all(abs(wb$within) < tol, na.rm = TRUE)) {
      message(
        "`",
        ex,
        "` does not vary within units; it is fitted between units only."
      )
      between_only <- c(between_only, ex)
    } else {
      d[[wn]] <- .z(wb$within)
      ex_terms <- c(ex_terms, wn)
    }
    d[[bn]] <- .z(wb$between)
    ex_terms <- c(ex_terms, bn)
  }

  # Fit and tidy the interaction terms -----------------------------------------
  fam <- if (!is.null(family)) family else .pick_family(d[[outcome]])
  fixed <- sprintf(
    "%s * (%s)",
    inputs$treatment,
    paste(c(ex_terms, if (!is.null(dag)) adjust), collapse = " + ")
  )
  rhs <- c(fixed, if (is.null(dag)) adjust, inputs$random)
  fit <- .fit_glmm(outcome, rhs, d, fam, ...)
  if (!isTRUE(fit$sdr$pdHess) || fit$fit$convergence != 0L) {
    warning(
      "The moderation fit for `",
      outcome,
      "` ",
      if (fit$fit$convergence != 0L) {
        "did not converge"
      } else {
        "has no positive-definite Hessian"
      },
      "; its standard errors and p-values are not reliable.",
      call. = FALSE
    )
  }
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
    ifelse(grepl("_between", tab$term), "between (vulnerable)", "adjustment")
  )
  # FDR across the robust (within) interaction terms only
  tab$p.fdr <- NA_real_
  w <- tab$component == "within (robust)"
  if (any(w)) {
    tab$p.fdr[w] <- pollin_fdr(tab$p.value[w])
  }

  # Negative control: the exposure slopes in the excluded arm ------------------
  negative <- if (ref == "excluded") {
    data.frame(
      term = ex_terms,
      estimate = co[ex_terms, 1L],
      std.error = co[ex_terms, 2L],
      p.value = co[ex_terms, 4L],
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  }

  structure(
    list(
      table = tab,
      negative_control = negative,
      between_only = between_only,
      adjustment = adjust,
      assumptions = .moderation_assumptions(adjust, dag, ref, between_only),
      fit = fit,
      outcome = outcome,
      design = x
    ),
    class = "pollin_moderation"
  )
}

#' @export
print.pollin_moderation <- function(x, ...) {
  cat("<pollin_moderation>  outcome:", x$outcome, " (treatment x exposure)\n")
  print(x$table, row.names = FALSE)
  if (!is.null(x$negative_control)) {
    cat("  negative control: exposure slopes in the excluded arm\n")
    print(x$negative_control, row.names = FALSE)
  }
  if (length(x$between_only)) {
    cat(
      "  between units only:",
      paste(x$between_only, collapse = ", "),
      "(no within-unit variation)\n"
    )
  }
  cat(
    "  adjusted for    :",
    if (length(x$adjustment)) paste(x$adjustment, collapse = ", ") else "none",
    "\n"
  )
  invisible(x)
}

#' The columns moderation adjusts for: the declared `adjust` columns, or with
#' a graph the union of each exposure's first minimal adjustment set; stops
#' when the graph identifies no set for an exposure
#' @noRd
#' @keywords internal
.moderation_adjustment <- function(x, dag, exposures) {
  roles <- x$data$roles
  if (is.null(dag)) {
    return(roles$adjust)
  }
  .check_class(dag, "pollin_dag", arg = "dag")
  if (!identical(dag$roles, roles)) {
    stop(
      "`dag` was built from a different design; rebuild it with ",
      "`pollin_dag(design, ...)`.",
      call. = FALSE
    )
  }
  sets <- dag$exposure_sets[exposures]
  none <- exposures[lengths(sets) == 0L]
  if (length(none)) {
    stop(
      "The graph identifies no adjustment set for the effect of ",
      paste0("`", none, "`", collapse = ", "),
      " on the outcome; its moderation term has no causal reading.",
      call. = FALSE
    )
  }
  setdiff(unique(unlist(lapply(sets, `[[`, 1L))), exposures)
}

#' The statements a moderation result rests on
#' @noRd
#' @keywords internal
.moderation_assumptions <- function(adjust, dag, ref, between_only) {
  set <- if (length(adjust)) .dag_set_label(adjust) else "{}"
  rows <- c(
    `between-unit moderation` = paste0(
      "No unmeasured unit-level variable associated with the exposure ",
      "changes the treatment contrast. A variable that shifts both arms ",
      "alike on the link scale cancels in the contrast and shows as a slope ",
      "in the excluded arm (the negative control)."
    ),
    `within-unit moderation` = if (length(between_only)) {
      paste0(
        "Not estimated for ",
        paste(between_only, collapse = ", "),
        ": no variation within units."
      )
    } else {
      paste0(
        "Robust to variables constant within a unit; not to a variable that ",
        "varies within units with the exposure."
      )
    },
    adjustment = if (is.null(dag)) {
      paste0(
        "The declared `adjust` columns ",
        set,
        ", as main effects; no graph was given."
      )
    } else {
      paste0(
        "The graph's adjustment set ",
        set,
        " for the exposures' effect on the outcome, as main effects and ",
        "treatment interactions. A graph cannot draw effect modification; ",
        "the set is borrowed from the exposures' effect and is valid for the ",
        "moderation term only if no unmeasured variable both changes the ",
        "treatment contrast and is associated with the exposure."
      )
    },
    `negative control` = if (ref == "excluded") {
      paste0(
        "Insects cannot reach the excluded plants, so an exposure's slope ",
        "there is confounding or a cage-correlated site effect. It detects ",
        "an equal shift of both arms, not a variable that changes the ",
        "contrast."
      )
    } else {
      "Not available: the design has no excluded arm."
    }
  )
  data.frame(
    assumption = names(rows),
    statement = unname(rows),
    stringsAsFactors = FALSE
  )
}
