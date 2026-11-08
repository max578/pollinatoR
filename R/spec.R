# Study specification ---------------------------------------------------------

#' Declare a pollinator-exclusion study
#'
#' Captures the structure of an open-versus-exclusion experiment in one object
#' so the downstream verbs share a single, validated description of the data.
#' The specification is deliberately general: outcomes, the exclusion treatment,
#' the nested random structure, optional exposures (the pollinator community),
#' optional adjustment covariates, an optional fixed grouping, and an optional
#' latent unit on which a random treatment slope is carried to absorb
#' unobserved between-unit heterogeneity (for example an unrecorded genotype).
#'
#' `pollin_spec()` is deprecated: [pollin_data()] and [pollin_design()]
#' declare the same study in two layers. It warns once per session and is
#' removed at version 0.2.0.
#'
#' @param data A data frame, one row per observational unit (e.g. one plant).
#' @param outcomes Character vector of outcome column names (yield components).
#' @param treatment Name of the treatment column (a factor with two or more
#'   levels, e.g. caged versus open).
#' @param control Treatment level that represents exclusion / the reference
#'   (e.g. `"Cage"`). The contribution is computed relative to this level.
#' @param random Character vector of grouping columns giving the nesting, outer
#'   to inner (e.g. `c("Site", "Distance", "Transect")`). Translated into nested
#'   random intercepts.
#' @param exposures Optional character vector of exposure columns (e.g. insect
#'   abundances) used by [pollin_moderation()].
#' @param adjust Optional character vector of adjustment covariates (e.g.
#'   landscape, climate) -- context, not the inferential target.
#' @param group Optional name of a fixed grouping column (e.g. region).
#' @param latent_unit Optional name of the unit carrying the random treatment
#'   slope (defaults to the outermost `random`). This absorbs unobserved
#'   unit-level modifiers of the contribution, such as genotype.
#'
#' @return An object of class `pollin_spec`.
#' @export
#' @examples
#' df <- data.frame(
#'   y = rpois(40, 5), trt = rep(c("Cage", "Open"), 20),
#'   Site = rep(letters[1:4], each = 10)
#' )
#' pollin_spec(df, outcomes = "y", treatment = "trt", control = "Cage",
#'             random = "Site")
pollin_spec <- function(
  data,
  outcomes,
  treatment,
  control,
  random,
  exposures = NULL,
  adjust = NULL,
  group = NULL,
  latent_unit = NULL
) {
  .warn_spec_deprecated()
  stopifnot(is.data.frame(data))
  cols <- names(data)
  need <- c(outcomes, treatment, random, exposures, adjust, group)
  miss <- setdiff(need, cols)
  if (length(miss)) {
    stop(
      "Columns not found in `data`: ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }
  trt_levels <- unique(as.character(data[[treatment]]))
  if (!control %in% trt_levels) {
    stop(
      "`control` level '",
      control,
      "' not present in `",
      treatment,
      "`.",
      call. = FALSE
    )
  }
  if (is.null(latent_unit)) {
    latent_unit <- random[1]
  }
  if (!latent_unit %in% random) {
    stop("`latent_unit` must be one of `random`.", call. = FALSE)
  }
  structure(
    list(
      data = data,
      outcomes = outcomes,
      treatment = treatment,
      control = control,
      treatment_levels = trt_levels,
      random = random,
      exposures = exposures,
      adjust = adjust,
      group = group,
      latent_unit = latent_unit,
      n = nrow(data)
    ),
    class = "pollin_spec"
  )
}

# Deprecation ------------------------------------------------------------------

#' Session state of the package (the once-per-session deprecation flag)
#' @noRd
#' @keywords internal
.pollin_env <- new.env(parent = emptyenv())

#' Warn, once per session, that `pollin_spec()` is deprecated
#' @noRd
#' @keywords internal
.warn_spec_deprecated <- function() {
  if (isTRUE(.pollin_env$spec_warned)) {
    return(invisible(FALSE))
  }
  .pollin_env$spec_warned <- TRUE
  warning(
    "`pollin_spec()` is deprecated; use `pollin_data()` + `pollin_design()`. ",
    "It is removed at pollinatoR 0.2.0. This warning is shown once per session.",
    call. = FALSE
  )
  invisible(TRUE)
}

# Build the random-effects term string for a glmmTMB formula: nested
# intercepts with a treatment slope on the latent unit (see `.nested_terms()`).
.random_terms <- function(spec) {
  .nested_terms(spec$random, spec$treatment, spec$latent_unit)
}

#' @export
print.pollin_spec <- function(x, ...) {
  cat("<pollin_spec>\n")
  cat("  observations :", x$n, "\n")
  cat("  outcomes     :", paste(x$outcomes, collapse = ", "), "\n")
  cat(
    "  treatment    :",
    x$treatment,
    "(control =",
    x$control,
    ";",
    length(x$treatment_levels),
    "levels)\n"
  )
  cat(
    "  random       :",
    paste(x$random, collapse = " / "),
    " [slope on:",
    x$latent_unit,
    "]\n"
  )
  if (!is.null(x$exposures)) {
    cat("  exposures    :", paste(x$exposures, collapse = ", "), "\n")
  }
  if (!is.null(x$adjust)) {
    cat("  adjust       :", paste(x$adjust, collapse = ", "), "\n")
  }
  if (!is.null(x$group)) {
    cat("  group        :", x$group, "\n")
  }
  invisible(x)
}
