# Target population (layer L2) ------------------------------------------------

#' Declare a population of strata for an estimand
#'
#' A population made of the levels (strata) of one column of the data, for
#' example distance bands, regions or varieties, each weighted by its area.
#' Passed to [pollin_estimand()] as `population`, it carries the estimand
#' from the trial to the population: within each stratum the expected
#' response under each arm comes from the fitted model, and the population
#' arm means are the area-weighted sums of the stratum arm means:
#'
#' `high_pop = sum(w * high_stratum)` and `low_pop = sum(w * low_stratum)`,
#'
#' with the weights `w` normalised to sum to 1. The estimand is formed from
#' the two population means, so on the ratio scale the population loss share
#' is `1 - sum(w * low_stratum) / sum(w * high_stratum)`: a ratio of
#' population means, not a weighted mean of the stratum loss shares. The two
#' agree when every stratum has the same expected response in the `high` arm,
#' or the same loss share; in general they differ.
#'
#' The model must be fitted with the strata as its `by` column,
#' `pollin_fit(estimand, by = <by>)`, and the population estimand is read with
#' [pollin_draws()].
#'
#' @param by Name of the column whose levels are the strata.
#' @param weights Named numeric vector of stratum areas, one positive value per
#'   level of `by`, named by level. Any unit of area serves: the weights are
#'   normalised to sum to 1 and the values given are kept as `area`.
#' @param label Optional short description of the population, for printing.
#'
#' @return An object of class `pollin_population`: a list holding `by`, the
#'   normalised `weights`, the `area` values as given, and `label`.
#' @seealso [pollin_estimand()], [pollin_draws()].
#' @export
#' @examples
#' pollin_population(
#'   "Region",
#'   c(north = 1200, south = 800),
#'   label = "two growing regions"
#' )
pollin_population <- function(by, weights, label = NULL) {
  if (!is.character(by) || length(by) != 1L || is.na(by) || !nzchar(by)) {
    stop("`by` must be one column name.", call. = FALSE)
  }
  .population_check_weights(weights)
  if (
    !is.null(label) &&
      (!is.character(label) || length(label) != 1L || is.na(label))
  ) {
    stop("`label` must be NULL or one string.", call. = FALSE)
  }
  structure(
    list(
      by = by,
      weights = weights / sum(weights),
      area = weights,
      label = label
    ),
    class = "pollin_population"
  )
}

#' @export
print.pollin_population <- function(x, ...) {
  cat("<pollin_population>\n")
  if (!is.null(x$label)) {
    cat("  label   :", x$label, "\n")
  }
  cat("  strata  :", .n_noun(length(x$weights), "level"), "of", x$by, "\n")
  cat(
    "  weights :",
    paste0(
      names(x$weights),
      " ",
      format(x$weights, digits = 3L),
      collapse = ", "
    ),
    "(area, normalised)\n"
  )
  invisible(x)
}

# Internal ---------------------------------------------------------------------

#' Population arm means from stratum arm-mean draws: the area-weighted sums
#' of the stratum columns, appended as a column `"population"`
#' @noRd
#' @keywords internal
.population_means <- function(high, low, population) {
  .check_class(population, "pollin_population", arg = "population")
  w <- population$weights
  absent <- setdiff(names(w), colnames(high))
  if (length(absent)) {
    stop(
      "The draws lack the strata ",
      paste(absent, collapse = ", "),
      " of the population; the model data must hold every stratum.",
      call. = FALSE
    )
  }
  if ("population" %in% colnames(high)) {
    stop(
      "A column of the draws is already named \"population\".",
      call. = FALSE
    )
  }
  add <- function(m) {
    out <- cbind(m, m[, names(w), drop = FALSE] %*% w)
    colnames(out)[ncol(out)] <- "population"
    out
  }
  list(high = add(high), low = add(low))
}

#' Stop unless the population's strata are exactly the levels of its `by`
#' column in the data
#' @noRd
#' @keywords internal
.population_check_data <- function(population, data) {
  by <- population$by
  if (!by %in% names(data)) {
    stop(
      "The population's `by` column '",
      by,
      "' is not in the data.",
      call. = FALSE
    )
  }
  lev <- unique(as.character(data[[by]]))
  lev <- lev[!is.na(lev)]
  if (!setequal(names(population$weights), lev)) {
    stop(
      "The population's weights must be named by the levels of '",
      by,
      "': ",
      paste(sort(lev), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Stop unless the fit was made with the population's `by` column
#' @noRd
#' @keywords internal
.population_check_fit_by <- function(population, fit_by) {
  if (!identical(fit_by, population$by)) {
    stop(
      "The estimand's population has strata in '",
      population$by,
      "'; fit with `pollin_fit(estimand, by = \"",
      population$by,
      "\")`.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Stop unless `weights` is a named vector of positive finite numbers, one per
#' distinct name
#' @noRd
#' @keywords internal
.population_check_weights <- function(weights) {
  nm <- names(weights)
  if (
    !is.numeric(weights) ||
      length(weights) == 0L ||
      is.null(nm) ||
      anyNA(nm) ||
      !all(nzchar(nm)) ||
      anyDuplicated(nm)
  ) {
    stop(
      "`weights` must be a numeric vector named by stratum, each name once.",
      call. = FALSE
    )
  }
  if (any(!is.finite(weights)) || any(weights <= 0)) {
    stop("`weights` must be positive finite areas.", call. = FALSE)
  }
  invisible(TRUE)
}

#' The population assumption for a `pollin_population`
#' @noRd
#' @keywords internal
.population_statement <- function(population) {
  w <- population$weights
  paste0(
    "The population is the strata of `",
    population$by,
    "` (",
    paste0(names(w), " ", format(w, digits = 3L), collapse = ", "),
    "), weighted by area. Within each stratum the expected response under ",
    "each arm comes from the fitted model and is taken to hold across the ",
    "stratum's area; the population estimand is formed from the area-weighted ",
    "arm means."
  )
}
