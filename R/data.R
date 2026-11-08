# Unit-of-record data (layer L0) ----------------------------------------------

#' Declare the unit-of-record table of a pollination study
#'
#' Pairs a data frame with a role map naming which columns are responses, the
#' pollination treatment, the units the treatment is applied within, and the
#' optional grouping, denominator, exposure and adjustment columns. The roles
#' are checked against the data when the object is built, and the structural
#' checks of [pollin_validate()] are run and stored, so every later step works
#' from one validated description of the table.
#'
#' `pollin_data()` records what the columns are. How the units are nested and
#' which treatment levels form the arms of the experiment are declared in the
#' design step that follows.
#'
#' @param data A data frame, one row per unit of record (for example one plant
#'   under one pollination treatment).
#' @param response Character vector of response column names (for example pod
#'   number, seed weight). Each must be numeric.
#' @param treatment Name of the pollination treatment column (for example caged
#'   versus open). It must have at least two distinct non-missing values.
#' @param unit Character vector of the columns identifying the units the
#'   treatment is applied within, outer to inner (for example
#'   `c("Site", "Plot")`).
#' @param groups Optional character vector of further grouping columns that are
#'   not part of the unit nesting (for example region, variety or year).
#' @param denominator Optional name of a numeric column giving what the response
#'   is counted out of (ovules, flowers, plants). Values must be non-negative.
#' @param exposures Optional character vector of exposure columns (for example
#'   insect abundances).
#' @param adjust Optional character vector of adjustment covariates (for example
#'   landscape or climate).
#'
#' @return An object of class `pollin_data`: a list holding `data`, the role map
#'   `roles`, the number of rows `n`, and the structural check table `checks`.
#' @seealso [pollin_validate()] for the checks, with expected levels and
#'   internal-consistency identities.
#' @export
#' @examples
#' df <- data.frame(
#'   pods = rpois(40, 30),
#'   trt = rep(c("Cage", "Open"), 20),
#'   Site = rep(letters[1:4], each = 10),
#'   bees = rpois(40, 5)
#' )
#' pd <- pollin_data(df, response = "pods", treatment = "trt", unit = "Site",
#'                   exposures = "bees")
#' pd
pollin_data <- function(
  data,
  response,
  treatment,
  unit,
  groups = NULL,
  denominator = NULL,
  exposures = NULL,
  adjust = NULL
) {
  roles <- list(
    response = response,
    treatment = treatment,
    unit = unit,
    groups = groups,
    denominator = denominator,
    exposures = exposures,
    adjust = adjust
  )
  .check_data_frame(data)
  .check_roles(roles, names(data))
  .check_role_columns(data, roles)

  x <- structure(
    list(data = data, roles = roles, n = nrow(data), checks = NULL),
    class = "pollin_data"
  )
  x$checks <- pollin_validate(x)$table
  x
}

#' @export
print.pollin_data <- function(x, ...) {
  r <- x$roles
  cat("<pollin_data>\n")
  cat("  rows        :", x$n, "\n")
  cat("  response    :", paste(r$response, collapse = ", "), "\n")
  cat(
    "  treatment   :",
    r$treatment,
    "(",
    paste(.treatment_levels(x$data[[r$treatment]]), collapse = ", "),
    ")\n"
  )
  cat("  unit        :", paste(r$unit, collapse = " / "), "\n")
  for (nm in c("groups", "denominator", "exposures", "adjust")) {
    if (!is.null(r[[nm]])) {
      cat(sprintf("  %-12s:", nm), paste(r[[nm]], collapse = ", "), "\n")
    }
  }
  n_warn <- sum(x$checks$status == "WARN")
  cat("  checks      :", n_warn, "warning(s); see `$checks`\n")
  invisible(x)
}

# Internal checks -------------------------------------------------------------

#' Distinct non-missing treatment values, in order of appearance
#' @noRd
#' @keywords internal
.treatment_levels <- function(x) {
  unique(as.character(x[!is.na(x)]))
}

#' Stop unless `data` is a data frame
#' @noRd
#' @keywords internal
.check_data_frame <- function(data) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Check the shape of the role map: character names, one treatment, at most
#' one denominator, every column present, no column in two roles
#' @noRd
#' @keywords internal
.check_roles <- function(roles, cols) {
  for (nm in names(roles)) {
    r <- roles[[nm]]
    if (!is.null(r) && (!is.character(r) || length(r) == 0L || anyNA(r))) {
      stop(
        "`",
        nm,
        "` must be a character vector of column names.",
        call. = FALSE
      )
    }
  }
  if (length(roles$treatment) != 1L) {
    stop("`treatment` must name one column.", call. = FALSE)
  }
  if (length(roles$denominator) > 1L) {
    stop("`denominator` must name one column.", call. = FALSE)
  }
  declared <- unlist(roles, use.names = FALSE)
  miss <- setdiff(declared, cols)
  if (length(miss)) {
    stop(
      "Columns not found in `data`: ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }
  dup <- unique(declared[duplicated(declared)])
  if (length(dup)) {
    stop(
      "Each column can take one role only; repeated: ",
      paste(dup, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Check column contents against their roles: numeric responses, a treatment
#' with two or more values, a numeric non-negative denominator
#' @noRd
#' @keywords internal
.check_role_columns <- function(data, roles) {
  is_num <- vapply(data[roles$response], is.numeric, logical(1))
  if (!all(is_num)) {
    stop(
      "Response columns must be numeric: ",
      paste(roles$response[!is_num], collapse = ", "),
      call. = FALSE
    )
  }
  n_lev <- length(.treatment_levels(data[[roles$treatment]]))
  if (n_lev < 2L) {
    stop(
      "`",
      roles$treatment,
      "` must have at least two distinct values; found ",
      n_lev,
      ".",
      call. = FALSE
    )
  }
  if (!is.null(roles$denominator)) {
    den <- data[[roles$denominator]]
    if (!is.numeric(den) || any(den < 0, na.rm = TRUE)) {
      stop(
        "`denominator` column '",
        roles$denominator,
        "' must be numeric and non-negative.",
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}
