# Structural data-integrity checks --------------------------------------------

#' Validate the structure and internal consistency of a study's data
#'
#' Runs aggregate structural checks against a declared table: missingness per
#' declared column, expected factor levels, and any user-supplied
#' internal-consistency identities (for example that a total equals a sum of
#' parts, or that a ratio equals its definition). Returns counts and pass/warn
#' flags only -- no raw rows -- so it is safe on sensitive data.
#' [pollin_data()] runs these checks without levels or identities when it
#' builds its object.
#'
#' @param x A [pollin_data()] object, or a [pollin_spec()] object.
#' @param expect_levels Optional named list mapping column names to the expected
#'   set of factor levels.
#' @param identities Optional named list of functions; each takes the data frame
#'   and returns a logical vector that is `TRUE` where the identity holds. The
#'   number of violations is reported.
#'
#' @return An object of class `pollin_validate` carrying a tidy `table`.
#' @export
#' @examples
#' df <- data.frame(y = c(1, 2, NA), trt = c("Cage", "Open", "Open"),
#'                  Site = c("a", "a", "b"))
#' pd <- pollin_data(df, response = "y", treatment = "trt", unit = "Site")
#' pollin_validate(pd, expect_levels = list(trt = c("Cage", "Open")))
pollin_validate <- function(x, expect_levels = NULL, identities = NULL) {
  declared <- .declared_columns(x)
  d <- x$data
  rows <- list()
  add <- function(check, status, detail = "") {
    rows[[length(rows) + 1L]] <<- data.frame(
      check = check,
      status = status,
      detail = as.character(detail),
      stringsAsFactors = FALSE
    )
  }
  add("rows", "INFO", nrow(d))
  add("cols", "INFO", ncol(d))
  for (cc in declared) {
    na <- sum(is.na(d[[cc]]))
    if (na > 0L) add(paste0("NA:", cc), "WARN", na)
  }
  if (!is.null(expect_levels)) {
    for (cc in names(expect_levels)) {
      got <- sort(unique(as.character(d[[cc]])))
      ok <- setequal(got, expect_levels[[cc]])
      add(
        paste0("levels:", cc),
        if (ok) "PASS" else "WARN",
        paste(got, collapse = ",")
      )
    }
  }
  if (!is.null(identities)) {
    for (nm in names(identities)) {
      v <- identities[[nm]](d)
      viol <- sum(!v, na.rm = TRUE)
      add(
        paste0("identity:", nm),
        if (viol == 0L) "PASS" else "WARN",
        paste0("violations=", viol)
      )
    }
  }
  structure(list(table = do.call(rbind, rows)), class = "pollin_validate")
}

#' @export
print.pollin_validate <- function(x, ...) {
  cat("<pollin_validate>\n")
  print(x$table, row.names = FALSE)
  invisible(x)
}

#' Declared columns of a `pollin_data` or `pollin_spec`, in role order
#' @noRd
#' @keywords internal
.declared_columns <- function(x) {
  if (inherits(x, "pollin_data")) {
    return(unique(unlist(x$roles, use.names = FALSE)))
  }
  if (inherits(x, "pollin_spec")) {
    return(unique(c(
      x$outcomes,
      x$treatment,
      x$random,
      x$exposures,
      x$adjust,
      x$group
    )))
  }
  stop("`x` must be a `pollin_data` or `pollin_spec` object.", call. = FALSE)
}
