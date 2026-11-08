# E3: yield-component compensation --------------------------------------------

#' Summarise yield-component trade-offs (compensation)
#'
#' Canola yield is the product of compensating components (pod number, seeds per
#' pod, seed size): a deficit in one is often offset by another, which a
#' component-by-component analysis can hide. This function returns the
#' correlation structure among the components across units and flags negative
#' associations as candidate compensation. It is a descriptive companion to the
#' per-component contribution models, not a substitute for a joint multivariate
#' model.
#'
#' @param components A numeric matrix or data frame, one column per yield
#'   component, one row per unit; or a [pollin_data()] object, whose response
#'   columns are taken as the components (at least two).
#' @param method Correlation method passed to [stats::cor()].
#'
#' @return An object of class `pollin_compensation` with the correlation matrix
#'   and a table of negatively-correlated component pairs.
#' @export
#' @examples
#' m <- data.frame(pods = c(10, 8, 6), seeds_per_pod = c(20, 24, 28))
#' pollin_compensation(m)
#'
#' # From a pollin_data object: its responses are the components
#' df <- data.frame(
#'   pods = c(10, 8, 6, 9, 7, 5),
#'   seeds_per_pod = c(20, 24, 28, 21, 25, 27),
#'   trt = rep(c("Cage", "Open"), 3),
#'   Site = rep(c("a", "b", "c"), each = 2)
#' )
#' pd <- pollin_data(df, c("pods", "seeds_per_pod"), "trt", unit = "Site")
#' pollin_compensation(pd)
pollin_compensation <- function(components, method = "spearman") {
  if (inherits(components, "pollin_data")) {
    .check_n_responses(components$roles$response)
    components <- components$data[components$roles$response]
  }
  components <- as.data.frame(components)
  num <- vapply(components, is.numeric, logical(1))
  components <- components[, num, drop = FALSE]
  if (ncol(components) < 2L) {
    stop("Need at least two numeric components.", call. = FALSE)
  }
  cm <- stats::cor(components, use = "pairwise.complete.obs", method = method)
  pairs <- which(lower.tri(cm), arr.ind = TRUE)
  tab <- data.frame(
    a = colnames(cm)[pairs[, 1L]],
    b = colnames(cm)[pairs[, 2L]],
    cor = cm[lower.tri(cm)],
    stringsAsFactors = FALSE
  )
  tab <- tab[order(tab$cor), , drop = FALSE]
  tab$compensation <- tab$cor < 0
  structure(
    list(cor = cm, pairs = tab, method = method),
    class = "pollin_compensation"
  )
}

#' @export
print.pollin_compensation <- function(x, ...) {
  cat("<pollin_compensation>  (", x$method, "correlations )\n")
  print(x$pairs, row.names = FALSE)
  neg <- sum(x$pairs$compensation)
  cat(
    "  ",
    .n_noun(neg, "negatively-correlated pair"),
    " -> candidate compensation\n",
    sep = ""
  )
  invisible(x)
}

# Internal checks -------------------------------------------------------------

#' Stop unless at least two responses are declared
#' @noRd
#' @keywords internal
.check_n_responses <- function(response) {
  if (length(response) < 2L) {
    stop(
      "Compensation needs at least two responses in the `pollin_data` ",
      "object; found ",
      length(response),
      ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
