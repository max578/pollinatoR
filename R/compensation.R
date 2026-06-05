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
#'   component, one row per unit.
#' @param method Correlation method passed to [stats::cor()].
#'
#' @return An object of class `pollin_compensation` with the correlation matrix
#'   and a table of negatively-correlated component pairs.
#' @export
#' @examples
#' m <- data.frame(pods = c(10, 8, 6), seeds_per_pod = c(20, 24, 28))
#' pollin_compensation(m)
pollin_compensation <- function(components, method = "spearman") {
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
  cat(sprintf(
    "  %d negatively-correlated pair(s) -> candidate compensation\n",
    neg
  ))
  invisible(x)
}
