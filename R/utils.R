# Internal helpers ------------------------------------------------------------

#' Within-between (Mundlak) decomposition of a numeric variable
#'
#' Splits a numeric vector into a between-group part (the group mean, broadcast
#' back to each row) and a within-group part (the deviation from that mean).
#' Used to separate the within-unit component of an exposure -- which is robust
#' to any unit-level confounder, observed or not -- from the between-unit
#' component.
#'
#' @param x Numeric vector.
#' @param group Grouping factor of the same length as `x`.
#'
#' @return A list with elements `within` and `between`, each the length of `x`.
#' @export
#' @examples
#' wb <- pollin_within_between(c(1, 3, 10, 12), c("a", "a", "b", "b"))
#' wb$between
#' wb$within
pollin_within_between <- function(x, group) {
  if (length(x) != length(group)) {
    stop("`x` and `group` must have the same length.", call. = FALSE)
  }
  group <- as.character(group)
  gm <- tapply(x, group, function(z) mean(z, na.rm = TRUE))
  between <- as.numeric(gm[group])
  list(within = x - between, between = between)
}

#' Benjamini-Hochberg false-discovery-rate adjustment
#'
#' Thin wrapper over [stats::p.adjust()] with the BH method, kept as a named
#' verb so the multiplicity step is explicit in analysis code.
#'
#' @param p Numeric vector of p-values.
#'
#' @return Numeric vector of FDR-adjusted p-values.
#' @export
#' @examples
#' pollin_fdr(c(0.001, 0.02, 0.30, 0.80))
pollin_fdr <- function(p) {
  stats::p.adjust(p, method = "BH")
}

# Heuristic family choice for an outcome. Returns a glmmTMB family object.
.pick_family <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(glmmTMB::nbinom2())
  }
  is_count <- all(abs(x - round(x)) < 1e-8) && all(x >= 0)
  in_unit <- all(x > 0 & x < 1)
  has_zero <- any(x == 0)
  if (in_unit) {
    return(glmmTMB::beta_family())
  }
  if (is_count) {
    return(glmmTMB::nbinom2())
  }
  if (has_zero) {
    return(glmmTMB::tweedie())
  }
  stats::Gamma(link = "log")
}

# Require a suggested package or stop with a helpful message.
.need <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      sprintf("Package '%s' is needed for this function. Install it.", pkg),
      call. = FALSE
    )
  }
  invisible(TRUE)
}
