# Uncertainty-aware dependence indices ----------------------------------------

#' Relative Interaction Index with a bootstrap interval
#'
#' The Relative Interaction Index (Armas et al. 2004) is widely reported for
#' pollinator dependence as a bare point estimate, despite being a ratio of two
#' correlated sums whose sampling distribution is awkward. This function returns
#' the index with a nonparametric bootstrap interval -- the minimum honest
#' upgrade over a point value. Prefer the model-based contribution
#' ([pollin_contribution()]) as the primary estimate and keep this for
#' comparability with the existing literature.
#'
#' @param open Numeric vector of the open-arm (pollinator-accessible) values.
#' @param cage Numeric vector of the exclusion-arm values.
#' @param n_boot Number of bootstrap resamples.
#' @param conf Confidence level.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return An object of class `pollin_rii`.
#' @export
#' @examples
#' set.seed(1)
#' pollin_rii(rpois(50, 8), rpois(50, 6), n_boot = 500, seed = 1)
pollin_rii <- function(open, cage, n_boot = 2000L, conf = 0.95, seed = NULL) {
  open <- open[is.finite(open)]
  cage <- cage[is.finite(cage)]
  rii <- function(o, c) {
    d <- mean(o) + mean(c)
    if (d == 0) {
      return(NA_real_)
    }
    (mean(o) - mean(c)) / d
  }
  est <- rii(open, cage)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  bo <- vapply(
    seq_len(n_boot),
    function(i) {
      rii(
        sample(open, length(open), replace = TRUE),
        sample(cage, length(cage), replace = TRUE)
      )
    },
    numeric(1)
  )
  a <- (1 - conf) / 2
  ci <- stats::quantile(bo, c(a, 1 - a), na.rm = TRUE)
  structure(
    list(
      rii = est,
      conf.low = unname(ci[1]),
      conf.high = unname(ci[2]),
      boot_sd = stats::sd(bo, na.rm = TRUE),
      n_open = length(open),
      n_cage = length(cage),
      n_boot = n_boot,
      conf = conf
    ),
    class = "pollin_rii"
  )
}

#' @export
print.pollin_rii <- function(x, ...) {
  cat(sprintf(
    "RII = %.3f  [%.0f%% boot CI: %.3f, %.3f]  (n_open=%d, n_cage=%d)\n",
    x$rii,
    100 * x$conf,
    x$conf.low,
    x$conf.high,
    x$n_open,
    x$n_cage
  ))
  invisible(x)
}

#' A menu of uncertainty-aware dependence indices
#'
#' Returns several dependence indices for an open-versus-exclusion contrast,
#' each with an interval where one is available: the Relative Interaction Index
#' (bootstrap interval), the log response ratio (closed-form interval, the
#' meta-analytic workhorse), and the proportional dependence ratio. Provided so
#' an analysis can report a defensible, interval-bearing index instead of a bare
#' ratio.
#'
#' @inheritParams pollin_rii
#'
#' @return A data frame of indices with estimates and intervals where available.
#' @export
#' @examples
#' set.seed(1)
#' pollin_index(rpois(50, 9), rpois(50, 6), n_boot = 500, seed = 1)
pollin_index <- function(open, cage, n_boot = 2000L, conf = 0.95, seed = NULL) {
  rii <- pollin_rii(open, cage, n_boot = n_boot, conf = conf, seed = seed)
  o <- open[is.finite(open)]
  c <- cage[is.finite(cage)]
  mo <- mean(o)
  mc <- mean(c)
  z <- stats::qnorm(1 - (1 - conf) / 2)
  if (mo > 0 && mc > 0) {
    lnrr <- log(mo / mc)
    se <- sqrt(
      stats::var(o) / (length(o) * mo^2) + stats::var(c) / (length(c) * mc^2)
    )
    lnrr_lo <- lnrr - z * se
    lnrr_hi <- lnrr + z * se
  } else {
    lnrr <- lnrr_lo <- lnrr_hi <- NA_real_
  }
  dep <- if (mo != 0) (mo - mc) / mo else NA_real_
  data.frame(
    index = c("RII", "lnRR", "dependence_ratio"),
    estimate = c(rii$rii, lnrr, dep),
    conf.low = c(rii$conf.low, lnrr_lo, NA_real_),
    conf.high = c(rii$conf.high, lnrr_hi, NA_real_),
    note = c(
      "bootstrap CI",
      "closed-form CI",
      "point; do not regress on open yield (see pollin_guard_spurious)"
    ),
    stringsAsFactors = FALSE
  )
}
