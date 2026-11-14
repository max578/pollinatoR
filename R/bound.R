# Partial identification under a cage effect (layer L4) ------------------------

#' Bound the loss share over a range of cage effects
#'
#' An exclusion device does more than exclude: a cage shades, shelters and
#' changes the microclimate, so the excluded arm measures the yield without
#' pollinators multiplied by an unknown cage factor `k`. If the cage
#' multiplies that yield by `k`, the observed loss share `D = 1 - low/high`
#' relates to the share net of the cage by
#'
#' `D_true = 1 - (1 - D) / k`.
#'
#' Without a sham arm `k` is not identified, but a range for it bounds the
#' share. The bound is taken per draw at both ends of the range: since
#' `D_true` rises with `k`, the lower bound is at `cage_factor[1]` and the
#' upper bound at `cage_factor[2]`. The result reports the median of each and
#' the outer interval, from the lower percentile at the low end to the upper
#' percentile at the high end.
#'
#' `cage_factor` has no default, because no value of it is neutral. A range
#' can come from a sham arm in another trial of the same crop and device, or
#' from a hand-pollinated arm grown under the device beside one grown without
#' it. A contrast of two devices that both exclude or both admit pollinators
#' does not estimate it. The range `c(1, 1)` assumes no cage effect and
#' returns the observed share.
#'
#' A design whose dependence reference is the sham arm is already net of the
#' cage and is refused. The bound assumes the reference arm carries no device.
#'
#' @param x A [pollin_draws()] object of a dependence estimand on the ratio
#'   scale with `basis = "loss"`.
#' @param cage_factor Two positive numbers, low then high: the range of the
#'   factor by which the device multiplies the yield without pollinators.
#'   Values above 1 mean the device raises it.
#' @param level The column of the draws to bound: a level of the fit's `by`,
#'   `"overall"` or `"population"`. The default is `"population"` when the
#'   draws hold it, else the only column.
#' @param conf_level Level of the outer interval.
#'
#' @return An object of class `pollin_bound`: a list holding `table` (one row
#'   with `level`, `k_low`, `k_high`, `observed` (the median observed share),
#'   `lower`, `upper` (the medians of the bounds), `conf.low` and `conf.high`
#'   (the outer interval) and `n`), the bound `draws` (a matrix with columns
#'   `lower` and `upper`), `cage_factor`, `conf_level` and the `estimand`.
#' @seealso [pollin_draws()]; a design with a `sham` arm and
#'   `reference = "sham"` estimates the share net of the cage directly.
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- data.frame(
#'   trt = rep(rep(c("Cage", "Open"), each = 10L), 8L),
#'   Site = rep(paste0("s", seq_len(8L)), each = 20L)
#' )
#' d$seed <- rgamma(nrow(d), 20, 20 / exp(1 + 0.3 * (d$trt == "Open")))
#' pd <- pollin_data(d, "seed", "trt", unit = "Site")
#' ds <- pollin_design(pd, c(excluded = "Cage", open = "Open"), "rcbd")
#' est <- pollin_estimand(ds, "dependence", "ratio", "loss")
#' dr <- pollin_draws(pollin_fit(est), n = 2000L, seed = 1L)
#' # A cage that raises the pollinator-free yield by 0 to 20 percent
#' pollin_bound(dr, cage_factor = c(1, 1.2))
#' }
pollin_bound <- function(x, cage_factor, level = NULL, conf_level = 0.95) {
  .check_class(x, "pollin_draws")
  .bound_check_estimand(x$estimand)
  .bound_check_factor(if (!missing(cage_factor)) cage_factor)
  .check_conf_level(conf_level)
  level <- .value_level(level, colnames(x$estimate))

  d_obs <- 1 - x$low[, level] / x$high[, level]
  draws <- cbind(
    lower = 1 - (1 - d_obs) / cage_factor[1L],
    upper = 1 - (1 - d_obs) / cage_factor[2L]
  )
  a <- (1 - conf_level) / 2
  table <- data.frame(
    level = level,
    k_low = cage_factor[1L],
    k_high = cage_factor[2L],
    observed = stats::median(d_obs),
    lower = stats::median(draws[, "lower"]),
    upper = stats::median(draws[, "upper"]),
    conf.low = stats::quantile(draws[, "lower"], a, names = FALSE),
    conf.high = stats::quantile(draws[, "upper"], 1 - a, names = FALSE),
    n = nrow(draws),
    stringsAsFactors = FALSE
  )
  structure(
    list(
      table = table,
      draws = draws,
      cage_factor = cage_factor,
      conf_level = conf_level,
      estimand = x$estimand
    ),
    class = "pollin_bound"
  )
}

#' @export
print.pollin_bound <- function(x, ...) {
  t <- x$table
  cat("<pollin_bound>  loss share net of a cage factor in [",
    format(t$k_low), ", ", format(t$k_high), "]\n", sep = "")
  cat(sprintf("  observed share : %.4f\n", t$observed))
  cat(sprintf("  bounds         : %.4f to %.4f (medians)\n", t$lower, t$upper))
  cat(sprintf(
    "  outer interval : %.4f to %.4f (%s%%)\n",
    t$conf.low,
    t$conf.high,
    format(100 * x$conf_level)
  ))
  cat("  draws          :", t$n, "\n")
  invisible(x)
}

# Internal checks -------------------------------------------------------------

#' Stop unless the estimand is a loss share whose reference is not the sham
#' @noRd
#' @keywords internal
.bound_check_estimand <- function(estimand) {
  if (
    estimand$family != "dependence" ||
      estimand$scale != "ratio" ||
      !identical(estimand$basis, "loss")
  ) {
    stop(
      "A cage bound applies to the loss share D = 1 - low/high; declare it ",
      "with `pollin_estimand(design, \"dependence\", \"ratio\", \"loss\")`.",
      call. = FALSE
    )
  }
  if (identical(estimand$design$reference, "sham")) {
    stop(
      "The estimand's reference is the sham arm, so its share is already net ",
      "of the cage effect; there is nothing to bound.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Stop unless `cage_factor` is two finite positive numbers, low then high
#' @noRd
#' @keywords internal
.bound_check_factor <- function(cage_factor) {
  if (
    !is.numeric(cage_factor) ||
      length(cage_factor) != 2L ||
      any(!is.finite(cage_factor)) ||
      any(cage_factor <= 0) ||
      cage_factor[1L] > cage_factor[2L]
  ) {
    stop(
      "`cage_factor` must be two positive numbers, low then high (for ",
      "example c(1, 1.2)). It has no default: no cage factor is neutral.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
