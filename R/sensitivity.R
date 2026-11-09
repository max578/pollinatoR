# Sensitivity to unmeasured confounding (layer L4) ----------------------------

#' Sensitivity of an effect to unmeasured confounding
#'
#' For a treatment contrast read by [pollin_effects()], returns the E-value
#' (VanderWeele and Ding 2017): the minimum strength of association, on the
#' risk ratio scale, that an unmeasured confounder would need with both the
#' treatment and the response to explain the contrast away. It is given for
#' the point estimate and for the interval limit closer to the null (1 when
#' the interval contains the null). The ratio of arm means is recovered from
#' the estimand's scale, and the E-value is that of the ratio or of its
#' inverse, whichever is at least 1.
#'
#' The bounding-factor result behind the E-value holds for a risk ratio or a
#' rate ratio. The ratio of arm means is a rate ratio for a count response
#' and a risk ratio for a binomial one, and only those are accepted. The
#' function refuses, each with its reason,
#'
#' * the difference scale, which carries no ratio;
#' * a direct effect, whose sensitivity also depends on unmeasured
#'   mediator-outcome confounding, which the E-value does not bound;
#' * a response fitted with a continuous family (for example Gamma or
#'   Gaussian), whose ratio of means is neither a risk nor a rate ratio.
#'
#' For a [pollin_moderation()] result no E-value is formed: an interaction
#' term is not a ratio of arm means. The output states, for each
#' treatment-by-exposure term, the additive bias on the link scale that would
#' move the estimate to zero and the bias that would move the interval limit
#' closer to zero onto it, beside the negative control of the moderation fit.
#'
#' @param x A [pollin_effects()] or [pollin_moderation()] object.
#' @param conf_level Level of the normal interval of a moderation term.
#'   Ignored for a `pollin_effects` object, whose interval is used as given.
#'
#' @return An object of class `pollin_sensitivity`: a list holding `type`
#'   (`"evalue"` or `"moderation"`) and `table`. For an E-value the table has
#'   `response`, any `by` column, `rr`, `rr.low`, `rr.high` (the ratio of arm
#'   means and its interval), `evalue` and `evalue.ci`. For a moderation
#'   result it has `term`, `component`, `estimate`, `conf.low`, `conf.high`,
#'   `bias.estimate` and `bias.ci` (link scale), with `negative_control` the
#'   moderation fit's negative-control table.
#' @references VanderWeele, T. J., & Ding, P. (2017). Sensitivity analysis in
#'   observational research: introducing the E-value. *Annals of Internal
#'   Medicine*, 167(4).
#' @seealso [pollin_bound()] for the cage effect, a bias that is not
#'   confounding.
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- data.frame(
#'   trt = rep(rep(c("Cage", "Open"), each = 10L), 8L),
#'   Site = rep(paste0("s", seq_len(8L)), each = 20L)
#' )
#' site <- rnorm(8L, 0, 0.3)[as.integer(factor(d$Site))]
#' d$pods <- rnbinom(
#'   nrow(d),
#'   mu = exp(3 + site + 0.3 * (d$trt == "Open")),
#'   size = 20
#' )
#' pd <- pollin_data(d, "pods", "trt", unit = "Site")
#' ds <- pollin_design(pd, c(excluded = "Cage", open = "Open"), "rcbd")
#' est <- pollin_estimand(ds, "dependence", "ratio", "loss")
#' pollin_sensitivity(pollin_effects(pollin_fit(est)))
#' }
pollin_sensitivity <- function(x, conf_level = 0.95) {
  if (inherits(x, "pollin_moderation")) {
    .check_conf_level(conf_level)
    return(.sensitivity_moderation(x, conf_level))
  }
  if (!inherits(x, "pollin_effects")) {
    stop(
      "`x` must be a `pollin_effects` or `pollin_moderation` object.",
      call. = FALSE
    )
  }
  e <- x$estimand
  tab <- x$table
  .sensitivity_check(e, tab)

  rr <- .sensitivity_rr(tab$estimate, e$scale, e$basis)
  lo <- .sensitivity_rr(tab$conf.low, e$scale, e$basis)
  hi <- .sensitivity_rr(tab$conf.high, e$scale, e$basis)
  near <- ifelse(lo > 1, lo, ifelse(hi < 1, hi, 1))
  out <- cbind(
    tab[c("response", x$by)],
    data.frame(
      rr = rr,
      rr.low = lo,
      rr.high = hi,
      evalue = .evalue(rr),
      evalue.ci = .evalue(near),
      row.names = NULL
    )
  )
  structure(
    list(type = "evalue", table = out, estimand = e),
    class = "pollin_sensitivity"
  )
}

#' @export
print.pollin_sensitivity <- function(x, ...) {
  if (x$type == "evalue") {
    cat("<pollin_sensitivity>  E-value of the ratio of arm means\n")
    print(x$table, row.names = FALSE, digits = 4L)
    cat(
      "  evalue: the association of an unmeasured confounder with both the",
      "treatment\n  and the response, on the risk ratio scale, that moves the",
      "estimate to the null;\n  evalue.ci: the same for the interval limit",
      "closer to the null.\n"
    )
    return(invisible(x))
  }
  cat("<pollin_sensitivity>  link-scale bias to the null, moderation terms\n")
  print(x$table, row.names = FALSE, digits = 4L)
  if (!is.null(x$negative_control)) {
    cat("  negative control: exposure slopes in the excluded arm\n")
    print(x$negative_control, row.names = FALSE, digits = 4L)
  }
  invisible(x)
}

# Internal ---------------------------------------------------------------------

#' Distributions whose ratio of means is a rate ratio or a risk ratio
#' @noRd
#' @keywords internal
.sensitivity_families <- c(
  "poisson",
  "nbinom1",
  "nbinom2",
  "binomial",
  "betabinomial"
)

#' Stop unless the effects carry a ratio of arm means of a count or binomial
#' response under a total effect
#' @noRd
#' @keywords internal
.sensitivity_check <- function(estimand, tab) {
  if (estimand$scale == "difference") {
    stop(
      "The difference scale carries no ratio, so it has no E-value; read ",
      "the effects on the `ratio`, `lnrr` or `rii` scale.",
      call. = FALSE
    )
  }
  if (identical(estimand$effect, "direct")) {
    stop(
      "A direct effect is refused: its bias also comes from unmeasured ",
      "confounding of the mediators and the response, which the E-value does ",
      "not bound. Use the total effect.",
      call. = FALSE
    )
  }
  bad <- unique(tab$distribution[
    !tab$distribution %in% .sensitivity_families
  ])
  if (length(bad)) {
    stop(
      "The E-value holds for a risk ratio or a rate ratio. A ratio of means ",
      "of a response fitted as ",
      paste(bad, collapse = ", "),
      " is neither; it is computed only for count and binomial responses.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' The ratio of arm means from a value on the estimand's scale
#' @noRd
#' @keywords internal
.sensitivity_rr <- function(v, scale, basis) {
  switch(
    scale,
    ratio = if (is.na(basis)) {
      v
    } else if (basis == "loss") {
      1 / (1 - v)
    } else {
      v + 1
    },
    rii = (1 + v) / (1 - v),
    lnrr = exp(v)
  )
}

#' The E-value of a risk ratio: RR + sqrt(RR (RR - 1)) after taking the
#' ratio or its inverse, whichever is at least 1
#' @noRd
#' @keywords internal
.evalue <- function(rr) {
  r <- ifelse(rr < 1, 1 / rr, rr)
  r + sqrt(r * (r - 1))
}

#' Link-scale bias that moves each moderation term, and its interval limit
#' closer to zero, onto zero
#' @noRd
#' @keywords internal
.sensitivity_moderation <- function(x, conf_level) {
  t <- x$table
  t <- t[t$component != "adjustment", , drop = FALSE]
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  lo <- t$estimate - z * t$std.error
  hi <- t$estimate + z * t$std.error
  out <- data.frame(
    term = t$term,
    component = t$component,
    estimate = t$estimate,
    conf.low = lo,
    conf.high = hi,
    bias.estimate = t$estimate,
    bias.ci = ifelse(lo > 0, lo, ifelse(hi < 0, hi, 0)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  structure(
    list(
      type = "moderation",
      table = out,
      negative_control = x$negative_control,
      conf_level = conf_level
    ),
    class = "pollin_sensitivity"
  )
}
