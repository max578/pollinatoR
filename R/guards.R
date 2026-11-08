# G4 guards (layer L5) --------------------------------------------------------

#' Guard against averaging bounded unit-level indices
#'
#' The Relative Interaction Index and the dependence ratio are bounded ratios.
#' The mean of per-unit (per-plant, per-plot, per-site) ratios is not the ratio
#' of the expected values: with unequal unit means the two differ, and the
#' mean of ratios has no causal reading as the trial's dependence. This guard
#' refuses the average unless `allow = TRUE`, and when allowed it reports the
#' average with that statement and, from a design, beside the ratio of means.
#'
#' Two inputs are accepted:
#'
#' * a [pollin_design()]: unit-level arm means are computed as in
#'   [pollin_rii()] (excluded versus open arm, outermost unit by default), the
#'   index is computed for each unit, and the ratio of means is computed from
#'   the arm means across units;
#' * a numeric vector of unit-level index values already computed. The ratio
#'   of means cannot be recovered from these and is reported as `NA`.
#'
#' The indices are `rii = (open - excluded) / (open + excluded)` and
#' `dependence = 1 - excluded / open`, as in [pollin_index()].
#'
#' @param x A [pollin_design()] with an `"excluded"` and an `"open"` arm, or a
#'   numeric vector of unit-level index values.
#' @param index `"rii"` or `"dependence"`. There is no default.
#' @param allow `FALSE` (the default) stops with an error; `TRUE` returns the
#'   average with the warning statement attached.
#' @param response,unit With a design, the response and the unit level, as in
#'   [pollin_rii()].
#'
#' @return With `allow = TRUE`, an object of class `pollin_bounded`: a list
#'   holding `index`, `n` (units averaged), `mean_of_units`, `ratio_of_means`
#'   (`NA` for a numeric vector), `per_unit` (a data frame, or `NULL` for a
#'   numeric vector) and the `note` printed with it.
#' @seealso [pollin_index()] for the index of the arm means with an interval.
#' @export
#' @examples
#' df <- data.frame(
#'   pods = rpois(48, rep(c(20, 30, 40, 60), each = 12)),
#'   trt = rep(c("Cage", "Open"), 24),
#'   Site = rep(letters[1:4], each = 12)
#' )
#' pd <- pollin_data(df, "pods", "trt", unit = "Site")
#' ds <- pollin_design(pd, c(excluded = "Cage", open = "Open"), "multi_site")
#' # Refused by default
#' try(pollin_guard_bounded(ds, index = "rii"))
#' # Allowed, with the statement and the ratio of means beside it
#' pollin_guard_bounded(ds, index = "rii", allow = TRUE)
pollin_guard_bounded <- function(
  x,
  index,
  allow = FALSE,
  response = NULL,
  unit = NULL
) {
  .check_choice(if (!missing(index)) index, "index", c("rii", "dependence"))
  .check_flag(allow, "allow")

  # Unit-level values ----------------------------------------------------------
  if (inherits(x, "pollin_design")) {
    m <- .unit_arm_means(x, response, unit)
    value <- .bounded_index(m$open, m$excluded, index)
    ratio <- .bounded_index(mean(m$open), mean(m$excluded), index)
    per_unit <- data.frame(
      unit = m$unit,
      open = m$open,
      excluded = m$excluded,
      value = value,
      stringsAsFactors = FALSE
    )
  } else {
    .check_unit_values(x, index)
    value <- x
    ratio <- NA_real_
    per_unit <- NULL
  }
  n_bad <- sum(!is.finite(value))
  if (n_bad > 0L) {
    warning(
      .n_noun(n_bad, "unit"),
      " with an undefined index left out of the average.",
      call. = FALSE
    )
  }
  value <- value[is.finite(value)]
  label <- if (index == "rii") "RII" else "dependence"

  # Refuse or report -----------------------------------------------------------
  note <- paste0(
    "Averaging bounded per-unit ratios is not the ratio of expected values; ",
    "the mean of ",
    .n_noun(length(value), paste("unit-level", label, "value")),
    " is not the trial's ",
    label,
    "."
  )
  if (!allow) {
    stop(
      note,
      if (is.finite(ratio)) {
        sprintf(" The %s of the arm means is %.4f.", label, ratio)
      },
      " Set `allow = TRUE` to average them anyway.",
      call. = FALSE
    )
  }
  structure(
    list(
      index = index,
      n = length(value),
      mean_of_units = mean(value),
      ratio_of_means = ratio,
      per_unit = per_unit,
      note = note
    ),
    class = "pollin_bounded"
  )
}

#' @export
print.pollin_bounded <- function(x, ...) {
  label <- if (x$index == "rii") "RII" else "dependence"
  cat("<pollin_bounded>  average of unit-level", label, "values\n")
  cat(sprintf(
    "  mean of %s : %.4f\n",
    .n_noun(x$n, "unit"),
    x$mean_of_units
  ))
  if (is.finite(x$ratio_of_means)) {
    cat(sprintf("  %s of arm means : %.4f\n", label, x$ratio_of_means))
  } else {
    cat("  ratio of means : not available from unit-level values alone\n")
  }
  cat(" ", x$note, "\n")
  invisible(x)
}

#' Guard against reading a between-unit slope as a within-unit one
#'
#' An exposure slope estimated across sites (between units) is an ecological
#' association: it need not hold for plants within a site, where unit-level
#' confounders are held fixed. This guard splits the exposure into its
#' between-unit and within-unit parts with [pollin_within_between()] and fits
#' `response ~ within + between` by least squares, reporting the pooled,
#' within-unit and between-unit slopes. It warns when
#'
#' * the exposure does not vary within units, so only the between-unit slope
#'   exists;
#' * only the between-unit slope is distinguishable from zero at `alpha`;
#' * the within-unit and between-unit slopes are both distinguishable from
#'   zero at `alpha` and disagree in sign.
#'
#' The guard works on the rows given; subset to one arm beforehand when the
#' slope of interest is within that arm.
#'
#' @param x A data frame, or a [pollin_data()] object.
#' @param exposure Name of the numeric exposure column (for example insect
#'   abundance).
#' @param response Name of the numeric response column. With a [pollin_data()]
#'   object it may be left `NULL` when one response is declared.
#' @param unit The unit column(s) the exposure is split by. Required for a data
#'   frame (several columns are combined); with a [pollin_data()] object, one of
#'   the unit columns, taken within its outer units, and `NULL` (the default)
#'   uses the outermost unit.
#' @param alpha Significance level for the third check.
#'
#' @return An object of class `pollin_ecological`: a list holding the slope
#'   table `slopes` (component, slope, standard error, p-value), the number of
#'   units `n_units`, the `flags` raised (a character vector, empty when none)
#'   and the column names used.
#' @seealso [pollin_moderation()] for the model-based within/between split.
#' @export
#' @examples
#' set.seed(1)
#' site <- rep(sprintf("s%02d", 1:10), each = 8)
#' site_bees <- rep(rnorm(10, 10, 3), each = 8)
#' bees <- site_bees + rnorm(80)
#' # Sites with more bees yield more, but plants within a site do not respond
#' pods <- 20 + 2 * site_bees + rnorm(80, 0, 2)
#' df <- data.frame(pods = pods, bees = bees, Site = site)
#' pollin_guard_ecological(df, "bees", "pods", unit = "Site")
pollin_guard_ecological <- function(
  x,
  exposure,
  response = NULL,
  unit = NULL,
  alpha = 0.05
) {
  # Columns from the input ----------------------------------------------------
  if (inherits(x, "pollin_data")) {
    roles <- x$roles
    response <- .resolve_response(response, roles$response)
    level <- .resolve_unit_level(unit, roles$unit)
    unit <- roles$unit[seq_len(match(level, roles$unit))]
    d <- x$data
  } else {
    .check_data_frame(x)
    d <- x
  }
  .check_ecological_columns(d, exposure, response, unit)
  .check_alpha(alpha)

  # Within and between slopes -------------------------------------------------
  key <- interaction(d[unit], sep = ":", drop = TRUE)
  ok <- is.finite(d[[exposure]]) & is.finite(d[[response]]) & !is.na(key)
  y <- d[[response]][ok]
  e <- d[[exposure]][ok]
  key <- droplevels(key[ok])
  n_units <- nlevels(key)
  if (n_units < 3L) {
    stop(
      "Need at least three units to estimate a between-unit slope; found ",
      n_units,
      ".",
      call. = FALSE
    )
  }
  wb <- pollin_within_between(e, key)
  within <- wb$within
  between <- wb$between
  no_within <- all(abs(within) < 1e-8 * max(1, abs(e)))
  pooled <- .slope_row(stats::lm(y ~ e), "e", "pooled")
  if (no_within) {
    fit <- stats::lm(y ~ between)
    slopes <- rbind(
      pooled,
      data.frame(
        component = "within",
        slope = NA_real_,
        std.error = NA_real_,
        p.value = NA_real_
      ),
      .slope_row(fit, "between", "between")
    )
  } else {
    fit <- stats::lm(y ~ within + between)
    slopes <- rbind(
      pooled,
      .slope_row(fit, "within", "within"),
      .slope_row(fit, "between", "between")
    )
  }

  # Flags ----------------------------------------------------------------------
  b_w <- slopes$slope[2L]
  b_b <- slopes$slope[3L]
  p_w <- slopes$p.value[2L]
  p_b <- slopes$p.value[3L]
  flags <- character(0)
  if (no_within) {
    flags <- c(
      flags,
      paste0(
        "`",
        exposure,
        "` does not vary within units; its slope is between-unit only and ",
        "cannot be read at plant level."
      )
    )
  } else {
    if (is.finite(p_b) && is.finite(p_w) && p_b < alpha && p_w >= alpha) {
      flags <- c(
        flags,
        sprintf(
          paste0(
            "Only the between-unit slope is distinguishable from zero ",
            "(p = %.3g between, %.3g within); the association is ",
            "between units, not within them."
          ),
          p_b,
          p_w
        )
      )
    }
    both <- is.finite(p_w) && is.finite(p_b) && p_w < alpha && p_b < alpha
    if (both && sign(b_w) * sign(b_b) < 0) {
      flags <- c(
        flags,
        sprintf(
          paste0(
            "The within-unit slope (%.4g) and the between-unit slope (%.4g) ",
            "disagree in sign."
          ),
          b_w,
          b_b
        )
      )
    }
  }
  for (f in flags) {
    warning(f, call. = FALSE)
  }
  structure(
    list(
      slopes = slopes,
      n_units = n_units,
      flags = flags,
      exposure = exposure,
      response = response,
      unit = unit
    ),
    class = "pollin_ecological"
  )
}

#' @export
print.pollin_ecological <- function(x, ...) {
  cat(
    "<pollin_ecological>  ",
    x$response,
    " ~ ",
    x$exposure,
    " split by ",
    paste(x$unit, collapse = ":"),
    " (",
    .n_noun(x$n_units, "unit"),
    ")\n",
    sep = ""
  )
  print(x$slopes, row.names = FALSE, digits = 4L)
  if (length(x$flags)) {
    cat(paste0("  ! ", x$flags, "\n"), sep = "")
  } else {
    cat("  within-unit and between-unit slopes agree\n")
  }
  invisible(x)
}

#' Guard against adjusting for a covariate downstream of treatment
#'
#' Adjusting a treatment contrast for a variable the treatment itself changes
#' (for example vegetative weight, or a pollinator visit count under a cage)
#' blocks part of the effect and can open a collider path. Whether a covariate
#' is downstream of treatment is a fact about the system that the data cannot
#' reveal. The declaration is the user's: `downstream` names the variables the
#' user holds to be affected by treatment, and the guard stops if any
#' adjustment covariate in `x` is among them.
#'
#' @param x A [pollin_data()], [pollin_design()] or [pollin_estimand()]
#'   object.
#' @param downstream Character vector of variables declared downstream of
#'   treatment. `character(0)` declares none. There is no default.
#'
#' @return `x`, invisibly, when no adjustment covariate is declared downstream.
#' @export
#' @examples
#' df <- data.frame(
#'   pods = rpois(40, 30),
#'   trt = rep(c("Cage", "Open"), 20),
#'   Site = rep(letters[1:4], each = 10),
#'   rain = rep(c(250, 300, 320, 410), each = 10),
#'   veg_weight = rnorm(40, 100, 10)
#' )
#' pd <- pollin_data(df, "pods", "trt", unit = "Site", adjust = "rain")
#' pollin_guard_posttreatment(pd, downstream = "veg_weight")
#'
#' pd2 <- pollin_data(df, "pods", "trt", unit = "Site",
#'                    adjust = c("rain", "veg_weight"))
#' try(pollin_guard_posttreatment(pd2, downstream = "veg_weight"))
pollin_guard_posttreatment <- function(x, downstream) {
  roles <- .roles_of(x)
  .check_downstream(if (!missing(downstream)) downstream)
  hit <- intersect(roles$adjust, downstream)
  if (length(hit)) {
    stop(
      "Adjustment ",
      if (length(hit) == 1L) "covariate " else "covariates ",
      paste0("`", hit, "`", collapse = ", "),
      if (length(hit) == 1L) " is" else " are",
      " declared downstream of `",
      roles$treatment,
      "`. Adjusting for a downstream covariate removes part of the treatment ",
      "effect; drop it from `adjust`.",
      call. = FALSE
    )
  }
  invisible(x)
}

# Internal helpers ------------------------------------------------------------

#' Bounded index from open and excluded values, elementwise; non-finite where
#' the denominator is zero
#' @noRd
#' @keywords internal
.bounded_index <- function(open, excluded, index) {
  out <- if (index == "rii") {
    (open - excluded) / (open + excluded)
  } else {
    1 - excluded / open
  }
  out[!is.finite(out)] <- NA_real_
  out
}

#' One row of a slope table from an lm fit
#' @noRd
#' @keywords internal
.slope_row <- function(fit, term, component) {
  cf <- summary(fit)$coefficients
  if (!term %in% rownames(cf)) {
    return(data.frame(
      component = component,
      slope = NA_real_,
      std.error = NA_real_,
      p.value = NA_real_
    ))
  }
  data.frame(
    component = component,
    slope = cf[term, 1L],
    std.error = cf[term, 2L],
    p.value = cf[term, 4L],
    row.names = NULL
  )
}

#' The role map of a pollin_data, pollin_design or pollin_estimand object
#' @noRd
#' @keywords internal
.roles_of <- function(x) {
  if (inherits(x, "pollin_estimand")) {
    x <- x$design
  }
  if (inherits(x, "pollin_design")) {
    x <- x$data
  }
  if (!inherits(x, "pollin_data")) {
    stop(
      "`x` must be a `pollin_data`, `pollin_design` or `pollin_estimand` ",
      "object.",
      call. = FALSE
    )
  }
  x$roles
}

# Internal checks -------------------------------------------------------------

#' Stop unless `value` is TRUE or FALSE
#' @noRd
#' @keywords internal
.check_flag <- function(value, arg) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop("`", arg, "` must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Stop unless `x` is a numeric vector of index values within the index bounds
#' @noRd
#' @keywords internal
.check_unit_values <- function(x, index) {
  if (!is.numeric(x) || length(x) == 0L) {
    stop(
      "`x` must be a `pollin_design` object or a numeric vector of ",
      "unit-level index values.",
      call. = FALSE
    )
  }
  v <- x[is.finite(x)]
  out <- if (index == "rii") abs(v) > 1 else v > 1
  if (any(out)) {
    stop(
      .n_noun(sum(out), "value"),
      " outside the bounds of the ",
      index,
      " index (",
      if (index == "rii") "-1 to 1" else "at most 1",
      ").",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Stop unless exposure, response and unit name columns of `d`, the first two
#' numeric
#' @noRd
#' @keywords internal
.check_ecological_columns <- function(d, exposure, response, unit) {
  for (nm in c("exposure", "response")) {
    v <- get(nm)
    if (!is.character(v) || length(v) != 1L) {
      stop("`", nm, "` must name one column.", call. = FALSE)
    }
  }
  if (!is.character(unit) || length(unit) == 0L) {
    stop("`unit` must name the unit column.", call. = FALSE)
  }
  miss <- setdiff(c(exposure, response, unit), names(d))
  if (length(miss)) {
    stop(
      "Columns not found in `x`: ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }
  is_num <- vapply(d[c(exposure, response)], is.numeric, logical(1))
  if (!all(is_num)) {
    stop(
      "Columns must be numeric: ",
      paste(c(exposure, response)[!is_num], collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Stop unless `alpha` is a single number in (0, 1)
#' @noRd
#' @keywords internal
.check_alpha <- function(alpha) {
  if (!is.numeric(alpha) || length(alpha) != 1L || !(alpha > 0 && alpha < 1)) {
    stop("`alpha` must be a single number between 0 and 1.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Stop unless `downstream` is a character vector (possibly empty) without NA
#' @noRd
#' @keywords internal
.check_downstream <- function(downstream) {
  if (!is.character(downstream) || anyNA(downstream)) {
    stop(
      "`downstream` must be a character vector of variable names; ",
      "use character(0) to declare none.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
