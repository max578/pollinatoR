# Uncertainty-aware dependence indices ----------------------------------------

#' Relative Interaction Index with a bootstrap interval
#'
#' The Relative Interaction Index (Armas et al. 2004) is widely reported for
#' pollinator dependence as a bare point estimate, despite being a ratio of two
#' correlated sums whose sampling distribution is awkward. This function returns
#' the index with a nonparametric bootstrap interval, the minimum upgrade over a
#' point value. Prefer the model-based contribution
#' ([pollin_contribution()]) as the primary estimate and keep this for
#' comparability with the existing literature.
#'
#' The arm values are given either as two vectors or as a [pollin_design()].
#' From a design, the values are unit-level arm means: the mean response of
#' each unit in the open arm and in the excluded arm, for the units carrying
#' both arms, in sorted unit order. The bootstrap then runs on those means as it
#' does on two vectors.
#'
#' @param open Numeric vector of the open-arm (pollinator-accessible) values,
#'   or a [pollin_design()] with an `"excluded"` and an `"open"` arm.
#' @param cage Numeric vector of the exclusion-arm values. Leave it empty when
#'   `open` is a design.
#' @param n_boot Number of bootstrap resamples.
#' @param conf Confidence level.
#' @param seed Optional integer seed for reproducibility.
#' @param response With a design, the response column to use. It may be left
#'   `NULL` when the data declare one response.
#' @param unit With a design, the unit column the arm means are taken over, one
#'   of the design's unit nesting. `NULL` (the default) uses the outermost
#'   unit; an inner unit is taken within its outer units.
#'
#' @return An object of class `pollin_rii`.
#' @export
#' @examples
#' set.seed(1)
#' pollin_rii(rpois(50, 8), rpois(50, 6), n_boot = 500, seed = 1)
#'
#' # From a design: one mean per site and arm
#' df <- data.frame(
#'   pods = rpois(48, 30),
#'   trt = rep(c("Cage", "Open"), 24),
#'   Site = rep(letters[1:4], each = 12)
#' )
#' pd <- pollin_data(df, "pods", "trt", unit = "Site")
#' ds <- pollin_design(pd, c(excluded = "Cage", open = "Open"), "multi_site")
#' pollin_rii(ds, n_boot = 200, seed = 1)
pollin_rii <- function(
  open,
  cage,
  n_boot = 2000L,
  conf = 0.95,
  seed = NULL,
  response = NULL,
  unit = NULL
) {
  if (inherits(open, "pollin_design")) {
    .check_no_cage(missing(cage))
    m <- .unit_arm_means(open, response, unit)
    open <- m$open
    cage <- m$excluded
  }
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
#' Given a [pollin_design()], the indices are computed from unit-level arm
#' means, as described in [pollin_rii()].
#'
#' @inheritParams pollin_rii
#'
#' @return A data frame of indices with estimates and intervals where available.
#' @export
#' @examples
#' set.seed(1)
#' pollin_index(rpois(50, 9), rpois(50, 6), n_boot = 500, seed = 1)
#'
#' # From a design with two responses: name the one to use
#' df <- data.frame(
#'   pods = rpois(48, 30),
#'   seed_wt = rgamma(48, 20, 5),
#'   trt = rep(c("Cage", "Open"), 24),
#'   Site = rep(letters[1:4], each = 12)
#' )
#' pd <- pollin_data(df, c("pods", "seed_wt"), "trt", unit = "Site")
#' ds <- pollin_design(pd, c(excluded = "Cage", open = "Open"), "multi_site")
#' pollin_index(ds, response = "pods", n_boot = 200, seed = 1)
pollin_index <- function(
  open,
  cage,
  n_boot = 2000L,
  conf = 0.95,
  seed = NULL,
  response = NULL,
  unit = NULL
) {
  if (inherits(open, "pollin_design")) {
    .check_no_cage(missing(cage))
    m <- .unit_arm_means(open, response, unit)
    open <- m$open
    cage <- m$excluded
  }
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

# Unit-level arm means from a design ------------------------------------------

#' Mean response per unit and arm, for the units carrying both the excluded
#' and the open arm. Returns `unit` (the unit labels, sorted), `open` and
#' `excluded` (numeric, aligned with `unit`), and the `response` and `level`
#' used.
#' @noRd
#' @keywords internal
.unit_arm_means <- function(design, response = NULL, unit = NULL) {
  .check_class(design, "pollin_design")
  .check_arm_pair(design$arms)
  roles <- design$data$roles
  response <- .resolve_response(response, roles$response)
  level <- .resolve_unit_level(unit, roles$unit)

  # Means per unit and arm -----------------------------------------------------
  d <- design$data$data
  nest <- roles$unit[seq_len(match(level, roles$unit))]
  key <- interaction(d[nest], sep = ":", drop = TRUE, lex.order = TRUE)
  y <- d[[response]]
  trt <- as.character(d[[roles$treatment]])
  arm_mean <- function(lev) {
    keep <- trt %in% lev & is.finite(y)
    m <- tapply(y[keep], droplevels(key[keep]), mean)
    m[!is.na(m)]
  }
  m_open <- arm_mean(design$arms[["open"]])
  m_excl <- arm_mean(design$arms[["excluded"]])

  # Units with both arms -------------------------------------------------------
  both <- intersect(levels(key), intersect(names(m_open), names(m_excl)))
  n_drop <- length(union(names(m_open), names(m_excl))) - length(both)
  if (n_drop > 0L) {
    warning(
      .n_noun(n_drop, "unit"),
      " of `",
      paste(nest, collapse = ":"),
      "` with only one of the open and excluded arms left out.",
      call. = FALSE
    )
  }
  if (length(both) < 2L) {
    stop(
      "Need at least two units of `",
      paste(nest, collapse = ":"),
      "` carrying both the open and the excluded arm; found ",
      length(both),
      ".",
      call. = FALSE
    )
  }
  list(
    unit = both,
    open = unname(as.numeric(m_open[both])),
    excluded = unname(as.numeric(m_excl[both])),
    response = response,
    level = level
  )
}

# Internal checks -------------------------------------------------------------

#' Stop when `cage` is supplied alongside a design
#' @noRd
#' @keywords internal
.check_no_cage <- function(cage_missing) {
  if (!cage_missing) {
    stop(
      "With a `pollin_design`, leave `cage` empty and name the response ",
      "with `response =`.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Stop unless the design has both an excluded and an open arm
#' @noRd
#' @keywords internal
.check_arm_pair <- function(arms) {
  absent <- setdiff(c("excluded", "open"), names(arms))
  if (length(absent)) {
    stop(
      "The design needs an excluded and an open arm; missing: ",
      paste(absent, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' One response column: the given name, or the only declared response
#' @noRd
#' @keywords internal
.resolve_response <- function(response, declared) {
  if (is.null(response)) {
    if (length(declared) > 1L) {
      stop(
        "The data declare several responses; name one with `response =`: ",
        paste(declared, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    return(declared)
  }
  if (
    !is.character(response) ||
      length(response) != 1L ||
      !response %in% declared
  ) {
    stop(
      "`response` must be one of the declared responses: ",
      paste(declared, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  response
}

#' One unit level: the given name, or the outermost unit
#' @noRd
#' @keywords internal
.resolve_unit_level <- function(unit, nesting) {
  if (is.null(unit)) {
    return(nesting[1L])
  }
  if (!is.character(unit) || length(unit) != 1L || !unit %in% nesting) {
    stop(
      "`unit` must be one of the unit columns: ",
      paste(nesting, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  unit
}
