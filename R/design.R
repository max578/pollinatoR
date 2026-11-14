# Experimental design (layer L1) ----------------------------------------------

#' Declare the design of a pollination experiment
#'
#' Names the role of each treatment level (the arms) and the layout of the
#' experiment. The layout fixes the random structure used when models are
#' fitted: one random intercept for each level of the unit nesting declared in
#' [pollin_data()], outer to inner, with an optional random treatment slope on
#' one unit level, giving a pollination effect that varies between those units.
#'
#' | Layout | Units, outer to inner | Treatment varies | Default slope |
#' |---|---|---|---|
#' | `"crd"` | experimental unit (and subsamples within it) | between units | none |
#' | `"rcbd"` | block, then any inner levels | within blocks | none |
#' | `"split_plot"` | block, whole plot, then any inner levels | within whole plots | none |
#' | `"matched_pairs"` | any outer levels, then the pair | within pairs | none |
#' | `"multi_site"` | site, then any inner levels | within sites | site |
#'
#' In a split-plot layout the pollination treatment is taken to be the
#' sub-plot factor (for example cages within variety plots). The design checks
#' that every arm appears in every block, whole plot, pair or site, and warns
#' with a count when it does not.
#'
#' @param x A [pollin_data()] object.
#' @param arms Named character vector mapping arm roles to treatment levels, for
#'   example `c(excluded = "Cage", open = "Open")`. Roles are `"excluded"`,
#'   `"open"`, `"supplemented"` and `"sham"`; every treatment level must be
#'   mapped. A sham arm carries a device identical to the exclusion device that
#'   admits pollinators (an open cage, a cage with its mesh removed); it needs
#'   an open arm beside it.
#' @param layout One of `"crd"`, `"rcbd"`, `"split_plot"`, `"matched_pairs"`
#'   or `"multi_site"`. There is no default.
#' @param reference The arm that dependence is measured against: `"open"` (the
#'   default), `"supplemented"` or `"sham"`. It must be one of the declared
#'   arms. Against `"sham"` the dependence is net of the cage effect.
#' @param slope The unit level carrying a random treatment slope. `NULL` (the
#'   default) uses the layout's default from the table above; `FALSE` fits no
#'   slope; a column name from the unit nesting puts the slope there.
#'
#' @return An object of class `pollin_design`: a list holding the
#'   [pollin_data()] object `data`, `arms`, `layout`, `reference`, `slope` and
#'   the random-effect terms `random`.
#' @export
#' @examples
#' df <- data.frame(
#'   pods = rpois(48, 30),
#'   trt = rep(c("Cage", "Open"), 24),
#'   Site = rep(letters[1:4], each = 12),
#'   Plot = rep(rep(c("p1", "p2"), each = 6), 4)
#' )
#' pd <- pollin_data(df, "pods", "trt", unit = c("Site", "Plot"))
#' pollin_design(pd, arms = c(excluded = "Cage", open = "Open"),
#'               layout = "multi_site")
pollin_design <- function(x, arms, layout, reference = "open", slope = NULL) {
  .check_class(x, "pollin_data")
  .check_layout(layout)
  .check_arms(arms, x)
  .check_reference(reference, arms)
  unit <- x$roles$unit
  .check_unit_depth(unit, layout)
  slope <- .resolve_slope(slope, unit, layout)

  # Structure of the treatment within the units --------------------------------
  d <- x$data
  trt <- x$roles$treatment
  if (layout == "crd") {
    n_mixed <- .n_units_varying(d, trt, unit)
    if (n_mixed > 0L) {
      stop(
        "In a `crd` layout the treatment is applied to whole units; found ",
        .n_noun(n_mixed, "innermost unit"),
        " carrying more than one treatment level.",
        call. = FALSE
      )
    }
  } else {
    key <- switch(
      layout,
      rcbd = unit[1L],
      split_plot = unit[1:2],
      matched_pairs = unit,
      multi_site = unit[1L]
    )
    n_short <- .n_units_missing_arm(d, trt, key, unname(arms))
    if (n_short > 0L) {
      warning(
        .n_noun(n_short, "unit"),
        " of `",
        paste(key, collapse = ":"),
        "` lacking at least one arm.",
        call. = FALSE
      )
    }
  }

  structure(
    list(
      data = x,
      arms = arms,
      layout = layout,
      reference = reference,
      slope = slope,
      random = .nested_terms(unit, trt, slope)
    ),
    class = "pollin_design"
  )
}

#' @export
print.pollin_design <- function(x, ...) {
  cat("<pollin_design>\n")
  cat("  layout    :", x$layout, "\n")
  cat(
    "  arms      :",
    paste0(names(x$arms), " = ", x$arms, collapse = ", "),
    "\n"
  )
  cat("  reference :", x$reference, "\n")
  cat("  units     :", paste(x$data$roles$unit, collapse = " / "), "\n")
  cat("  random    :", paste(x$random, collapse = " + "), "\n")
  invisible(x)
}

# Random-effect terms ---------------------------------------------------------

#' Nested random terms, outer to inner: a random intercept per level, or a
#' random treatment slope on the `slope` level.
#' @noRd
#' @keywords internal
.nested_terms <- function(unit, treatment, slope = NULL) {
  terms <- character(0)
  for (i in seq_along(unit)) {
    key <- paste(unit[seq_len(i)], collapse = ":")
    if (identical(unit[i], slope)) {
      terms <- c(terms, sprintf("(%s | %s)", treatment, key))
    } else {
      terms <- c(terms, sprintf("(1 | %s)", key))
    }
  }
  terms
}

# Internal checks -------------------------------------------------------------

#' A count with its noun, pluralised: "1 unit", "2 units"
#' @noRd
#' @keywords internal
.n_noun <- function(n, noun) {
  paste(n, if (n == 1L) noun else paste0(noun, "s"))
}

.arm_roles <- c("excluded", "open", "supplemented", "sham")
.layouts <- c("crd", "rcbd", "split_plot", "matched_pairs", "multi_site")

#' Stop unless `x` inherits from `class`
#' @noRd
#' @keywords internal
.check_class <- function(x, class, arg = "x") {
  if (!inherits(x, class)) {
    stop("`", arg, "` must be a `", class, "` object.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Stop unless `layout` is one of the supported layouts
#' @noRd
#' @keywords internal
.check_layout <- function(layout) {
  if (!is.character(layout) || length(layout) != 1L || !layout %in% .layouts) {
    stop(
      "`layout` must be one of: ",
      paste(.layouts, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Check the arm map: named by known roles, each role and level once, at least
#' two arms, and every treatment level in the data mapped to an arm
#' @noRd
#' @keywords internal
.check_arms <- function(arms, x) {
  if (!is.character(arms) || is.null(names(arms)) || anyNA(arms)) {
    stop(
      "`arms` must be a named character vector, for example ",
      "c(excluded = \"Cage\", open = \"Open\").",
      call. = FALSE
    )
  }
  bad <- setdiff(names(arms), .arm_roles)
  if (length(bad)) {
    stop(
      "Unknown arm roles: ",
      paste(bad, collapse = ", "),
      ". Use ",
      paste(.arm_roles, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (anyDuplicated(names(arms)) || anyDuplicated(arms)) {
    stop(
      "Each arm role and each treatment level may appear once.",
      call. = FALSE
    )
  }
  if (length(arms) < 2L) {
    stop("A design needs at least two arms.", call. = FALSE)
  }
  lev <- .treatment_levels(x$data[[x$roles$treatment]])
  absent <- setdiff(arms, lev)
  if (length(absent)) {
    stop(
      "Arm levels not found in `",
      x$roles$treatment,
      "`: ",
      paste(absent, collapse = ", "),
      call. = FALSE
    )
  }
  unmapped <- setdiff(lev, arms)
  if (length(unmapped)) {
    stop(
      "Treatment levels with no arm: ",
      paste(unmapped, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Stop unless `reference` is open, supplemented or sham and is a declared arm
#' @noRd
#' @keywords internal
.check_reference <- function(reference, arms) {
  ok <- c("open", "supplemented", "sham")
  if (
    !is.character(reference) || length(reference) != 1L || !reference %in% ok
  ) {
    stop(
      "`reference` must be \"open\", \"supplemented\" or \"sham\".",
      call. = FALSE
    )
  }
  if (!reference %in% names(arms)) {
    stop("The reference arm '", reference, "' is not in `arms`.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Stop when the unit nesting is too shallow for the layout
#' @noRd
#' @keywords internal
.check_unit_depth <- function(unit, layout) {
  if (layout == "split_plot" && length(unit) < 2L) {
    stop(
      "A `split_plot` layout needs `unit` to name the block and the whole ",
      "plot, outer to inner.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Resolve the slope argument to a unit name or NULL
#' @noRd
#' @keywords internal
.resolve_slope <- function(slope, unit, layout) {
  if (is.null(slope)) {
    return(if (layout == "multi_site") unit[1L] else NULL)
  }
  if (isFALSE(slope)) {
    return(NULL)
  }
  if (!is.character(slope) || length(slope) != 1L || !slope %in% unit) {
    stop(
      "`slope` must be NULL, FALSE or one of the unit columns: ",
      paste(unit, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  slope
}

#' Number of innermost units carrying more than one treatment level
#' @noRd
#' @keywords internal
.n_units_varying <- function(d, treatment, unit) {
  key <- interaction(d[unit], drop = TRUE)
  n_lev <- tapply(
    d[[treatment]],
    key,
    function(z) length(unique(z[!is.na(z)]))
  )
  sum(n_lev > 1L)
}

#' Number of units (given by the `key` columns) lacking at least one arm level
#' @noRd
#' @keywords internal
.n_units_missing_arm <- function(d, treatment, key, levels) {
  k <- interaction(d[key], drop = TRUE)
  has_all <- tapply(
    as.character(d[[treatment]]),
    k,
    function(z) all(levels %in% z)
  )
  sum(!has_all)
}
