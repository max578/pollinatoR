# Causal estimand (layer L2) --------------------------------------------------

#' Declare the estimand of a pollination experiment
#'
#' Names the contrast, the scale it is reported on and the basis of any share,
#' for a declared [pollin_design()]. Two estimand families are available:
#'
#' * `"deficit"` (E1): supplemented versus open, the yield gap left by the
#'   pollination the crop received.
#' * `"dependence"` (E2): excluded versus the design's reference arm (open, or
#'   supplemented), the yield owed to pollinators.
#'
#' Writing `high` for the expected response in the higher arm (open or
#' supplemented) and `low` for the lower one, the scales are:
#'
#' | Scale | Value |
#' |---|---|
#' | `"difference"` | `high - low` |
#' | `"ratio"`, `basis = "loss"` | `1 - low / high`, the share of `high` lost without the contrast |
#' | `"ratio"`, `basis = "gain"` | `high / low - 1`, the gain over `low` |
#' | `"rii"` | `(high - low) / (high + low)`, equal to `D / (2 - D)` for the loss share `D` |
#' | `"lnrr"` | `log(high / low)` |
#'
#' A loss share and a gain share answer different questions: an excluded arm
#' yielding 75% of the open arm is a 25% loss share and a 33% gain. `basis` has
#' no default; every analysis states it. It changes the value only on the ratio
#' scale and is recorded on every scale.
#'
#' @param design A [pollin_design()] object.
#' @param family `"deficit"` or `"dependence"`. There is no default.
#' @param scale One of `"difference"`, `"ratio"`, `"rii"` or `"lnrr"`. There
#'   is no default.
#' @param basis `"loss"` or `"gain"`. There is no default.
#' @param population The population the estimand refers to. Only `"trial"`
#'   (the units in the experiment) is available.
#'
#' @return An object of class `pollin_estimand`: a list holding `design`,
#'   `family`, `scale`, `basis`, `population`, the `contrast` (the treatment
#'   levels of the `high` and `low` arms) and the `assumptions` table.
#' @seealso [pollin_assumptions()] for the assumptions the estimand rests on.
#' @export
#' @examples
#' df <- data.frame(
#'   pods = rpois(48, 30),
#'   trt = rep(c("Cage", "Open"), 24),
#'   Site = rep(letters[1:4], each = 12)
#' )
#' pd <- pollin_data(df, "pods", "trt", unit = "Site")
#' ds <- pollin_design(pd, c(excluded = "Cage", open = "Open"), "multi_site")
#' pollin_estimand(ds, family = "dependence", scale = "ratio", basis = "loss")
pollin_estimand <- function(
  design,
  family,
  scale,
  basis,
  population = "trial"
) {
  .check_class(design, "pollin_design", arg = "design")
  .check_choice(family, "family", c("deficit", "dependence"))
  .check_choice(scale, "scale", c("difference", "ratio", "rii", "lnrr"))
  .check_choice(basis, "basis", c("loss", "gain"))
  .check_choice(population, "population", "trial")

  arms <- design$arms
  roles <- if (family == "deficit") {
    c(high = "supplemented", low = "open")
  } else {
    c(high = design$reference, low = "excluded")
  }
  absent <- setdiff(roles, names(arms))
  if (length(absent)) {
    stop(
      "The `",
      family,
      "` estimand needs the ",
      paste(absent, collapse = " and "),
      " arm; the design has ",
      paste(names(arms), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  structure(
    list(
      design = design,
      family = family,
      scale = scale,
      basis = basis,
      population = population,
      contrast = c(
        high = unname(arms[roles["high"]]),
        low = unname(arms[roles["low"]])
      ),
      assumptions = .assumptions(design$layout, family, scale, basis)
    ),
    class = "pollin_estimand"
  )
}

#' List the assumptions an estimand rests on
#'
#' Returns the identification assumptions for the estimand's layout, family
#' and scale. Each is a statement about the experiment that the data cannot
#' confirm on their own.
#'
#' @param x A [pollin_estimand()] object.
#'
#' @return A data frame with columns `assumption` and `statement`.
#' @export
#' @examples
#' df <- data.frame(
#'   pods = rpois(48, 30),
#'   trt = rep(c("Cage", "Open"), 24),
#'   Site = rep(letters[1:4], each = 12)
#' )
#' pd <- pollin_data(df, "pods", "trt", unit = "Site")
#' ds <- pollin_design(pd, c(excluded = "Cage", open = "Open"), "multi_site")
#' est <- pollin_estimand(ds, "dependence", "ratio", "loss")
#' pollin_assumptions(est)
pollin_assumptions <- function(x) {
  .check_class(x, "pollin_estimand")
  x$assumptions
}

#' @export
print.pollin_estimand <- function(x, ...) {
  label <- if (x$family == "deficit") "deficit (E1)" else "dependence (E2)"
  cat("<pollin_estimand>\n")
  cat("  family     :", label, "\n")
  cat("  contrast   :", x$contrast[["high"]], "vs", x$contrast[["low"]], "\n")
  cat("  scale      :", x$scale, "\n")
  cat("  basis      :", x$basis, "\n")
  cat("  population :", x$population, "\n")
  cat("  layout     :", x$design$layout, "\n")
  cat(
    "  assumptions:",
    nrow(x$assumptions),
    "; see `pollin_assumptions()`\n"
  )
  invisible(x)
}

# Internal ---------------------------------------------------------------------

#' Stop unless `value` is one string from `choices`
#' @noRd
#' @keywords internal
.check_choice <- function(value, arg, choices) {
  if (!is.character(value) || length(value) != 1L || !value %in% choices) {
    stop(
      "`",
      arg,
      "` must be one of: ",
      paste0("\"", choices, "\"", collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Assumption table for a layout, family and scale
#' @noRd
#' @keywords internal
.assumptions <- function(layout, family, scale, basis) {
  a <- c(
    consistency = "Each arm is one well-defined treatment, applied the same way on every unit.",
    excludability = if (family == "dependence") {
      "The exclusion device changes the response only by excluding pollinators (no cage effect on microclimate, pests or handling)."
    } else {
      "Hand supplementation changes the response only through the pollen added (no damage from handling)."
    },
    no_interference = "One unit's treatment does not change another unit's response, for example by displacing pollinators to neighbouring units."
  )
  a["assignment"] <- switch(
    layout,
    crd = "Arms were randomised to units.",
    rcbd = "Arms were randomised to units within each block.",
    split_plot = "Arms were randomised to sub-plots within each whole plot.",
    matched_pairs = "Arms were randomised within each pair, and the pair members are alike apart from the arm.",
    multi_site = "Arms were randomised within each site; the contrast is read within sites."
  )
  if (scale %in% c("ratio", "rii", "lnrr")) {
    a[
      "positive_means"
    ] <- "The expected response is positive in both arms, so the ratio is defined."
  }
  a[
    "population"
  ] <- "The estimand refers to the units in the trial; carrying it to other fields or seasons is a further assumption."
  data.frame(
    assumption = names(a),
    statement = unname(a),
    stringsAsFactors = FALSE
  )
}
