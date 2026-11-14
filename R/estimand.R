# Causal estimand (layer L2) --------------------------------------------------

#' Declare the estimand of a pollination experiment
#'
#' Names the contrast, the scale it is reported on and the basis of any share,
#' for a declared [pollin_design()]. Three estimand families are available:
#'
#' * `"deficit"` (E1): supplemented versus open, the yield gap left by the
#'   pollination the crop received.
#' * `"dependence"` (E2): excluded versus the design's reference arm (open,
#'   supplemented or sham), the yield owed to pollinators. Against a sham arm
#'   it is net of the cage effect.
#' * `"cage"`: sham versus open, the effect of the exclusion device itself on
#'   plants that pollinators still reach. It has no basis, and its ratio scale
#'   is the cage factor `sham / open`; scales `"difference"`, `"ratio"` and
#'   `"lnrr"`.
#'
#' Every estimand is a total effect unless `effect = "direct"` is asked for.
#' A direct effect holds the named `mediators` fixed (for example vegetative
#' weight measured at harvest, which pollination itself can change); it is
#' never the default, and with a `dag` it is refused when the graph leaves a
#' path from the mediators to the outcome open, as an unmeasured common cause
#' of the mediator and the yield does.
#'
#' Writing `high` for the expected response in the higher arm (open or
#' supplemented or sham) and `low` for the lower one, the scales are:
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
#' @param family `"deficit"`, `"dependence"` or `"cage"`. There is no
#'   default.
#' @param scale One of `"difference"`, `"ratio"`, `"rii"` or `"lnrr"`. There
#'   is no default.
#' @param basis `"loss"` or `"gain"`. There is no default, except for the
#'   `"cage"` family, which takes none (leave it out or pass `NULL`).
#' @param population The population the estimand refers to: `"trial"` (the
#'   default, the units in the experiment) or a [pollin_population()] object,
#'   strata of one column weighted by area. A population estimand is fitted
#'   with the strata as `by` and read with [pollin_draws()].
#' @param effect `"total"` (the default) or `"direct"`.
#' @param mediators For `effect = "direct"`, the declared `adjust` columns held
#'   fixed. Required for a direct effect, refused for a total one.
#' @param dag Optional [pollin_dag()] object built from `design`. With a graph,
#'   a total effect requires the unadjusted contrast to be valid, and a direct
#'   effect requires each mediator to lie on a path from the treatment to the
#'   outcome and the declared `adjust` columns to contain a direct-effect
#'   adjustment set of the graph.
#'
#' @return An object of class `pollin_estimand`: a list holding `design`,
#'   `family`, `scale`, `basis` (`NA` for the cage family), `population`,
#'   `effect`, `mediators`, `adjustment` (the columns a direct effect is
#'   adjusted for), `dag`, the `contrast` (the treatment levels of the `high`
#'   and `low` arms) and the `assumptions` table.
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
  population = "trial",
  effect = "total",
  mediators = NULL,
  dag = NULL
) {
  .check_class(design, "pollin_design", arg = "design")
  .check_choice(family, "family", c("deficit", "dependence", "cage"))
  if (family == "cage") {
    .check_choice(scale, "scale", c("difference", "ratio", "lnrr"))
    if (!missing(basis) && !is.null(basis)) {
      stop(
        "The `cage` family has no basis: its ratio is the cage factor ",
        "sham / open. Leave `basis` out.",
        call. = FALSE
      )
    }
    basis <- NA_character_
  } else {
    .check_choice(scale, "scale", c("difference", "ratio", "rii", "lnrr"))
    .check_choice(basis, "basis", c("loss", "gain"))
  }
  .check_population(population, design$data$data)
  .check_choice(effect, "effect", c("total", "direct"))
  adjustment <- .check_effect(design, family, effect, mediators, dag)

  arms <- design$arms
  roles <- switch(
    family,
    deficit = c(high = "supplemented", low = "open"),
    dependence = c(high = design$reference, low = "excluded"),
    cage = c(high = "sham", low = "open")
  )
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
      effect = effect,
      mediators = mediators,
      adjustment = adjustment,
      dag = dag,
      contrast = c(
        high = unname(arms[roles["high"]]),
        low = unname(arms[roles["low"]])
      ),
      assumptions = .assumptions(
        design$layout,
        family,
        scale,
        basis,
        population,
        reference = design$reference,
        sham = "sham" %in% names(arms),
        effect = effect,
        mediators = mediators,
        adjustment = adjustment,
        dag = dag
      )
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
  cat("<pollin_estimand>\n")
  cat("  family     :", .family_label(x$family), "\n")
  cat("  contrast   :", x$contrast[["high"]], "vs", x$contrast[["low"]], "\n")
  cat("  effect     :", .effect_label(x$effect, x$mediators), "\n")
  cat("  scale      :", x$scale, "\n")
  cat("  basis      :", .basis_label(x$basis), "\n")
  cat("  population :", .population_label(x$population), "\n")
  cat("  layout     :", x$design$layout, "\n")
  cat(
    "  assumptions:",
    nrow(x$assumptions),
    "; see `pollin_assumptions()`\n"
  )
  invisible(x)
}

# Internal ---------------------------------------------------------------------

#' The family in words, for printing
#' @noRd
#' @keywords internal
.family_label <- function(family) {
  switch(
    family,
    deficit = "deficit (E1)",
    dependence = "dependence (E2)",
    cage = "cage effect (sham versus open)"
  )
}

#' The basis in words, for printing; the cage family has none
#' @noRd
#' @keywords internal
.basis_label <- function(basis) {
  if (is.na(basis)) "none (cage factor)" else basis
}

#' The effect type in words, for printing
#' @noRd
#' @keywords internal
.effect_label <- function(effect, mediators = NULL) {
  if (identical(effect, "direct")) {
    return(paste0(
      "direct, holding ",
      paste(mediators, collapse = ", "),
      " fixed"
    ))
  }
  "total"
}

#' Check the effect type against the mediators and the graph; returns the
#' columns a direct effect is adjusted for (`NULL` for a total effect)
#' @noRd
#' @keywords internal
.check_effect <- function(design, family, effect, mediators, dag) {
  roles <- design$data$roles
  if (!is.null(dag)) {
    .check_class(dag, "pollin_dag", arg = "dag")
    if (!identical(dag$roles, roles)) {
      stop(
        "`dag` was built from a different design; rebuild it with ",
        "`pollin_dag(design, ...)`.",
        call. = FALSE
      )
    }
  }
  if (effect == "total") {
    if (!is.null(mediators)) {
      stop(
        "`mediators` are named only for a direct effect; adjusting for a ",
        "mediator under the total effect would report neither. Use ",
        "`effect = \"direct\"`.",
        call. = FALSE
      )
    }
    if (
      !is.null(dag) &&
        !dagitty::isAdjustmentSet(
          dag$graph,
          character(0),
          exposure = dag$treatment,
          outcome = dag$outcome
        )
    ) {
      stop(
        "The graph gives the treatment a common cause with the outcome, so ",
        "the unadjusted contrast is biased. Arms are randomised; check the ",
        "edges into `",
        dag$treatment,
        "`.",
        call. = FALSE
      )
    }
    return(NULL)
  }
  if (family == "cage") {
    stop("The cage effect is reported as a total effect.", call. = FALSE)
  }
  if (!is.character(mediators) || length(mediators) == 0L) {
    stop(
      "A direct effect needs `mediators`: the declared `adjust` columns held ",
      "fixed.",
      call. = FALSE
    )
  }
  absent <- setdiff(mediators, roles$adjust)
  if (length(absent)) {
    stop(
      "Mediators must be declared `adjust` columns in `pollin_data()`: ",
      paste(absent, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (is.null(dag)) {
    return(roles$adjust)
  }
  cls <- dag$classes[dag$classes$target == dag$treatment, , drop = FALSE]
  not_med <- setdiff(mediators, cls$node[cls$class == "mediator"])
  if (length(not_med)) {
    stop(
      "Not on a path from the treatment to the outcome in the graph: ",
      paste(not_med, collapse = ", "),
      ". A direct effect holds mediators fixed.",
      call. = FALSE
    )
  }
  sets <- dag$sets$direct
  if (length(sets) == 0L) {
    g <- dag$graph
    culprit <- intersect(
      dagitty::latents(g),
      intersect(
        unlist(lapply(mediators, function(m) dagitty::ancestors(g, m))),
        dagitty::ancestors(g, dag$outcome)
      )
    )
    stop(
      "The graph identifies no direct effect: holding ",
      paste(mediators, collapse = ", "),
      " fixed opens a path to the outcome that no measured column closes",
      if (length(culprit)) {
        paste0(" (through the latent ", paste(culprit, collapse = ", "), ")")
      },
      ". Report the total effect.",
      call. = FALSE
    )
  }
  fits <- Filter(function(z) all(z %in% roles$adjust), sets)
  if (length(fits) == 0L) {
    stop(
      "The declared `adjust` columns contain no direct-effect set of the ",
      "graph; declare one of: ",
      paste(vapply(sets, .dag_set_label, character(1)), collapse = " or "),
      ".",
      call. = FALSE
    )
  }
  fits[[which.min(lengths(fits))]]
}

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

#' Stop unless `population` is "trial" or a `pollin_population` whose strata
#' are the levels of its column in the data
#' @noRd
#' @keywords internal
.check_population <- function(population, data) {
  if (inherits(population, "pollin_population")) {
    return(.population_check_data(population, data))
  }
  if (!identical(population, "trial")) {
    stop(
      "`population` must be \"trial\" or a `pollin_population` object.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' The population in words, for printing
#' @noRd
#' @keywords internal
.population_label <- function(population) {
  if (!inherits(population, "pollin_population")) {
    return(population)
  }
  paste0(
    if (!is.null(population$label)) paste0(population$label, ": "),
    length(population$weights),
    if (length(population$weights) == 1L) " stratum" else " strata",
    " of ",
    population$by,
    ", weighted by area"
  )
}

#' Assumption table for a layout, family and scale
#' @noRd
#' @keywords internal
.assumptions <- function(
  layout,
  family,
  scale,
  basis,
  population = "trial",
  reference = "open",
  sham = FALSE,
  effect = "total",
  mediators = NULL,
  adjustment = NULL,
  dag = NULL
) {
  a <- c(
    consistency = "Each arm is one well-defined treatment, applied the same way on every unit.",
    excludability = if (family == "cage") {
      "The sham device matches the exclusion device in every effect but excluding pollinators."
    } else if (family == "dependence" && reference == "sham") {
      "The sham device matches the exclusion device in every effect but excluding pollinators, and the cage changes the response by the same factor with and without pollinators; the contrast is then net of the cage effect."
    } else if (family == "dependence" && sham) {
      "The exclusion device changes the response only by excluding pollinators. The design has a sham arm, so this can be checked: the cage effect is sham versus open, and this contrast includes it."
    } else if (family == "dependence") {
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
  a["population"] <- if (inherits(population, "pollin_population")) {
    .population_statement(population)
  } else {
    "The estimand refers to the units in the trial; carrying it to other fields or seasons is a further assumption."
  }
  a["effect"] <- .effect_statement(effect, mediators, dag)
  if (!is.null(dag) || effect == "direct") {
    a["adjustment"] <- .adjustment_statement(effect, adjustment, dag)
  }
  data.frame(
    assumption = names(a),
    statement = unname(a),
    stringsAsFactors = FALSE
  )
}

#' The effect-type assumption
#' @noRd
#' @keywords internal
.effect_statement <- function(effect, mediators, dag) {
  if (effect == "total") {
    return(paste0(
      "Total effect: the contrast carries every path from the arm to the ",
      "response, including any through vegetative growth or other ",
      "intermediates."
    ))
  }
  m <- paste(mediators, collapse = ", ")
  paste0(
    "Direct effect holding ",
    m,
    " fixed: the arm does not change how ",
    m,
    " acts on the response (no arm-by-mediator interaction)",
    if (is.null(dag)) {
      paste0(
        ", and no unmeasured variable causes both ",
        m,
        " and the response. No graph was given, so neither is checked."
      )
    } else {
      "; the graph closes every other path through the mediators."
    }
  )
}

#' The adjustment assumption: the set used and what it closes
#' @noRd
#' @keywords internal
.adjustment_statement <- function(effect, adjustment, dag) {
  if (effect == "total") {
    return(paste0(
      "The graph implies the unadjusted contrast is unbiased: no common cause ",
      "of the arm and the outcome, as randomisation provides."
    ))
  }
  set <- .dag_set_label(adjustment)
  if (is.null(dag)) {
    return(paste0(
      "Adjusted for ",
      set,
      ", the declared adjust columns; without a graph the set is the ",
      "analyst's declaration."
    ))
  }
  paste0(
    "Adjusted for ",
    set,
    ", a direct-effect set of the graph: it blocks every path from the arm ",
    "to the outcome except the direct one."
  )
}
