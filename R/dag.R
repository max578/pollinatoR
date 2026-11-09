# Causal graph (layer L2) ------------------------------------------------------

#' Attach a causal graph to a design
#'
#' Reads a directed acyclic graph drawn in `dagitty` syntax against the columns
#' of a [pollin_design()] and derives what the graph implies for the analysis:
#' the adjustment sets for the treatment contrast (total and direct effects),
#' the adjustment set for each declared exposure, and the class of every
#' observed node. Pass the result to [pollin_estimand()] as `dag`.
#'
#' Nodes are column names. The treatment column must be a node, and so must
#' every declared `adjust` and `exposures` column. `outcome` names the node that
#' stands for every response column. Any other node that is not a column of the
#' data must be declared latent (`fertility [latent]`); a latent node never
#' enters an adjustment set. Draw an unmeasured common cause as a latent node,
#' not as a bidirected edge.
#'
#' Each observed node is classed against a target (the treatment, then each
#' exposure), with `S` the first minimal adjustment set for the target's total
#' effect:
#'
#' | Class | Meaning |
#' |---|---|
#' | `"mediator"` | a descendant of the target and an ancestor of the outcome |
#' | `"descendant"` | a descendant of the target, not on a path to the outcome |
#' | `"confounder"` | in at least one minimal adjustment set |
#' | `"collider"` | adding it to `S` makes the set invalid |
#' | `"instrument"` | adding it to `S` keeps the set valid; it is associated with the target and, given `S`, independent of the outcome once the target's own effects are removed |
#' | `"precision"` | adding it to `S` keeps the set valid; it is associated with the outcome given `S` and the target |
#' | `"neutral"` | none of the above |
#' | `"pre-treatment"` | not a descendant of the target, when no valid set exists |
#'
#' A graph cannot draw effect modification. The set returned for an exposure
#' is the set for the exposure's effect on the outcome; [pollin_moderation()]
#' borrows it for the exposure-by-treatment term, which is valid only if no
#' unmeasured variable both changes the treatment contrast and is associated
#' with the exposure.
#'
#' @param x A [pollin_design()] or [pollin_estimand()] object.
#' @param dag A `dagitty` object or a string in `dagitty` syntax.
#' @param outcome The node that stands for the response columns.
#'
#' @return An object of class `pollin_dag`: a list holding the `graph`, the
#'   `outcome` node, the `treatment` column, the data `roles`, the `sets` for
#'   the treatment (`total` and `direct`, each a list of character vectors; an
#'   empty list means no valid set exists), the `exposure_sets` and the
#'   `classes` table (`target`, `node`, `class`).
#' @seealso [pollin_dag_check()] to test the graph's implied independencies.
#' @export
#' @examplesIf requireNamespace("dagitty", quietly = TRUE)
#' df <- data.frame(
#'   pods = rpois(48, 30),
#'   trt = rep(c("Cage", "Open"), 24),
#'   Site = rep(letters[1:4], each = 12),
#'   veg = rnorm(48, 50, 5)
#' )
#' pd <- pollin_data(df, "pods", "trt", unit = "Site", adjust = "veg")
#' ds <- pollin_design(pd, c(excluded = "Cage", open = "Open"), "multi_site")
#' g <- pollin_dag(ds, "dag { trt -> veg -> yield; trt -> yield }", "yield")
#' g
pollin_dag <- function(x, dag, outcome) {
  if (inherits(x, "pollin_estimand")) {
    x <- x$design
  }
  .check_class(x, "pollin_design")
  if (!requireNamespace("dagitty", quietly = TRUE)) {
    stop(
      "`pollin_dag()` needs dagitty; install it with ",
      "install.packages(\"dagitty\").",
      call. = FALSE
    )
  }
  g <- .dag_parse(dag)
  roles <- x$data$roles
  trt <- roles$treatment
  .dag_check_nodes(g, outcome, roles, names(x$data$data))

  sets <- list(
    total = .dag_sets(g, trt, outcome, "total"),
    direct = .dag_sets(g, trt, outcome, "direct")
  )
  exposure_sets <- lapply(
    stats::setNames(roles$exposures, roles$exposures),
    function(ex) .dag_sets(g, ex, outcome, "total")
  )
  targets <- c(trt, roles$exposures)
  classes <- do.call(
    rbind,
    lapply(targets, function(t) .dag_classes(g, t, outcome))
  )

  structure(
    list(
      graph = g,
      outcome = outcome,
      treatment = trt,
      roles = roles,
      sets = sets,
      exposure_sets = exposure_sets,
      classes = classes
    ),
    class = "pollin_dag"
  )
}

#' @export
print.pollin_dag <- function(x, ...) {
  fmt <- function(s) {
    if (length(s) == 0L) {
      return("none (not identified)")
    }
    paste(vapply(s, .dag_set_label, character(1)), collapse = " or ")
  }
  cat("<pollin_dag>  ", x$treatment, " -> ", x$outcome, "\n", sep = "")
  lat <- dagitty::latents(x$graph)
  if (length(lat)) {
    cat("  latent        :", paste(lat, collapse = ", "), "\n")
  }
  cat("  total effect  :", fmt(x$sets$total), "\n")
  cat("  direct effect :", fmt(x$sets$direct), "\n")
  for (ex in names(x$exposure_sets)) {
    cat("  exposure", ex, ":", fmt(x$exposure_sets[[ex]]), "\n")
  }
  print(x$classes, row.names = FALSE)
  invisible(x)
}

#' Test the independencies a causal graph implies
#'
#' Tests each conditional independence the graph implies among the columns of
#' `data`, with [dagitty::localTests()] (`type = "cis"`, linear). Pass a table
#' at the level the covariates vary: a site table for site-level covariates.
#' Nodes that are not columns of `data` are treated as unobserved. The result
#' is reported, never gated: with few units a test detects only a strong
#' violation, and a graph that passes is not thereby confirmed.
#'
#' @param dag A [pollin_dag()] object.
#' @param data A data frame whose numeric columns are nodes of the graph.
#' @param conf_level Confidence level of the intervals.
#'
#' @return A data frame with columns `test`, `estimate`, `conf.low`,
#'   `conf.high`, `p.value` and `n` (the complete rows used).
#' @export
#' @examplesIf requireNamespace("dagitty", quietly = TRUE)
#' df <- data.frame(
#'   pods = rpois(48, 30),
#'   trt = rep(c("Cage", "Open"), 24),
#'   Site = rep(letters[1:4], each = 12),
#'   rain = rep(rnorm(4, 300, 40), each = 12),
#'   forest = rep(runif(4), each = 12)
#' )
#' pd <- pollin_data(df, "pods", "trt", unit = "Site", adjust = "rain")
#' ds <- pollin_design(pd, c(excluded = "Cage", open = "Open"), "multi_site")
#' g <- pollin_dag(ds, "dag { trt -> yield; rain -> yield; forest -> yield }",
#'                 "yield")
#' site <- unique(df[c("Site", "rain", "forest")])
#' pollin_dag_check(g, site)
pollin_dag_check <- function(dag, data, conf_level = 0.95) {
  .check_class(dag, "pollin_dag", arg = "dag")
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  g <- dag$graph
  nodes <- names(g)
  obs <- intersect(nodes, names(data))
  bad <- obs[!vapply(data[obs], is.numeric, logical(1))]
  if (length(bad)) {
    stop(
      "Graph nodes in `data` must be numeric: ",
      paste(bad, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  dagitty::latents(g) <- setdiff(nodes, obs)
  cis <- dagitty::impliedConditionalIndependencies(g)
  if (length(cis) == 0L) {
    stop(
      "The graph implies no testable independence among the columns of ",
      "`data`.",
      call. = FALSE
    )
  }
  d <- data[stats::complete.cases(data[obs]), obs, drop = FALSE]
  lt <- dagitty::localTests(
    g,
    data = d,
    type = "cis",
    tests = cis,
    conf.level = conf_level,
    abbreviate.names = FALSE
  )
  lo <- grep("%$", names(lt), value = TRUE)
  data.frame(
    test = rownames(lt),
    estimate = lt$estimate,
    conf.low = lt[[lo[1L]]],
    conf.high = lt[[lo[2L]]],
    p.value = lt$p.value,
    n = nrow(d),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

# Internal ---------------------------------------------------------------------

#' A `dagitty` DAG from an object or a string; stops on any other graph type
#' or on an edge that is not directed
#' @noRd
#' @keywords internal
.dag_parse <- function(dag) {
  g <- if (inherits(dag, "dagitty")) {
    dag
  } else if (is.character(dag) && length(dag) == 1L) {
    dagitty::dagitty(dag)
  } else {
    stop(
      "`dag` must be a dagitty object or a string in dagitty syntax.",
      call. = FALSE
    )
  }
  if (dagitty::graphType(g) != "dag") {
    stop(
      "`dag` must be a directed acyclic graph (`dag { ... }`).",
      call. = FALSE
    )
  }
  e <- dagitty::edges(g)
  if (any(e$e != "->")) {
    stop(
      "Every edge must be directed (`->`); draw an unmeasured common cause ",
      "as a latent node (`u [latent]; u -> a; u -> b`).",
      call. = FALSE
    )
  }
  g
}

#' Stop unless the graph names the treatment, the outcome and every declared
#' adjust and exposure column, and every other non-column node is latent
#' @noRd
#' @keywords internal
.dag_check_nodes <- function(g, outcome, roles, columns) {
  nodes <- names(g)
  if (!is.character(outcome) || length(outcome) != 1L) {
    stop("`outcome` must be one node name.", call. = FALSE)
  }
  need <- c(roles$treatment, outcome, roles$adjust, roles$exposures)
  missing <- setdiff(need, nodes)
  if (length(missing)) {
    stop(
      "Nodes missing from the graph: ",
      paste(missing, collapse = ", "),
      ". The graph must name the treatment, the outcome and every declared ",
      "adjust and exposure column.",
      call. = FALSE
    )
  }
  if (outcome == roles$treatment) {
    stop("`outcome` cannot be the treatment column.", call. = FALSE)
  }
  stray <- setdiff(nodes, c(columns, outcome, dagitty::latents(g)))
  if (length(stray)) {
    stop(
      "Nodes that are neither columns of the data nor latent: ",
      paste(stray, collapse = ", "),
      ". Declare an unmeasured node latent (`",
      stray[1L],
      " [latent]`).",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Minimal adjustment sets as a list of sorted character vectors; an empty
#' list means no valid set, `list(character(0))` means no adjustment is needed
#' @noRd
#' @keywords internal
.dag_sets <- function(g, exposure, outcome, effect) {
  s <- dagitty::adjustmentSets(
    g,
    exposure = exposure,
    outcome = outcome,
    effect = effect
  )
  lapply(unname(unclass(s)), function(z) sort(as.character(unlist(z))))
}

#' A set in braces, for printing
#' @noRd
#' @keywords internal
.dag_set_label <- function(s) {
  paste0("{", paste(s, collapse = ", "), "}")
}

#' The graph with the target's outgoing edges removed
#' @noRd
#' @keywords internal
.dag_cut_out <- function(g, target) {
  e <- dagitty::edges(g)
  e <- e[e$v != target, , drop = FALSE]
  lat <- dagitty::latents(g)
  body <- c(
    paste0(names(g), ifelse(names(g) %in% lat, " [latent]", "")),
    paste(e$v, "->", e$w)
  )
  dagitty::dagitty(paste0("dag { ", paste(body, collapse = "; "), " }"))
}

#' Class of every observed node other than the target and the outcome
#' @noRd
#' @keywords internal
.dag_classes <- function(g, target, outcome) {
  nodes <- setdiff(names(g), c(target, outcome, dagitty::latents(g)))
  if (length(nodes) == 0L) {
    return(data.frame(
      target = character(0),
      node = character(0),
      class = character(0),
      stringsAsFactors = FALSE
    ))
  }
  desc <- setdiff(dagitty::descendants(g, target), target)
  anc_y <- dagitty::ancestors(g, outcome)
  sets <- .dag_sets(g, target, outcome, "total")
  in_any <- unique(unlist(sets))
  g_cut <- .dag_cut_out(g, target)
  s <- if (length(sets)) sets[[1L]] else character(0)
  valid_with <- function(v) {
    dagitty::isAdjustmentSet(g, c(s, v), exposure = target, outcome = outcome)
  }
  cls <- vapply(
    nodes,
    function(v) {
      if (v %in% desc) {
        return(if (v %in% anc_y) "mediator" else "descendant")
      }
      if (length(sets) == 0L) {
        return("pre-treatment")
      }
      if (v %in% in_any) {
        return("confounder")
      }
      if (!valid_with(v)) {
        return("collider")
      }
      if (
        !dagitty::dseparated(g, v, target, s) &&
          dagitty::dseparated(g_cut, v, outcome, s)
      ) {
        return("instrument")
      }
      if (!dagitty::dseparated(g, v, outcome, c(s, target))) {
        return("precision")
      }
      "neutral"
    },
    character(1)
  )
  data.frame(
    target = target,
    node = nodes,
    class = unname(cls),
    stringsAsFactors = FALSE
  )
}
