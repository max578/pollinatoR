#' pollinatoR: uncertainty-aware pollinator-dependence analytics
#'
#' Tools for estimating the contribution of insect pollinators to crop yield
#' from pollinator-exclusion experiments (open versus caged, and related
#' contrasts), and for relating that contribution to the pollinator community.
#' The methods are general to open-versus-exclusion designs and not specific to
#' any one crop or study.
#'
#' The package is organised around a small set of verbs. [pollin_spec()]
#' declares the study; [pollin_contribution()] estimates the per-outcome
#' contribution as a model-based treatment effect with uncertainty;
#' [pollin_moderation()] asks which exposures drive the contribution while
#' guarding against site-level confounders through a within-between
#' decomposition; [pollin_index()] and [pollin_rii()] give uncertainty-aware
#' dependence indices; [pollin_guard_spurious()] flags the spurious correlation
#' that arises when an index is regressed on a yield axis sharing a term; and
#' [pollin_validate()] runs structural integrity checks.
#'
#' @keywords internal
#' @importFrom stats lm coef sd quantile p.adjust as.formula reformulate
#'   predict terms median complete.cases
"_PACKAGE"
