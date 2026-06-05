# pollinatoR 0.0.0.9000

* Initial development scaffold.
* `pollin_spec()` declares an exclusion-experiment study (outcomes, treatment,
  nested random structure, exposures, adjustment, latent grouping).
* `pollin_contribution()` estimates the per-outcome pollinator contribution
  (open minus exclusion) as a model-based average treatment effect with
  propagated uncertainty, carrying a random treatment slope on the latent unit
  to absorb unobserved between-unit (e.g. genotype) heterogeneity.
* `pollin_moderation()` estimates how the contribution varies with exposures
  (treatment-by-exposure conditional effects) using a within-between
  decomposition so the within-unit estimate is robust to site-level confounders.
* `pollin_index()` / `pollin_rii()` provide uncertainty-aware dependence indices,
  including the Relative Interaction Index with a bootstrap interval.
* `pollin_guard_spurious()` tests an index-versus-yield slope against a
  permutation null that isolates the mechanical (shared-term) coupling.
* `pollin_validate()` runs structural and internal-consistency integrity checks.
* `pollin_compensation()` summarises yield-component trade-offs.
