# pollinatoR (development version)

* New `pollin_data()` declares the unit-of-record table: response, treatment,
  unit nesting, and optional grouping, denominator, exposure and adjustment
  columns. Roles are checked against the data and the structural checks are
  stored with the object.
* New `pollin_design()` declares the arms (excluded, open, supplemented), the
  layout (`crd`, `rcbd`, `split_plot`, `matched_pairs`, `multi_site`) and the
  reference arm. The layout sets the random structure; a random treatment
  slope sits on the site in `multi_site` and is opt-in elsewhere. Units lacking
  an arm are counted in a warning.
* `pollin_validate()` now takes a `pollin_data` object. Its first argument is
  renamed from `spec` to `x`; a `pollin_spec` is still accepted.
* Removed six unused packages from Suggests (DHARMa, emmeans, brms, boot, fst,
  data.table).

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
