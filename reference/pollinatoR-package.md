# pollinatoR: uncertainty-aware pollinator-dependence analytics

Tools for estimating the contribution of insect pollinators to crop
yield from pollinator-exclusion experiments (open versus caged, and
related contrasts), and for relating that contribution to the pollinator
community. The methods are general to open-versus-exclusion designs and
not specific to any one crop or study.

## Details

The package is organised around a small set of verbs.
[`pollin_data()`](https://max578.github.io/pollinatoR/reference/pollin_data.md)
and
[`pollin_design()`](https://max578.github.io/pollinatoR/reference/pollin_design.md)
declare the study;
[`pollin_estimand()`](https://max578.github.io/pollinatoR/reference/pollin_estimand.md)
names the contrast and its scale;
[`pollin_fit()`](https://max578.github.io/pollinatoR/reference/pollin_fit.md)
and
[`pollin_effects()`](https://max578.github.io/pollinatoR/reference/pollin_effects.md)
estimate it;
[`pollin_contribution()`](https://max578.github.io/pollinatoR/reference/pollin_contribution.md)
estimates the per-outcome contribution as a model-based treatment effect
with uncertainty;
[`pollin_moderation()`](https://max578.github.io/pollinatoR/reference/pollin_moderation.md)
asks which exposures drive the contribution while guarding against
site-level confounders through a within-between decomposition;
[`pollin_index()`](https://max578.github.io/pollinatoR/reference/pollin_index.md)
and
[`pollin_rii()`](https://max578.github.io/pollinatoR/reference/pollin_rii.md)
give uncertainty-aware dependence indices;
[`pollin_guard_spurious()`](https://max578.github.io/pollinatoR/reference/pollin_guard_spurious.md)
flags the spurious correlation that arises when an index is regressed on
a yield axis sharing a term; and
[`pollin_validate()`](https://max578.github.io/pollinatoR/reference/pollin_validate.md)
runs structural integrity checks.

## See also

Useful links:

- <https://github.com/max578/pollinatoR>

- <https://max578.github.io/pollinatoR/>

- Report bugs at <https://github.com/max578/pollinatoR/issues>

## Author

**Maintainer**: Max Moldovan <max.moldovan@gmail.com>
([ORCID](https://orcid.org/0000-0001-9680-8474))

Authors:

- Max Moldovan <max.moldovan@gmail.com>
  ([ORCID](https://orcid.org/0000-0001-9680-8474))
