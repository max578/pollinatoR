# goodpractice dispositions, 0.2.0 candidate

Run 2026-11-09 on `stage2-value` with goodpractice 1.2.0, all checks except coverage and
`rcmdcheck` (run separately: coverage 98% with `NOT_CRAN=true`, every file at 90% or above;
`R CMD check --as-cran` 0 E / 0 W / 2 NOTEs, new submission with a development version and the
licence stub, both known).

The flags fall in the classes dispositioned for 0.1.0 in `goodpractice_S1.md`, with the same
reasons:

| Check | Where | Disposition |
|---|---|---|
| tidyverse_no_missing | `draws.R:68`, `effects.R:69`, `guards.R:49,360`, `index.R:47` and two more | Deliberate: an argument with no default (no silent choice) gets a tailored error. `draws.R` follows the same rule for `method` |
| line length | roxygen tables in `design.R`, `estimand.R` | Markdown tables in roxygen cannot wrap; code lines formatted by Air |
| tidyverse_r_file_names | `R/pollinatoR-package.R` | The roxygen convention for the package-level doc file |
| tidyverse_test_file_names | `compensation.R`, `contribution.R`, `index.R`, `market.R`, `moderation.R` and four more | Tests are organised by layer and purpose; `market.R` and `value.R` share `test-value.R` |
| tidyverse_export_order | `moderation.R:4` | Pre-existing layout, golden-covered, unchanged |

Closed since 0.1.0: `roxygen2_duplicate_params` (the `@inheritParams` pass), the `.need()` dead
code and the `spec.R` export-order flags (file removed with `pollin_spec()`).
