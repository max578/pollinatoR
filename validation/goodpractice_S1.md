# goodpractice at 0.1.0 candidate (de2872f, 2026-11-08)

Run: `goodpractice::gp(".")` with `NOT_CRAN=true`, all checks except `covr`/`rcmdcheck`
(run separately: coverage 96.9 %; `R CMD check --as-cran` 0 E / 0 W / 2 NOTEs).
URL checks skipped (no network in the run). Each failed check, with its disposition:

| Check | Where | Disposition |
|---|---|---|
| complexity_unused_internal | `R/utils.R:67` `.need()` | Dead code predating Stage 1 (its callers left with the six dropped Suggests). Flagged for Max; not deleted. |
| complexity_function_length, cyclocomp | `pollin_guard_ecological` (26), `pollin_effects`, `pollin_design`, the spec path of `pollin_contribution` | Accepted for 0.1.0; the bodies are sectioned with banners. Candidates for splitting at 0.2.0. |
| tidyverse_no_missing | `effects.R:74`, `index.R:57,156`, `guards.R:57,363` | Deliberate: arguments with no default (no silent choice) get a tailored error via `missing()`. |
| lintr_assignment_linter / undesirable_operator | `validate.R:32` `<<-` | Pre-existing, golden-covered; a closure accumulator. Kept to keep the golden check table untouched. |
| line length (lintr, tidyverse) | mostly roxygen tables in `design.R`, `estimand.R` | Markdown tables in roxygen cannot wrap; code lines formatted by Air. |
| lintr_paste_linter | `cat(..., )` / `paste()` calls | Not inspected line by line (27 lines); style only, no behaviour change. Review at 0.2.0. |
| lintr_fixed_regex_linter, nested_ifelse | `moderation.R:102-114` | Pre-existing moderation code covered by the golden test; left unchanged at 0.1.0. |
| lintr_unnecessary_lambda | `utils.R:25` | Pre-existing; `na.rm = TRUE` needs the wrapper. False positive. |
| roxygen2_duplicate_params | `contribution.R`, `fit.R`, `moderation.R` | Accepted; `@inheritParams` consolidation at 0.2.0. |
| tidyverse_r_file_names | `R/pollinatoR-package.R` | The roxygen/usethis convention for the package-level doc file. |
| tidyverse_test_file_names | several | Tests are organised by layer (`test-data.R`, `test-design.R`, ...) and by purpose (`test-golden.R`, `test-ports.R`). |
| tidyverse_export_order | `moderation.R:4`, `spec.R:112,127` | Pre-existing layout; left unchanged at 0.1.0. |

Side effect of the run: goodpractice installed `RTMB` from CRAN into the user library
while preparing its dependency checks.
