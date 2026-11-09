# make_golden_v0.2.0.R -- Capture the draws, population, value and total
# outputs of 0.2.0.
#
# Run once, from the package root, on the 0.2.0 code (before any Stage 3
# change under R/):
#   Rscript --vanilla validation/golden/make_golden_v0.2.0.R
# Writes tests/testthat/_golden/golden_v0.2.0.rds. The file is the oracle for
# Stage 3 (criterion C1) and is never regenerated from later code.

pkgload::load_all(".", quiet = TRUE)
source("tests/testthat/helper-golden.R")

out <- golden_capture_v020()

attr(out, "provenance") <- list(
  package_version = as.character(utils::packageVersion("pollinatoR")),
  git_head = system("git rev-parse HEAD", intern = TRUE),
  glmmTMB = as.character(utils::packageVersion("glmmTMB")),
  TMB = as.character(utils::packageVersion("TMB")),
  marginaleffects = as.character(utils::packageVersion("marginaleffects")),
  R = R.version.string,
  created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
)

dest <- file.path("tests", "testthat", "_golden", "golden_v0.2.0.rds")
if (file.exists(dest)) {
  stop("`", dest, "` exists; the golden file is never overwritten.", call. = FALSE)
}
saveRDS(out, dest)
message("Wrote ", dest)
