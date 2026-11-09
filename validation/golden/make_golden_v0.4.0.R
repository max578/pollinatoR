# make_golden_v0.4.0.R -- Capture the outputs that 0.4.0 changes: the joint
# delta effects on the 0.1.0 inputs, and the draws, population, value and
# total outputs of 0.2.0 now drawn from the joint covariance.
#
# Run once, from the package root, on the 0.4.0 code after the joint-interval
# tests pass:
#   Rscript --vanilla validation/golden/make_golden_v0.4.0.R
# Writes tests/testthat/_golden/golden_v0.4.0.rds, never overwritten.

pkgload::load_all(".", quiet = TRUE)
source("tests/testthat/helper-golden.R")

out <- suppressWarnings(golden_capture_v040())

attr(out, "provenance") <- list(
  package_version = as.character(utils::packageVersion("pollinatoR")),
  git_head = system("git rev-parse HEAD", intern = TRUE),
  glmmTMB = as.character(utils::packageVersion("glmmTMB")),
  TMB = as.character(utils::packageVersion("TMB")),
  Matrix = as.character(utils::packageVersion("Matrix")),
  R = R.version.string,
  created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
)

dest <- file.path("tests", "testthat", "_golden", "golden_v0.4.0.rds")
if (file.exists(dest)) {
  stop("`", dest, "` exists; the golden file is never overwritten.", call. = FALSE)
}
saveRDS(out, dest)
message("Wrote ", dest)
