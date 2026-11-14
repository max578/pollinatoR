# make_golden_v0.5.0_decompose.R -- Capture the 0.5.0 whole-yield
# decomposition of the test fixture, so later releases cannot move it.
#
# Run once, from the package root:
#   Rscript --vanilla validation/golden/make_golden_v0.5.0_decompose.R
# Writes tests/testthat/_golden/golden_v0.5.0_decompose.rds, never overwritten.

pkgload::load_all(".", quiet = TRUE)
source("tests/testthat/helper-decompose.R")
source("tests/testthat/helper-golden.R")

out <- golden_capture_v050_decompose()

attr(out, "provenance") <- list(
  package_version = as.character(utils::packageVersion("pollinatoR")),
  git_head = system("git rev-parse HEAD", intern = TRUE),
  glmmTMB = as.character(utils::packageVersion("glmmTMB")),
  TMB = as.character(utils::packageVersion("TMB")),
  Matrix = as.character(utils::packageVersion("Matrix")),
  R = R.version.string,
  created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
)

dest <- file.path("tests", "testthat", "_golden", "golden_v0.5.0_decompose.rds")
if (file.exists(dest)) {
  stop("`", dest, "` exists; the golden file is never overwritten.", call. = FALSE)
}
saveRDS(out, dest)
message("Wrote ", dest)
