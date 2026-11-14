# make_golden_v0.4.0_cate.R -- Capture the 0.4.0 CATE outputs that
# Stage 4 (0.5.0) must leave unchanged: the shipped example on the ratio and
# difference scales.
#
# Run once, from the package root of a checkout of v0.4.0 (2a5748b) with this
# script and tests/testthat/helper-golden.R copied in (S4 Step 5):
#   Rscript --vanilla validation/golden/make_golden_v0.4.0_cate.R
# Writes tests/testthat/_golden/golden_v0.4.0_cate.rds, never overwritten.

pkgload::load_all(".", quiet = TRUE)
source("tests/testthat/helper-golden.R")

out <- suppressWarnings(suppressMessages(golden_capture_v040_cate()))

attr(out, "provenance") <- list(
  package_version = as.character(utils::packageVersion("pollinatoR")),
  git_head = system("git rev-parse HEAD", intern = TRUE),
  glmmTMB = as.character(utils::packageVersion("glmmTMB")),
  TMB = as.character(utils::packageVersion("TMB")),
  Matrix = as.character(utils::packageVersion("Matrix")),
  R = R.version.string,
  created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
)

dest <- file.path("tests", "testthat", "_golden", "golden_v0.4.0_cate.rds")
if (file.exists(dest)) {
  stop("`", dest, "` exists; the golden file is never overwritten.", call. = FALSE)
}
saveRDS(out, dest)
message("Wrote ", dest)
