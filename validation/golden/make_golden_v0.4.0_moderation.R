# make_golden_v0.4.0_moderation.R -- Capture the 0.4.0 moderation outputs that
# Stage 4 (0.5.0) must leave unchanged: the shipped example and the three-arm
# block trial, both responses. Only cleanly converged fits are captured.
#
# Run once, from the package root, on the 0.4.0 code (S4 Step 1):
#   Rscript --vanilla validation/golden/make_golden_v0.4.0_moderation.R
# Writes tests/testthat/_golden/golden_v0.4.0_moderation.rds, never overwritten.

pkgload::load_all(".", quiet = TRUE)
source("tests/testthat/helper-golden.R")

out <- suppressWarnings(suppressMessages(golden_capture_v040_moderation()))

attr(out, "provenance") <- list(
  package_version = as.character(utils::packageVersion("pollinatoR")),
  git_head = system("git rev-parse HEAD", intern = TRUE),
  glmmTMB = as.character(utils::packageVersion("glmmTMB")),
  TMB = as.character(utils::packageVersion("TMB")),
  Matrix = as.character(utils::packageVersion("Matrix")),
  R = R.version.string,
  created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
)

dest <- file.path("tests", "testthat", "_golden", "golden_v0.4.0_moderation.rds")
if (file.exists(dest)) {
  stop("`", dest, "` exists; the golden file is never overwritten.", call. = FALSE)
}
saveRDS(out, dest)
message("Wrote ", dest)
