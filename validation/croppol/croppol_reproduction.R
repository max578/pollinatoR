# CropPol reproduction, Stage 1 (criterion C4).
#
# LICENCE: this file and everything in validation/croppol/ are products of the
# CropPol database, released under CC BY-NC-SA 4.0 for non-commercial research
# use (https://creativecommons.org/licenses/by-nc-sa/4.0/). The MIT licence of
# the pollinatoR package does not cover this directory, and `^validation$` in
# .Rbuildignore keeps it out of the built package. See LICENCE.md here.
# Data source: Allen-Perkins et al. (2022). CropPol: A dynamic, open and global
# database on crop pollination. Ecology 103(3), e3614.
# https://doi.org/10.1002/ecy.3614
#
# No CropPol data are stored in this repository. The field-level table is read
# at run time from GitHub, pinned to one commit, checked against its MD5, and
# cached only under tempdir().
#
# Run from the package root:
#
#   Rscript --vanilla validation/croppol/croppol_reproduction.R
#
# Writes validation/croppol/croppol_reproduction_result.csv.
#
# Target (a number we did not produce) -----------------------------------------
#
# Zou, Y., Xiao, H., Bianchi, F. J. J. A., Jauker, F., Luo, S., & van der Werf,
# W. (2017). Wild pollinators enhance oilseed rape yield in small-holder
# farming systems in China. BMC Ecology 17, 6.
# https://doi.org/10.1186/s12898-017-0116-1
# Results: "Oilseed rape plants in closed cages had 38% (±4%) lower seed set"
# than plants in open cages; seed set is the number of seeds per pod; "Means and
# standard errors are reported throughout the paper", so 4% is a standard error.
#
# Pass rule (fixed before the first run) ----------------------------------------
#
# Estimand: dependence (E2), excluded versus open, ratio scale, loss basis,
# D = 1 - mean(excluded) / mean(open), on seeds per pod.
# Published interval: the 95% normal interval implied by the published value
# and standard error, 0.38 +/- 1.96 * 0.04 = [0.3016, 0.4584].
# PASS when our point estimate lies inside that interval; FAIL otherwise.
# The narrower +/- 1 SE band [0.34, 0.42] is reported alongside, not used for
# the verdict.

suppressMessages(devtools::load_all(".", quiet = TRUE))

# Read CropPol at run time -----------------------------------------------------

commit <- "febec7dcc088ccf1bd721c41ce8dda87a2dcba1a"
url <- paste0(
  "https://raw.githubusercontent.com/ibartomeus/CropPol/",
  commit,
  "/Final_Data/CropPol_field_level_data.csv"
)
md5_expected <- "63521cdbf6080c4e79da97ae2c907ff1"
cache <- file.path(tempdir(), "CropPol_field_level_data.csv")
if (!file.exists(cache)) {
  utils::download.file(url, cache, mode = "wb", quiet = TRUE)
}
md5 <- unname(tools::md5sum(cache))
if (!identical(md5, md5_expected)) {
  stop("CropPol file MD5 ", md5, " differs from ", md5_expected, ".")
}
cp <- utils::read.csv(cache, stringsAsFactors = FALSE, encoding = "UTF-8")

# The study: one row per field, open-cage and closed-cage seeds per pod ---------

study <- "Yi_Zou_Brassica_napus_China_2015"
z <- cp[cp$study_id == study, ]
stopifnot(
  nrow(z) > 0L,
  all(z$crop == "Brassica napus"),
  all(z$yield2_units == "seed set")
)
z <- z[!is.na(z$yield2) & !is.na(z$yield_treatments_no_pollinators2), ]
n_fields <- nrow(z)

# Long table: two rows per field (open, excluded). Only the study, field and
# seed-set columns are kept; contact and credit columns are not carried.
long <- data.frame(
  field = rep(z$site_id, 2L),
  cage = rep(c("open", "excluded"), each = n_fields),
  seeds_per_pod = c(z$yield2, z$yield_treatments_no_pollinators2),
  stringsAsFactors = FALSE
)

# Layered pipeline ---------------------------------------------------------------

pd <- pollin_data(
  long,
  response = "seeds_per_pod",
  treatment = "cage",
  unit = "field"
)
ds <- pollin_design(
  pd,
  arms = c(excluded = "excluded", open = "open"),
  layout = "matched_pairs"
)
est <- pollin_estimand(ds, family = "dependence", scale = "ratio", basis = "loss")
fit <- pollin_fit(est)
eff <- pollin_effects(fit)
print(pd)
print(ds)
print(est)
print(fit)
print(eff)

# Compare with the published interval ---------------------------------------------

published <- 0.38
published_se <- 0.04
pub_low <- published - 1.96 * published_se
pub_high <- published + 1.96 * published_se
row <- eff$table[1L, ]
verdict <- if (row$estimate >= pub_low && row$estimate <= pub_high) {
  "PASS"
} else {
  "FAIL"
}
inside_1se <- row$estimate >= published - published_se &&
  row$estimate <= published + published_se

# Descriptive ratio of raw means, for the report (not the verdict)
raw_ratio_loss <- 1 - mean(z$yield_treatments_no_pollinators2) / mean(z$yield2)

result <- data.frame(
  source = "CropPol field-level data",
  commit = commit,
  md5 = md5,
  study_id = study,
  n_fields = n_fields,
  response = "seeds per pod (CropPol yield2, 'seed set')",
  estimand = "dependence, ratio scale, loss basis: 1 - excluded/open",
  model = paste(
    "seeds_per_pod ~",
    paste(fit$rhs, collapse = " + "),
    "|",
    fit$families[["seeds_per_pod"]]$family,
    "(",
    fit$families[["seeds_per_pod"]]$link,
    ")"
  ),
  estimate = row$estimate,
  conf.low = row$conf.low,
  conf.high = row$conf.high,
  interval = paste0(100 * eff$conf_level, "% ", eff$interval),
  raw_mean_ratio_loss = raw_ratio_loss,
  published = published,
  published_se = published_se,
  published_low_95 = pub_low,
  published_high_95 = pub_high,
  inside_published_1se = inside_1se,
  verdict = verdict,
  r_version = R.version.string,
  pollinatoR = as.character(utils::packageVersion("pollinatoR")),
  glmmTMB = as.character(utils::packageVersion("glmmTMB")),
  marginaleffects = as.character(utils::packageVersion("marginaleffects")),
  run_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  stringsAsFactors = FALSE
)
utils::write.csv(
  result,
  file.path("validation", "croppol", "croppol_reproduction_result.csv"),
  row.names = FALSE
)
cat("\nVerdict:", verdict, "\n")
cat(sprintf(
  "Ours %.4f [%.4f, %.4f]; published %.2f, 95%% [%.4f, %.4f]\n",
  row$estimate,
  row$conf.low,
  row$conf.high,
  published,
  pub_low,
  pub_high
))
