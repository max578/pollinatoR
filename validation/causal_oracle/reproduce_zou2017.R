# Independent oracle for the sham-referenced dependence contrast (criterion C7).
#
# Reproduces the treatment t values and degrees of freedom of Table 1 of
# Zou et al. (2017), BMC Ecology 17:6, doi:10.1186/s12898-017-0116-1 (CC BY
# 4.0), from the plant-level data of its Additional file 1, through
# pollin_design(reference = "sham") and pollin_effects(). See README.md.
#
# Run from the package root:
#
#   Rscript --vanilla validation/causal_oracle/reproduce_zou2017.R [xlsx]
#
# With no argument the supplement is downloaded from its published URL to a
# temporary file. Either way its SHA-256 is checked before it is read, and no
# copy of the data is kept. Exits with status 1 on any missed target.

suppressMessages(devtools::load_all(".", quiet = TRUE))

xlsx_url <- paste0(
  "https://static-content.springer.com/esm/",
  "art%3A10.1186%2Fs12898-017-0116-1/MediaObjects/",
  "12898_2017_116_MOESM1_ESM.xlsx"
)
xlsx_sha256 <- "b8b627e2f8e7136ea0277ef81b296cc0db214881919fa049ac03ea4b5824d588"
xlsx_sheet <- "Plant yield components"

# Input ------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1L) {
  xlsx <- args[1L]
  if (!file.exists(xlsx)) {
    stop("File not found: ", xlsx, call. = FALSE)
  }
} else {
  xlsx <- tempfile(fileext = ".xlsx")
  utils::download.file(xlsx_url, xlsx, mode = "wb", quiet = TRUE)
}
hash <- unname(tools::sha256sum(xlsx))
if (!identical(hash, xlsx_sha256)) {
  stop("SHA-256 mismatch: ", hash, call. = FALSE)
}

raw <- as.data.frame(readxl::read_excel(xlsx, sheet = xlsx_sheet))
expected_names <- c(
  "Site", "Treatment (O=Open cage, C=Closed Cage)", "Sample", "Pod", "Flower",
  "Seed set", "1000 Seed weight", "Straw Biomass"
)
stopifnot(identical(names(raw), expected_names))
d <- data.frame(
  field = sprintf("f%02d", as.integer(raw[[1L]])),
  cage = raw[[2L]],
  seed_set = raw[["Seed set"]],
  sqrt_pods = sqrt(raw[["Pod"]]),
  sqrt_flowers = sqrt(raw[["Flower"]]),
  log_tgw = log(raw[["1000 Seed weight"]]),
  sqrt_straw = sqrt(raw[["Straw Biomass"]]),
  stringsAsFactors = FALSE
)
stopifnot(setequal(unique(d$cage), c("C", "O")))

# Published targets (Table 1; README) -----------------------------------------

targets <- data.frame(
  response = c("seed_set", "sqrt_pods", "sqrt_flowers", "log_tgw", "sqrt_straw"),
  table_row = c(
    "Seed set", "Number of pods", "Number of flowers",
    "Thousand seed weight", "Plant straw biomass"
  ),
  t_published = c(7.86, -1.06, -2.99, -3.25, -3.72),
  df_published = c(113L, 115L, 115L, 111L, 115L),
  stringsAsFactors = FALSE
)

# The package route --------------------------------------------------------------
# The open-top cage admits pollinators under the same device as the closed cage:
# a sham arm. There is no uncaged arm, so the contrast is the dependence net of
# the cage, closed against open-top, on each transformed response.

pd <- pollin_data(d, targets$response, "cage", unit = "field")
ds <- pollin_design(pd, c(excluded = "C", sham = "O"), "rcbd",
                    reference = "sham")
est <- pollin_estimand(ds, "dependence", "difference", "loss")
fit <- pollin_fit(est, family = stats::gaussian())
eff <- pollin_effects(fit)$table
eff <- eff[match(targets$response, eff$response), ]

# The authors' model: nlme::lme, field random intercept, fitted by ML, whose
# degrees of freedom for a within-field treatment are n - fields - 1.
n_obs <- vapply(targets$response, function(r) sum(!is.na(d[[r]])), integer(1))
n_field <- vapply(targets$response, function(r) {
  length(unique(d$field[!is.na(d[[r]])]))
}, integer(1))

res <- data.frame(
  targets,
  estimate = eff$estimate,
  std.error = eff$std.error,
  t_package = eff$estimate / eff$std.error,
  df_design = n_obs - n_field - 1L,
  effect = eff$effect,
  stringsAsFactors = FALSE
)
res$t_gap <- res$t_package - res$t_published
res$t_met <- abs(res$t_gap) <= 0.005 + 1e-9
res$df_met <- res$df_design == res$df_published
res$gated <- res$response == "log_tgw"

# Cross-check, not a target: nlme's ML fit of the same model. Its coefficient
# must equal the package estimate; its standard error is conditional on the
# variance parameters, where glmmTMB's comes from the full observed
# information, so the two t values differ by a fraction of a percent.
res$coef_nlme_ml <- NA_real_
res$t_nlme_ml <- NA_real_
for (i in seq_len(nrow(res))) {
  m <- nlme::lme(stats::reformulate("cage", res$response[i]),
                 random = ~ 1 | field, data = d, na.action = stats::na.omit,
                 method = "ML")
  res$coef_nlme_ml[i] <- summary(m)$tTable[2L, "Value"]
  res$t_nlme_ml[i] <- summary(m)$tTable[2L, "t-value"]
}
res$coef_met <- abs(res$estimate - res$coef_nlme_ml) <= 1e-5 * abs(res$coef_nlme_ml)

print(res[setdiff(names(res), "table_row")], digits = 5, row.names = FALSE)
met <- all(res$t_met[res$gated]) && all(res$df_met) && all(res$coef_met) &&
  all(res$effect == "total")
cat("\nGated row (thousand seed weight) and all degrees of freedom and",
    "coefficients met:", met, "\n")
cat("Data SHA-256:", hash, "\n")
cat("pollinatoR git HEAD:", system("git rev-parse HEAD", intern = TRUE), "\n")
if (!met) {
  quit(status = 1L)
}
