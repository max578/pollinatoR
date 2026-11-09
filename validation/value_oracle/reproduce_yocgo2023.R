# Independent oracle for the dependence channel (criterion C4).
#
# Reproduces the insect-pollination economic value (IPEV) of Rwandan crop
# production in Table 1 of Yocgo et al. (2023), Scientific Reports 13:20108,
# doi:10.1038/s41598-023-46936-w (CC BY 4.0), by valuing every crop-year of
# the published supplement through pollin_value(channel = "dependence") and
# adding the crops of a year with pollin_guard_doublecount(). See README.md.
#
# Run from the package root:
#
#   Rscript --vanilla validation/value_oracle/reproduce_yocgo2023.R [xlsx]
#
# With no argument the supplement is downloaded from its published URL to a
# temporary file. Either way its SHA-256 is checked before it is read, and no
# copy of the data is kept. Exits with status 1 on any missed target.

suppressMessages(devtools::load_all(".", quiet = TRUE))

xlsx_url <- paste0(
  "https://static-content.springer.com/esm/",
  "art%3A10.1038%2Fs41598-023-46936-w/MediaObjects/",
  "41598_2023_46936_MOESM1_ESM.xlsx"
)
xlsx_sha256 <- "3ee00fa06239b498d477eb2c4079dab8d90c70d21c07bc2ae7b7769cfebbf9f0"
xlsx_sheet <- "Rwanda_cleaned crop data and ca"
source_label <- "Yocgo et al. 2023 supplement"

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
  stop(
    "SHA-256 mismatch: expected ",
    xlsx_sha256,
    ", got ",
    hash,
    ". The supplement is not the file the targets were set against.",
    call. = FALSE
  )
}

tab <- as.data.frame(readxl::read_excel(
  xlsx,
  sheet = xlsx_sheet,
  skip = 1L,
  na = c("", "NA"),
  .name_repair = "unique_quiet"
))
if (length(args) == 0L) {
  unlink(xlsx) # the download is not kept
}
need <- c(
  "FAO_Crop",
  "Year",
  "Crop_type",
  "GroupDi",
  "Mean_Di",
  "Production.tons",
  "Value.USD/ton"
)
miss <- setdiff(need, names(tab))
if (length(miss)) {
  stop("Columns not found: ", paste(miss, collapse = ", "), call. = FALSE)
}
tab$D <- as.numeric(tab$Mean_Di)
tab$D[is.na(tab$D)] <- 0 # indirect and non-dependent crops, as in the paper
if (any(tab$D < 0 | tab$D >= 1)) {
  stop("`Mean_Di` outside [0, 1).", call. = FALSE)
}
tab$Q <- tab[["Production.tons"]]
tab$P <- tab[["Value.USD/ton"]]

# Valuation through the package ----------------------------------------------

# A minimal dependence / ratio / loss estimand from the public API. Its data
# only declare the contrast; the arm means come from the draws below.
toy <- data.frame(
  y = c(1, 1, 1, 1),
  trt = c("Cage", "Open", "Cage", "Open"),
  Site = c("a", "a", "b", "b")
)
estimand <- pollin_estimand(
  pollin_design(
    pollin_data(toy, "y", "trt", unit = "Site"),
    c(excluded = "Cage", open = "Open"),
    "multi_site"
  ),
  family = "dependence",
  scale = "ratio",
  basis = "loss"
)

# Degenerate draws whose arm means give D = 1 - low/high: high = 1,
# low = 1 - D on n_draws identical draws.
n_draws <- 4L
draws_for <- function(D) {
  high <- matrix(1, nrow = n_draws, ncol = 1L, dimnames = list(NULL, "overall"))
  low <- matrix(1 - D, nrow = n_draws, ncol = 1L, dimnames = list(NULL, "overall"))
  pollinatoR:::.new_pollin_draws(
    high = high,
    low = low,
    estimand = estimand,
    response = "y",
    by = NULL,
    method = "simulation",
    engine = "oracle",
    n = n_draws,
    seed = NULL,
    population = NULL
  )
}

value_row <- function(i) {
  r <- tab[i, ]
  pollin_value(
    draws_for(r$D),
    pollin_market(
      quantity = r$Q,
      price = r$P,
      unit = "t",
      currency = "USD",
      year = r$Year,
      source = source_label,
      crop = r$FAO_Crop
    ),
    channel = "dependence",
    product = "crop"
  )
}
values <- lapply(seq_len(nrow(tab)), value_row)
names(values) <- tab$FAO_Crop

# Per-crop value: every draw must be the same number.
one_value <- function(d) {
  if (any(d != d[1L])) {
    stop("Degenerate draws gave unequal value draws.", call. = FALSE)
  }
  d[1L]
}
tab$ipev <- vapply(values, function(v) one_value(v$draws), numeric(1))

# Sum of the crops in a subset through the double-count guard. A single
# crop is its own value; the guard needs two or more.
guard_total <- function(idx) {
  if (length(idx) == 1L) {
    return(tab$ipev[idx])
  }
  total <- do.call(pollin_guard_doublecount, values[idx])
  one_value(total$draws)
}

# Targets (Table 1 and text of the article) -----------------------------------

years <- c(2014L, 2020L)
published <- list(
  ipev = c("2014" = 186242200, "2020" = 119515500),
  ev = c("2014" = 2550993000, "2020" = 2787960000),
  vr = c("2014" = 7.3, "2020" = 4.3),
  groups = data.frame(
    crop_type = c(
      "Pulses",
      "Vegetables",
      "Nuts/oilseeds",
      "Fruits",
      "Stimulants",
      "Spices/herbs"
    ),
    y2014 = c(12.2e6, 153.0e6, 4.1e6, 10.8e6, 5.4e6, 767.5e3),
    y2020 = c(15.2e6, 94.5e6, 4.9e6, 3.0e6, 1.3e6, 655.2e3),
    unit = c(1e6, 1e6, 1e6, 1e6, 1e6, 1e3),
    stringsAsFactors = FALSE
  ),
  change = -35.8
)

rows <- list()
# A target passes when |reproduced - published| <= tolerance. The slack of
# 1e-9 relative to the tolerance absorbs binary floating-point error only: the
# Spices/herbs 2020 row is an exact tie (655.25 K printed as 655.2 K), whose
# difference is half a unit in the last printed place up to that error.
add <- function(target, published, reproduced, tolerance, scale = 1) {
  diff <- reproduced - published
  slack <- 1e-9 * tolerance
  rows[[length(rows) + 1L]] <<- data.frame(
    target = target,
    published = published / scale,
    reproduced = reproduced / scale,
    abs_diff = abs(diff) / scale,
    tolerance = tolerance / scale,
    pass = abs(diff) <= tolerance + slack,
    stringsAsFactors = FALSE
  )
}

ipev_year <- numeric(0)
for (yr in years) {
  y <- as.character(yr)
  in_year <- which(tab$Year == yr)
  ipev <- guard_total(in_year)
  ev <- sum(tab$P[in_year] * tab$Q[in_year])
  ipev_year[y] <- ipev
  add(paste("IPEV", y, "total (USD)"), published$ipev[[y]], ipev, 1)
  add(paste("EV", y, "total (USD)"), published$ev[[y]], ev, 1)
  add(paste("VR", y, "(%)"), published$vr[[y]], 100 * ipev / ev, 0.05)
  g <- published$groups
  for (k in seq_len(nrow(g))) {
    idx <- in_year[
      tab$GroupDi[in_year] == "Direct" & tab$Crop_type[in_year] == g$crop_type[k]
    ]
    add(
      paste0("IPEV ", y, " Direct ", g$crop_type[k], if (g$unit[k] == 1e6) " (M USD)" else " (K USD)"),
      g[[paste0("y", y)]][k],
      guard_total(idx),
      0.05 * g$unit[k],
      scale = g$unit[k]
    )
  }
}
add(
  "IPEV change 2014 to 2020 (%)",
  published$change,
  100 * (ipev_year[["2020"]] / ipev_year[["2014"]] - 1),
  0.05
)

# Per-row agreement with the authors' own IPEV column, a diagnostic only:
# that column is an output of the article, not an input of this script.
authors <- suppressWarnings(as.numeric(tab$IPEV.USD))
authors[is.na(authors)] <- 0
add(
  "Largest per-row |IPEV - authors' IPEV.USD| (USD)",
  0,
  max(abs(tab$ipev - authors)),
  0.01
)

# Report -----------------------------------------------------------------------

out <- do.call(rbind, rows)
fmt <- out
for (col in c("published", "reproduced", "abs_diff", "tolerance")) {
  fmt[[col]] <- formatC(out[[col]], format = "f", digits = 4, big.mark = ",")
}
fmt$pass <- ifelse(out$pass, "yes", "NO")
options(width = 200L)
cat("Supplement SHA-256:", hash, "\n")
cat("Crop-years valued:", nrow(tab), "\n\n")
print(fmt, row.names = FALSE, right = FALSE)
cat("\n")
if (all(out$pass)) {
  cat("C4 verdict: PASS (", nrow(out), " of ", nrow(out), " targets met)\n", sep = "")
} else {
  cat(
    "C4 verdict: FAIL (",
    sum(!out$pass),
    " of ",
    nrow(out),
    " targets missed)\n",
    sep = ""
  )
  quit(save = "no", status = 1L)
}
