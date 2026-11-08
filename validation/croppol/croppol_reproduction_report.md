# CropPol reproduction (criterion C4)

**Verdict: PASS.** The dependence of oilseed rape seed set on insect pollination,
estimated by the pollinatoR pipeline from CropPol, is 0.395 (95% delta interval
0.317 to 0.465). The published value is 0.38 with standard error 0.04, whose 95%
interval is 0.302 to 0.458. Our point estimate lies inside that interval, and
also inside the narrower one-standard-error band (0.34 to 0.42).

This directory is CC BY-NC-SA 4.0; see `LICENCE.md`.

## Source

| Item | Value |
|---|---|
| Database | CropPol field-level table, `Final_Data/CropPol_field_level_data.csv` |
| URL (pinned) | `https://raw.githubusercontent.com/ibartomeus/CropPol/febec7dcc088ccf1bd721c41ce8dda87a2dcba1a/Final_Data/CropPol_field_level_data.csv` |
| Commit | `febec7dcc088ccf1bd721c41ce8dda87a2dcba1a` (2023-05-07, the latest commit touching this file when the script was written) |
| MD5 | `63521cdbf6080c4e79da97ae2c907ff1`, checked by the script on every run |
| Study | `Yi_Zou_Brassica_napus_China_2015`, 17 fields with both arms |
| Data citation | Allen-Perkins, A., et al. (2022). CropPol: A dynamic, open and global database on crop pollination. *Ecology*, *103*(3), e3614. <https://doi.org/10.1002/ecy.3614> |

The script reads the table at run time and caches it under `tempdir()` only. No
copy of the data is kept in the repository.

## Licence evidence

| Source | Statement |
|---|---|
| GitHub `LICENSE`, <https://github.com/ibartomeus/CropPol/blob/master/LICENSE> | MIT License, "Copyright (c) 2020 BartomeusLab"; worded for "this software and associated documentation files" |
| Repository metadata `Metadata/data/metadata/biblio.csv`, `dataspice.json` | `license: "CC By"` |
| Repository `Metadata/docs/EML_CP.xml` | `intellectualRights`: "CC-by"; the abstract in the same file states the CC BY-NC-SA terms below |
| Zenodo concept record <https://doi.org/10.5281/zenodo.4311291>, all five versions (v1.0.0 to v2.1.0) | licence id `other-open` |
| Data paper abstract, <https://doi.org/10.1002/ecy.3614> | "This data set is released for non-commercial use only. Credits should be given to this paper (i.e., proper citation), and the products generated with this database should be shared under the same license terms (CC BY-NC-SA)." |

The sources conflict. This work follows the strictest reading, CC BY-NC-SA, for
non-commercial research use, with this directory licensed to match and kept out
of the built package by `^validation$` in `.Rbuildignore`. A query to the
maintainers on which statement governs the data has been drafted.

## Published number

Zou, Y., Xiao, H., Bianchi, F. J. J. A., Jauker, F., Luo, S., & van der Werf, W.
(2017). Wild pollinators enhance oilseed rape yield in small-holder farming
systems in China. *BMC Ecology*, *17*, 6.
<https://doi.org/10.1186/s12898-017-0116-1> (citation checked against Crossref).

Results section: "Oilseed rape plants in closed cages had 38% (±4%) lower seed
set" than plants in open cages. Seed set is defined in the Methods as the number
of seeds per pod. The paper states "Means and standard errors are reported
throughout the paper", so 4% is read as a standard error.

| | Value |
|---|---|
| Published estimate | 0.38 (loss share, closed versus open) |
| Published standard error | 0.04 |
| Published 95% interval (0.38 +/- 1.96 x 0.04) | 0.3016 to 0.4584 |
| Published +/- 1 SE band | 0.34 to 0.42 |

The paper does not print a confidence interval; the 95% interval is the normal
interval implied by its estimate and standard error. That rule, and the PASS
condition (our point estimate inside the 95% interval), were written into the
script before the first run.

## Our estimate

| | Value |
|---|---|
| Pipeline | `pollin_data()` -> `pollin_design(layout = "matched_pairs")` -> `pollin_estimand("dependence", "ratio", "loss")` -> `pollin_fit()` -> `pollin_effects()` |
| Model | `seeds_per_pod ~ cage + (1 \| field)`, Gamma (log link), `glmmTMB` 1.1.14, family chosen by the package |
| Estimate | 0.3954 |
| 95% delta interval | 0.3166 to 0.4650 |
| Raw ratio of field means, 1 - mean(excluded)/mean(open) | 0.386 (descriptive, not used for the verdict) |
| Verdict | **PASS** (inside the published 95% interval and the +/- 1 SE band) |

Run: R 4.6.1, pollinatoR 0.0.0.9000 (branch `stage1-crop-core`),
marginaleffects 0.32.0. Machine-readable result:
`croppol_reproduction_result.csv`.

## Differences in method

| Aspect | Zou et al. (2017) | This reproduction |
|---|---|---|
| Fields | 18 | 17 (the fields in CropPol with both arms) |
| Unit of analysis | Plants (eight caged plants per field, open and closed alternating) | Field means of seeds per pod as stored in CropPol, one open and one excluded value per field |
| Model | Linear mixed model (`nlme`), treatment fixed, field random; Gaussian for seed set | Gamma GLMM with log link (`glmmTMB`), treatment fixed, field random intercept |
| Effect scale | Percentage difference; the exact formula is not stated in the text read | Ratio of G-computed arm means, 1 - excluded/open |
| Uncertainty | Standard error | Delta interval on log(open/excluded), fixed effects only; `marginaleffects` warns that random-effect uncertainty is not included |

The two analyses share data but not data level, model or effect formula. The
agreement shows the pipeline recovers a published dependence value from the
same study's data; it does not test the pipeline against plant-level data.

## Notes from the run

- `marginaleffects` warned that for this model type it takes into account only
  the uncertainty in fixed-effect parameters. The package documents this for
  the delta interval.
- The run timestamp in the CSV carries the machine's time zone label.
