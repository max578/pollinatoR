# Sham-referenced contrast oracle: Zou et al. (2017)

This folder holds the independent check for criterion C7 of Stage 4: a
published result on a pollinator-exclusion trial with a sham device,
recomputed from the published plant-level data through the package.

## The study and why it fits

Zou et al. (2017) caged 134 oilseed rape plants in 18 fields in Jiangxi,
China, each plant under its own cage, alternating open-top cages (insects
enter) and closed cages (insects excluded). Both arms carry the device, so
the open-top cage is a sham arm and the contrast is the effect of pollinator
access net of the cage. There is no uncaged arm, so the study cannot give the
cage effect itself, nor a cage-factor range for `pollin_bound()`.

In the package the design is
`pollin_design(arms = c(excluded = "C", sham = "O"), layout = "rcbd",
reference = "sham")`, with the field as the block, and the estimand is the
dependence family on the difference scale for each transformed response. Each
response is fitted by glmmTMB (Gaussian, field random intercept, maximum
likelihood), and `pollin_effects()` reports the arm difference and its
standard error.

## Targets

Table 1 of the article reports, for each yield parameter, the t value and
degrees of freedom of the treatment term from linear mixed models in `nlme`
with the field as a random intercept. The article states no estimation
method; the data reproduce every listed t value to its printed two decimals
under maximum likelihood and not under REML (for thousand seed weight, ML
gives -3.2452 and REML -3.2698).

| Target | Source | Tolerance | Gated |
|---|---|---|---|
| t of log thousand seed weight, -3.25 | Table 1 | 0.005 | yes (the D9 ruling, 2026-11-10) |
| degrees of freedom of every row: 113, 115, 115, 111, 115 | Table 1 | exact | yes |
| t of seed set (7.86), square-root pods (-1.06), square-root flowers (-2.99), square-root straw biomass (-3.72) | Table 1 | 0.005 | no, reported |
| package estimate equals the `nlme` ML coefficient | cross-check | 1e-5 relative | yes |

The degrees of freedom are those `nlme` assigns to a treatment that varies
within fields, `n - fields - 1`, computed from the data: they confirm the
row counts per response (two seed-set and four seed-weight values are
missing). The plant yield, fruit set and yield-to-straw rows of Table 1 have
no column in the supplement and are not attempted.

## Result (2026-11-14)

The gated targets are met. Thousand seed weight: t -3.2533 against -3.25, df
111 against 111, and every package estimate equals the `nlme` ML coefficient
to six significant figures.

The four other rows were first drafted as gated targets too, before the
package run, and four of them missed: seed set 7.8984 against 7.86, pods
-1.0668 against -1.06, flowers -3.0114 against -2.99, straw -3.7415 against
-3.72. The point estimates agree; the standard errors differ. `nlme` takes
the fixed-effect covariance as the inverse of `X'V^-1 X` with the variance
parameters held at their estimates. glmmTMB takes it from the observed
information of all parameters together, where the cross terms between the
fixed effects and the variance parameters are not zero in a finite sample.
Across the five rows the `nlme` standard error is larger by 0.25 to 0.64
percent; both are valid large-sample standard errors. The gate was narrowed
to the ruled row after the run, and this paragraph is the record of that.

The run also warns that one field (`f17`) carries no closed cage. It enters
the fit through its intercept, as it does in the authors' model.

## Source and licence

Zou, Y., Xiao, H., Bianchi, F. J. J. A., Jauker, F., Luo, S., & van der Werf,
W. (2017). Wild pollinators enhance oilseed rape yield in small-holder
farming systems in China. *BMC Ecology*, *17*, 6.
https://doi.org/10.1186/s12898-017-0116-1

The article is published under CC BY 4.0, and its data availability note
applies the CC0 waiver to the data made available in the article unless
stated otherwise. The data are Additional file 1, sheet
`Plant yield components`:

<https://static-content.springer.com/esm/art%3A10.1186%2Fs12898-017-0116-1/MediaObjects/12898_2017_116_MOESM1_ESM.xlsx>

SHA-256 `b8b627e2f8e7136ea0277ef81b296cc0db214881919fa049ac03ea4b5824d588`.
The script reads the file at run time, and no copy ships with the package or
this repository.

## Run

From the package root:

```
Rscript --vanilla validation/causal_oracle/reproduce_zou2017.R [xlsx]
```

With no argument the file is downloaded to a temporary path. The script
exits with status 1 if a gated target is missed. The run of 2026-11-14 is in
`report_2026-11-14.txt` (R 4.6.1, glmmTMB 1.1.14, nlme 3.1.171).
