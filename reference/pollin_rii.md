# Relative Interaction Index with a bootstrap interval

The Relative Interaction Index (Armas et al. 2004) is widely reported
for pollinator dependence as a bare point estimate, despite being a
ratio of two correlated sums whose sampling distribution is awkward.
This function returns the index with a nonparametric bootstrap interval,
the minimum upgrade over a point value. Prefer the model-based
contribution
([`pollin_contribution()`](https://max578.github.io/pollinatoR/reference/pollin_contribution.md))
as the primary estimate and keep this for comparability with the
existing literature.

## Usage

``` r
pollin_rii(
  open,
  cage,
  n_boot = 2000L,
  conf = 0.95,
  seed = NULL,
  response = NULL,
  unit = NULL
)
```

## Arguments

- open:

  Numeric vector of the open-arm (pollinator-accessible) values, or a
  [`pollin_design()`](https://max578.github.io/pollinatoR/reference/pollin_design.md)
  with an `"excluded"` and an `"open"` arm.

- cage:

  Numeric vector of the exclusion-arm values. Leave it empty when `open`
  is a design.

- n_boot:

  Number of bootstrap resamples.

- conf:

  Confidence level.

- seed:

  Optional integer seed for reproducibility.

- response:

  With a design, the response column to use. It may be left `NULL` when
  the data declare one response.

- unit:

  With a design, the unit column the arm means are taken over, one of
  the design's unit nesting. `NULL` (the default) uses the outermost
  unit; an inner unit is taken within its outer units.

## Value

An object of class `pollin_rii`.

## Details

The arm values are given either as two vectors or as a
[`pollin_design()`](https://max578.github.io/pollinatoR/reference/pollin_design.md).
From a design, the values are unit-level arm means: the mean response of
each unit in the open arm and in the excluded arm, for the units
carrying both arms, in sorted unit order. The bootstrap then runs on
those means as it does on two vectors.

## Examples

``` r
set.seed(1)
pollin_rii(rpois(50, 8), rpois(50, 6), n_boot = 500, seed = 1)
#> RII = 0.160  [95% boot CI: 0.105, 0.219]  (n_open=50, n_cage=50)

# From a design: one mean per site and arm
df <- data.frame(
  pods = rpois(48, 30),
  trt = rep(c("Cage", "Open"), 24),
  Site = rep(letters[1:4], each = 12)
)
pd <- pollin_data(df, "pods", "trt", unit = "Site")
ds <- pollin_design(pd, c(excluded = "Cage", open = "Open"), "multi_site")
pollin_rii(ds, n_boot = 200, seed = 1)
#> RII = -0.023  [95% boot CI: -0.081, 0.021]  (n_open=4, n_cage=4)
```
