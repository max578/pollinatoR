# Guard against averaging bounded unit-level indices

The Relative Interaction Index and the dependence ratio are bounded
ratios. The mean of per-unit (per-plant, per-plot, per-site) ratios is
not the ratio of the expected values: with unequal unit means the two
differ, and the mean of ratios has no causal reading as the trial's
dependence. This guard refuses the average unless `allow = TRUE`, and
when allowed it reports the average with that statement and, from a
design, beside the ratio of means.

## Usage

``` r
pollin_guard_bounded(x, index, allow = FALSE, response = NULL, unit = NULL)
```

## Arguments

- x:

  A
  [`pollin_design()`](https://max578.github.io/pollinatoR/reference/pollin_design.md)
  with an `"excluded"` and an `"open"` arm, or a numeric vector of
  unit-level index values.

- index:

  `"rii"` or `"dependence"`. There is no default.

- allow:

  `FALSE` (the default) stops with an error; `TRUE` returns the average
  with the warning statement attached.

- response:

  With a design, the response column to use. It may be left `NULL` when
  the data declare one response.

- unit:

  With a design, the unit column the arm means are taken over, one of
  the design's unit nesting. `NULL` (the default) uses the outermost
  unit; an inner unit is taken within its outer units.

## Value

With `allow = TRUE`, an object of class `pollin_bounded`: a list holding
`index`, `n` (units averaged), `mean_of_units`, `ratio_of_means` (`NA`
for a numeric vector), `per_unit` (a data frame, or `NULL` for a numeric
vector) and the `note` printed with it.

## Details

Two inputs are accepted:

- a
  [`pollin_design()`](https://max578.github.io/pollinatoR/reference/pollin_design.md):
  unit-level arm means are computed as in
  [`pollin_rii()`](https://max578.github.io/pollinatoR/reference/pollin_rii.md)
  (excluded versus open arm, outermost unit by default), the index is
  computed for each unit, and the ratio of means is computed from the
  arm means across units;

- a numeric vector of unit-level index values already computed. The
  ratio of means cannot be recovered from these and is reported as `NA`.

The indices are `rii = (open - excluded) / (open + excluded)` and
`dependence = 1 - excluded / open`, as in
[`pollin_index()`](https://max578.github.io/pollinatoR/reference/pollin_index.md).

## See also

[`pollin_index()`](https://max578.github.io/pollinatoR/reference/pollin_index.md)
for the index of the arm means with an interval.

## Examples

``` r
df <- data.frame(
  pods = rpois(48, rep(c(20, 30, 40, 60), each = 12)),
  trt = rep(c("Cage", "Open"), 24),
  Site = rep(letters[1:4], each = 12)
)
pd <- pollin_data(df, "pods", "trt", unit = "Site")
ds <- pollin_design(pd, c(excluded = "Cage", open = "Open"), "multi_site")
# Refused by default
try(pollin_guard_bounded(ds, index = "rii"))
#> Error : Averaging bounded per-unit ratios is not the ratio of expected values; the mean of 4 unit-level RII values is not the trial's RII. The RII of the arm means is -0.0248. Set `allow = TRUE` to average them anyway.
# Allowed, with the statement and the ratio of means beside it
pollin_guard_bounded(ds, index = "rii", allow = TRUE)
#> <pollin_bounded>  average of unit-level RII values
#>   mean of 4 units : -0.0297
#>   RII of arm means : -0.0248
#>   Averaging bounded per-unit ratios is not the ratio of expected values; the mean of 4 unit-level RII values is not the trial's RII. 
```
