# Declare the unit-of-record table of a pollination study

Pairs a data frame with a role map naming which columns are responses,
the pollination treatment, the units the treatment is applied within,
and the optional grouping, denominator, exposure and adjustment columns.
The roles are checked against the data when the object is built, and the
structural checks of
[`pollin_validate()`](https://max578.github.io/pollinatoR/reference/pollin_validate.md)
are run and stored, so every later step works from one validated
description of the table.

## Usage

``` r
pollin_data(
  data,
  response,
  treatment,
  unit,
  groups = NULL,
  denominator = NULL,
  exposures = NULL,
  adjust = NULL
)
```

## Arguments

- data:

  A data frame, one row per unit of record (for example one plant under
  one pollination treatment).

- response:

  Character vector of response column names (for example pod number,
  seed weight). Each must be numeric.

- treatment:

  Name of the pollination treatment column (for example caged versus
  open). It must have at least two distinct non-missing values.

- unit:

  Character vector of the columns identifying the units the treatment is
  applied within, outer to inner (for example `c("Site", "Plot")`).

- groups:

  Optional character vector of further grouping columns that are not
  part of the unit nesting (for example region, variety or year).

- denominator:

  Optional name of a numeric column giving what the response is counted
  out of (ovules, flowers, plants). Values must be non-negative.

- exposures:

  Optional character vector of exposure columns (for example insect
  abundances).

- adjust:

  Optional character vector of adjustment covariates (for example
  landscape or climate).

## Value

An object of class `pollin_data`: a list holding `data`, the role map
`roles`, the number of rows `n`, and the structural check table
`checks`.

## Details

`pollin_data()` records what the columns are. How the units are nested
and which treatment levels form the arms of the experiment are declared
in the design step that follows.

## See also

[`pollin_validate()`](https://max578.github.io/pollinatoR/reference/pollin_validate.md)
for the checks, with expected levels and internal-consistency
identities.

## Examples

``` r
df <- data.frame(
  pods = rpois(40, 30),
  trt = rep(c("Cage", "Open"), 20),
  Site = rep(letters[1:4], each = 10),
  bees = rpois(40, 5)
)
pd <- pollin_data(df, response = "pods", treatment = "trt", unit = "Site",
                  exposures = "bees")
pd
#> <pollin_data>
#>   rows        : 40 
#>   response    : pods 
#>   treatment   : trt ( Cage, Open )
#>   unit        : Site 
#>   exposures   : bees 
#>   checks      : 0 warnings; see `$checks`
```
