# Declare the design of a pollination experiment

Names the role of each treatment level (the arms) and the layout of the
experiment. The layout fixes the random structure used when models are
fitted: one random intercept for each level of the unit nesting declared
in
[`pollin_data()`](https://max578.github.io/pollinatoR/reference/pollin_data.md),
outer to inner, with an optional random treatment slope on one unit
level, giving a pollination effect that varies between those units.

## Usage

``` r
pollin_design(x, arms, layout, reference = "open", slope = NULL)
```

## Arguments

- x:

  A
  [`pollin_data()`](https://max578.github.io/pollinatoR/reference/pollin_data.md)
  object.

- arms:

  Named character vector mapping arm roles to treatment levels, for
  example `c(excluded = "Cage", open = "Open")`. Roles are `"excluded"`,
  `"open"` and `"supplemented"`; every treatment level must be mapped.

- layout:

  One of `"crd"`, `"rcbd"`, `"split_plot"`, `"matched_pairs"` or
  `"multi_site"`. There is no default.

- reference:

  The arm that dependence is measured against: `"open"` (the default) or
  `"supplemented"`. It must be one of the declared arms.

- slope:

  The unit level carrying a random treatment slope. `NULL` (the default)
  uses the layout's default from the table above; `FALSE` fits no slope;
  a column name from the unit nesting puts the slope there.

## Value

An object of class `pollin_design`: a list holding the
[`pollin_data()`](https://max578.github.io/pollinatoR/reference/pollin_data.md)
object `data`, `arms`, `layout`, `reference`, `slope` and the
random-effect terms `random`.

## Details

|  |  |  |  |
|----|----|----|----|
| Layout | Units, outer to inner | Treatment varies | Default slope |
| `"crd"` | experimental unit (and subsamples within it) | between units | none |
| `"rcbd"` | block, then any inner levels | within blocks | none |
| `"split_plot"` | block, whole plot, then any inner levels | within whole plots | none |
| `"matched_pairs"` | any outer levels, then the pair | within pairs | none |
| `"multi_site"` | site, then any inner levels | within sites | site |

In a split-plot layout the pollination treatment is taken to be the
sub-plot factor (for example cages within variety plots). The design
checks that every arm appears in every block, whole plot, pair or site,
and warns with a count when it does not.

## Examples

``` r
df <- data.frame(
  pods = rpois(48, 30),
  trt = rep(c("Cage", "Open"), 24),
  Site = rep(letters[1:4], each = 12),
  Plot = rep(rep(c("p1", "p2"), each = 6), 4)
)
pd <- pollin_data(df, "pods", "trt", unit = c("Site", "Plot"))
pollin_design(pd, arms = c(excluded = "Cage", open = "Open"),
              layout = "multi_site")
#> <pollin_design>
#>   layout    : multi_site 
#>   arms      : excluded = Cage, open = Open 
#>   reference : open 
#>   units     : Site / Plot 
#>   random    : (trt | Site) + (1 | Site:Plot) 
```
