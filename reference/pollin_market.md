# Declare the market inputs a value rests on

The production and price that
[`pollin_value()`](https://max578.github.io/pollinatoR/reference/pollin_value.md)
multiplies a ratio estimand by. `quantity` is current production, the
output of the crop as it is pollinated now (for example area times
yield), in `unit`; `price` is the price per `unit` in `currency`, in the
prices of `year`. The value is `ratio x price x quantity`.

## Usage

``` r
pollin_market(
  quantity,
  price,
  unit,
  currency,
  year,
  source,
  quantity_draws = NULL,
  price_draws = NULL,
  crop = NULL
)
```

## Arguments

- quantity:

  Current production, one positive number in `unit`.

- price:

  Price per `unit`, one positive number in `currency`.

- unit:

  Unit of `quantity`, one string (for example `"t"`).

- currency:

  Currency of `price`, one string (for example `"AUD"`).

- year:

  Price year, one whole number.

- source:

  Where the inputs come from, free text printed with every value.

- quantity_draws, price_draws:

  Optional numeric vectors of positive draws of `quantity` and `price`.

- crop:

  Optional name of the crop, one string. It enters the attribution of
  each value, which
  [`pollin_guard_doublecount()`](https://max578.github.io/pollinatoR/reference/pollin_guard_doublecount.md)
  reads; a value without a crop is taken to be of the same crop as any
  other.

## Value

An object of class `pollin_market`: a list holding `quantity`, `price`,
`unit`, `currency`, `year`, `source`, `quantity_draws` and `price_draws`
(`NULL` when not given), and `crop`.

## Details

Each input is a scalar, and may also be given as draws: a vector of
positive values expressing its uncertainty. When draws are given they
are used in place of the scalar, paired by index with the draws of the
estimand, so their number must equal the number of estimand draws. The
scalar is then the point value reported beside them. The package
propagates the draws given; it does not model or forecast prices.

## See also

[`pollin_value()`](https://max578.github.io/pollinatoR/reference/pollin_value.md).

## Examples

``` r
pollin_market(
  quantity = 12000,
  price = 650,
  unit = "t",
  currency = "AUD",
  year = 2024,
  source = "Illustrative inputs",
  crop = "canola"
)
#> <pollin_market>
#>   crop     : canola 
#>   quantity : 12,000 t 
#>   price    : 650 AUD per t 
#>   year     : 2024 
#>   source   : Illustrative inputs 

# Price uncertainty as draws
set.seed(1)
pollin_market(
  12000,
  650,
  "t",
  "AUD",
  2024,
  "Illustrative inputs",
  price_draws = rlnorm(4000, log(650), 0.1)
)
#> <pollin_market>
#>   quantity : 12,000 t 
#>   price    : 650 AUD per t (4000 draws) 
#>   year     : 2024 
#>   source   : Illustrative inputs 
```
