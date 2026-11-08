# Benjamini-Hochberg false-discovery-rate adjustment

Thin wrapper over
[`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) with the BH
method, kept as a named verb so the multiplicity step is explicit in
analysis code.

## Usage

``` r
pollin_fdr(p)
```

## Arguments

- p:

  Numeric vector of p-values.

## Value

Numeric vector of FDR-adjusted p-values.

## Examples

``` r
pollin_fdr(c(0.001, 0.02, 0.30, 0.80))
#> [1] 0.004 0.040 0.400 0.800
```
