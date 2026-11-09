# Dependence-channel oracle: Yocgo et al. (2023)

This folder holds the independent check for criterion C4: a published
dependence-ratio valuation, recomputed from its published inputs through the
package and matched within its stated rounding.

## What is reproduced

Yocgo et al. (2023) value the contribution of insect pollinators to Rwandan
crop production, 2014 to 2020, as

- `IPEV = sum(D_i x P_i x Q_i)`, the insect-pollination economic value;
- `EV = sum(P_i x Q_i)`, the economic value of production;
- `VR = 100 x IPEV / EV`, the vulnerability ratio (%),

with `D_i` the midpoint of each crop's dependence class (Klein et al., 2007:
0.95, 0.65, 0.25, 0.05, with indirect and non-dependent crops at 0), `Q_i` the
FAOSTAT production in tonnes and `P_i` the price in USD per tonne.

`reproduce_yocgo2023.R` values each of the 224 crop-years (32 crops by 7
years) through the package. For each crop it builds draws of a dependence
estimand on the ratio scale with the loss basis whose arm means give the
published `D` (`high = 1`, `low = 1 - D`, on identical draws), declares the
market with `pollin_market()` and values it with
`pollin_value(channel = "dependence")`. The crops of a year, and the crops of
each crop-type group, are added with `pollin_guard_doublecount()`. Distinct
crops do not overlap, so this exercises the guard's summing path. The
results are compared with the article's Table 1 for 2014 and 2020:

| Target | Tolerance |
|---|---|
| IPEV and EV grand totals (printed to the dollar) | less than 1 USD |
| VR (printed to 0.1 %) | 0.05 percentage points |
| IPEV of the directly dependent crops by crop type (printed to 0.1 M or 0.1 K USD) | half a unit in the last printed place, inclusive |
| Change in IPEV, 2014 to 2020 (text, -35.8 %) | 0.05 percentage points |

The Spices/herbs row for 2020 is an exact rounding tie (655.25 K printed as
655.2 K), so the group comparison is inclusive (`<=`), with a slack of
`1e-9` of the tolerance for binary floating-point error. As a diagnostic the
script also reports the largest per-row difference from the authors' own
`IPEV.USD` column. That column is an output of the article and is not used
as an input.

The script exits with status 1 if any target is missed.

## Source and licence

Yocgo, R. E. E., Hitimana, I., Hakizimana, M., & Birachi, E. A. (2023).
Insect pollinators can unlock an annual monetary value of more than US $100
million from crop production in Rwanda. *Scientific Reports*, *13*, 20108.
https://doi.org/10.1038/s41598-023-46936-w

The article is published under CC BY 4.0. The inputs are read from its
Supplementary Information 1, sheet `Rwanda_cleaned crop data and ca`:

<https://static-content.springer.com/esm/art%3A10.1038%2Fs41598-023-46936-w/MediaObjects/41598_2023_46936_MOESM1_ESM.xlsx>

SHA-256 `3ee00fa06239b498d477eb2c4079dab8d90c70d21c07bc2ae7b7769cfebbf9f0`.
The script reads the file at run time, and no copy ships with the package or
this repository.

## How to run

From the package root, with `devtools` and `readxl` installed:

```sh
# Download the supplement from its published URL to a temporary file
Rscript --vanilla validation/value_oracle/reproduce_yocgo2023.R

# Or use a copy already on disk
Rscript --vanilla validation/value_oracle/reproduce_yocgo2023.R path/to/41598_2023_46936_MOESM1_ESM.xlsx
```

Either way the file's SHA-256 is checked before it is read, and the script
stops on a mismatch. A downloaded file is deleted once read.

## What the oracle tests, and what it does not

It tests the dependence channel's arithmetic against an authority the
package did not produce: that `pollin_value()` forms `D = 1 - low/high` from
the arm means and multiplies it by price and quantity draw by draw, and that
`pollin_guard_doublecount()` adds values of distinct crops correctly.

It does not test the estimation of `D` from trial data. The paper's `D` are
class midpoints, entered here as degenerate draws, so the fitting, the
draws of the arm means and the propagation of their uncertainty are outside
its reach. Nor does it test prices: the supplement's price is the authors'
gross production value divided by production, so `P x Q` equals the FAOSTAT
gross production value by construction, and the EV rows confirm that
identity rather than an independent price.

The estimand used is a minimal synthetic declaration (four rows, two
treatments) that carries the family, scale and basis only. The draws are
built with the internal constructor `.new_pollin_draws()`, since the public
`pollin_draws()` draws from a fitted model.
