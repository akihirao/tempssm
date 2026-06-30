# Extract the level component as a time series

Extract the level component as a time series

## Usage

``` r
get_level_ts(res, ci = FALSE, ci_level = 0.95)
```

## Arguments

- res:

  An object of class `"tempssm"` returned by
  [`tempssm()`](https://akihirao.github.io/tempssm/reference/tempssm.md).

- ci:

  Logical; if TRUE, pointwise confidence intervals are returned.

- ci_level:

  Numeric confidence level between 0 and 1 (default: 0.95).

## Value

A univariate `ts` object of the smoothed level component (in degrees
Celsius). If `ci = TRUE`, a multivariate `ts` object with columns
`level`, `lwr`, and `upr` is returned.

## Examples

``` r
if (FALSE) { # \dontrun{
data(niigata_sst)
res <- tempssm(niigata_sst)
level_ts <- get_level_ts(res)
} # }
```
