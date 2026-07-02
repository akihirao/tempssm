# Extract the seasonal component as a time series

Extract the seasonal component as a time series

## Usage

``` r
get_season_ts(
  res,
  ci = FALSE,
  ci_level = 0.95,
  estimate = c("smoothed", "filtered")
)
```

## Arguments

- res:

  An object of class `"tempssm"` returned by
  [`tempssm()`](https://akihirao.github.io/tempssm/reference/tempssm.md).

- ci:

  Logical; if TRUE, pointwise confidence intervals are returned.

- ci_level:

  Numeric confidence level between 0 and 1 (default: 0.95).

- estimate:

  Character scalar specifying the state estimate to return. Use
  `"smoothed"` (the default) for estimates conditional on all
  observations, or `"filtered"` for estimates conditional on
  observations up to each time point. Filtered confidence intervals are
  not currently supported, so `estimate = "filtered"` requires
  `ci = FALSE`.

## Value

A univariate `ts` object of the selected seasonal estimate (in degrees
Celsius). If `ci = TRUE`, a multivariate `ts` object with columns
`season`, `lwr`, and `upr` is returned.

## Details

The seasonal component represents recurrent intra-year variability
captured by seasonal dummy state components in the state space model.

## Examples

``` r
if (FALSE) { # \dontrun{
data(niigata_sst)
res <- tempssm(niigata_sst)
season_ts <- get_season_ts(res)
} # }
```
