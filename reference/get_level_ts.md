# Extract the level component as a time series

Extract the level component as a time series

## Usage

``` r
get_level_ts(
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

A univariate `ts` object of the selected level estimate (in degrees
Celsius). If `ci = TRUE`, a multivariate `ts` object with columns
`level`, `lwr`, and `upr` is returned.

## Details

Filtered states condition on observations up to each time point, whereas
smoothed states condition on all observations in the fitted series. In
both cases, model parameters are those estimated from the complete input
series; requesting filtered states does not refit the model
sequentially.

## Examples

``` r
if (FALSE) { # \dontrun{
data(niigata_sst)
res <- tempssm(niigata_sst)
level_ts <- get_level_ts(res)
filtered_level <- get_level_ts(res, estimate = "filtered")
} # }
```
