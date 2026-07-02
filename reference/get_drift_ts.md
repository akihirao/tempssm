# Extract the drift (slope) component as a time series

Extract the drift (slope) component as a time series

## Usage

``` r
get_drift_ts(
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

A univariate `ts` object of the selected drift estimate (in degrees
Celsius per year). If `ci = TRUE`, a multivariate `ts` object with
columns `drift`, `lwr`, and `upr` is returned.

## Details

The drift component is scaled to represent change per year. For example,
monthly data (frequency = 12) are multiplied by 12.

## Examples

``` r
if (FALSE) { # \dontrun{
data(niigata_sst)
res <- tempssm(niigata_sst)
drift <- get_drift_ts(res)
} # }
```
