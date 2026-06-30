# Extract the smoothed drift (slope) component as a time series

Extract the smoothed drift (slope) component as a time series

## Usage

``` r
get_drift_ts(res, ci = FALSE, ci_level = 0.95)
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

A univariate `ts` object of the smoothed drift component (in degrees
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
