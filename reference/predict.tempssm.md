# Forecast from a fitted tempssm model

Produces forecasts after the end of the observed sample from a fitted
`tempssm` model without exogenous variables. By default, the method
returns the forecast for the next single time point.

## Usage

``` r
# S3 method for class 'tempssm'
predict(
  object,
  n.ahead = 1L,
  interval = c("none", "confidence", "prediction"),
  level = 0.95,
  ...
)
```

## Arguments

- object:

  A fitted object of class `"tempssm"` returned by
  [`tempssm()`](https://akihirao.github.io/tempssm/reference/tempssm.md).

- n.ahead:

  Positive integer giving the forecast horizon. The default is one time
  point.

- interval:

  Character scalar specifying the interval type. One of `"none"`,
  `"confidence"`, or `"prediction"`.

- level:

  Numeric scalar between 0 and 1 specifying the confidence level used
  when `interval` is not `"none"`.

- ...:

  Additional arguments passed to the KFAS prediction method via
  [`stats::predict()`](https://rdrr.io/r/stats/predict.html).

## Value

For `interval = "none"`, a univariate `ts` object of point forecasts.
When intervals are requested, a multivariate `ts` object with columns
`fit`, `lwr`, and `upr`.

## Details

Forecasts begin at the first regular time point after the end of
`object$temp_data`. They account for uncertainty in the estimated state
conditional on the fitted model parameters, but not parameter estimation
uncertainty.

Models fitted with exogenous variables are not yet supported by this
method because future values of those variables must be supplied
explicitly. Such models produce an informative error rather than
assuming zero-valued future covariates.

This method produces forecasts beyond the end of the sample. It does not
return in-sample one-step-ahead predictions.

## Examples

``` r
if (FALSE) { # \dontrun{
data(niigata_sst)
res <- tempssm(niigata_sst)

# Forecast the next time point
predict(res)

# Forecast the next 12 time points with prediction intervals
predict(res, n.ahead = 12, interval = "prediction", level = 0.95)
} # }
```
