# AIC method for tempssm objects (S3 method)

Compute the Akaike Information Criterion (AIC) for a model fitted by
[`tempssm()`](https://akihirao.github.io/tempssm/reference/tempssm.md).
This method extends the generic
[`AIC`](https://rdrr.io/r/stats/AIC.html) function.

## Usage

``` r
# S3 method for class 'tempssm'
AIC(object, ..., k = 2)
```

## Arguments

- object:

  An object of class `"tempssm"`, typically returned by
  [`tempssm()`](https://akihirao.github.io/tempssm/reference/tempssm.md).

- ...:

  Additional arguments passed to the generic
  [`AIC()`](https://rdrr.io/r/stats/AIC.html) function. These are
  currently ignored but are included for compatibility with the generic
  interface.

- k:

  Numeric penalty coefficient for the number of parameters. This
  argument is included for compatibility with
  [`AIC`](https://rdrr.io/r/stats/AIC.html) but is not used in the
  tempssm method, where the standard AIC definition (`k = 2`) is
  applied.

## Value

A numeric scalar giving the AIC of the fitted `tempssm` model.

## Examples

``` r
if (FALSE) { # \dontrun{
data(niigata_sst)
res <- tempssm(niigata_sst)
aic <- AIC(res)
} # }
```
