# AIC method for tempssm objects (S3 method)

Compute the Akaike Information Criterion (AIC) for a model fitted by
[`tempssm()`](https://akihirao.github.io/tempssm/reference/tempssm.md).
This method extends the generic
[`AIC`](https://rdrr.io/r/stats/AIC.html) function.

## Usage

``` r
# S3 method for class 'tempssm'
AIC(object, ..., k = 2, marginal = NULL)
```

## Arguments

- object:

  An object of class `"tempssm"`, typically returned by
  [`tempssm()`](https://akihirao.github.io/tempssm/reference/tempssm.md).

- ...:

  Additional arguments passed to the generic
  [`AIC()`](https://rdrr.io/r/stats/AIC.html) function. These are
  currently ignored.

- k:

  Numeric penalty coefficient for the number of parameters. The default
  `k = 2` gives the standard AIC definition.

- marginal:

  Logical scalar or `NULL`. If `NULL`, use the likelihood setting stored
  when the model was fitted. Set to `TRUE` or `FALSE` to evaluate the
  fitted parameters with the marginal or diffuse likelihood,
  respectively. An explicit value does not refit the model.

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
