# Extract the Akaike Information Criterion (AIC)

Compute the Akaike Information Criterion (AIC) for a fitted `tempssm`
model using the model log-likelihood and the number of estimated
parameters.

## Usage

``` r
get_aic(res, marginal = NULL)
```

## Arguments

- res:

  An object of class `"tempssm"` returned by
  [`tempssm()`](https://akihirao.github.io/tempssm/reference/tempssm.md).

- marginal:

  Logical scalar or `NULL`. If `NULL`, use the likelihood setting stored
  when the model was fitted. Set to `TRUE` or `FALSE` to evaluate the
  fitted parameters with the marginal or diffuse likelihood,
  respectively. An explicit value does not refit the model.

## Value

A numeric scalar representing the AIC of the fitted model.

## Details

The number of parameters is determined from the optimization results
stored in the fitted model. If exogenous variables are included, their
coefficients are added to the parameter count.

## Examples

``` r
if (FALSE) { # \dontrun{
data(niigata_sst)
res <- tempssm(niigata_sst)
aic <- get_aic(res)
} # }
```
