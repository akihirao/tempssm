# Summary method for tempssm objects

Provides a concise summary of a fitted linear Gaussian state-space model
estimated by
[`tempssm()`](https://akihirao.github.io/tempssm/reference/tempssm.md).

## Usage

``` r
# S3 method for class 'tempssm'
summary(object, ...)
```

## Arguments

- object:

  An object of class `"tempssm"` returned by
  [`tempssm()`](https://akihirao.github.io/tempssm/reference/tempssm.md).

- ...:

  Additional arguments (currently not used).

## Value

An object of class `"summary.tempssm"`, implemented as a named list with
components `call`, `logLik`, `k`, `AIC`, `convergence`, `variances`,
`coef_ar`, `exogenous`, and `exogenous_coef`.

## Examples

``` r
if (FALSE) { # \dontrun{
data(niigata_sst)

# fit model
res <- tempssm(niigata_sst)

# compute summary
s <- summary(res)

s

# access components programmatically
s$AIC
s$variances
} # }
```
