# Log-likelihood method for tempssm objects (S3 method)

Log-likelihood method for tempssm objects (S3 method)

## Usage

``` r
# S3 method for class 'tempssm'
logLik(object, ...)
```

## Arguments

- object:

  An object of class `"tempssm"`, typically returned by
  [`tempssm()`](https://akihirao.github.io/tempssm/reference/tempssm.md).

- ...:

  Additional arguments passed to the generic
  [`logLik()`](https://rdrr.io/r/stats/logLik.html) function. These are
  currently ignored but are included for compatibility with the generic
  interface.

## Value

An object of class `"logLik"` containing the numeric log-likelihood,
with `df` and `nobs` attributes.

## See also

[`AIC.tempssm`](https://akihirao.github.io/tempssm/reference/AIC.tempssm.md)
for computing AIC,
[`get_aic`](https://akihirao.github.io/tempssm/reference/get_aic.md) as
a convenience wrapper.

## Examples

``` r
if (FALSE) { # \dontrun{
data(niigata_sst)

# fit model
res <- tempssm(niigata_sst)

# extract log-likelihood
ll <- logLik(res)

ll

# access attributes
attr(ll, "df") # number of parameters
attr(ll, "nobs") # number of observations
} # }
```
