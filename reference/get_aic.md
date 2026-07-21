# Deprecated AIC helper for tempssm objects

This function is retained for backward compatibility but is deprecated.
AIC is intentionally not computed for `tempssm` objects.

## Usage

``` r
get_aic(res, marginal = NULL)
```

## Arguments

- res:

  An object of class `"tempssm"` returned by
  [`tempssm()`](https://akihirao.github.io/tempssm/reference/tempssm.md).

- marginal:

  Logical scalar or `NULL` accepted for backward compatibility. This
  argument is ignored.

## Value

This function always raises an error.

## Examples

``` r
if (FALSE) { # \dontrun{
data(niigata_sst)
res <- tempssm(niigata_sst)
logLik(res)
} # }
```
