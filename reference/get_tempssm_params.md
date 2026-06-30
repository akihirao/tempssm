# Extract estimated parameters in the fitted models

Extract estimated parameters in the fitted models

## Usage

``` r
get_tempssm_params(res)
```

## Arguments

- res:

  An object of class `"tempssm"` returned by
  [`tempssm()`](https://akihirao.github.io/tempssm/reference/tempssm.md).

## Value

A `list` object of the estimated parameters.

## Examples

``` r
if (FALSE) { # \dontrun{
data(niigata_sst)
res <- tempssm(niigata_sst)
params <- get_tempssm_params(res)
} # }
```
