# Compare AIC across fitted tempssm models

Compare AIC across fitted tempssm models

## Usage

``` r
compare_tempssm_aic(models)
```

## Arguments

- models:

  A named list containing at least two fitted and converged `tempssm`
  objects. Names must be non-empty and unique and are used as model
  labels in the returned table.

## Value

A `tibble` ordered by increasing AIC with columns `model`, `logLik`,
`df`, `nobs`, `observed_n`, `start`, `end`, `frequency`, `likelihood`,
`AIC`, `delta_AIC`, and `weight`.

## Details

AIC comparisons require all models to use the same response
observations. This function verifies response length, frequency,
complete time index, values, and missing-value positions before
calculating the table. Model structure may otherwise differ, including
autoregressive order, inclusion of a seasonal component, and inclusion
or number of exogenous variables.

All models must also have been fitted using the same likelihood type.
The function uses the `marginal` setting stored in each fitted model and
does not provide an argument for changing the likelihood after fitting.
Models fitted with diffuse and marginal likelihoods are therefore
rejected when supplied together.

In the current tempssm model, changing `use_season` changes the diffuse
state structure, and adding exogenous variables adds diffuse regression
states. Changing only the AR order adds non-diffuse AR states and does
not ordinarily change the diffuse rank. When comparing models whose
seasonal or exogenous structure differs, fitting every model with
`marginal = TRUE` is recommended because the marginal likelihood
includes the correction associated with diffuse initialization.
Comparisons using `marginal = FALSE` remain available when that setting
is shared by every model.

Akaike weights are calculated as normalized relative likelihoods from
`delta_AIC`.

## See also

[`tempssm`](https://akihirao.github.io/tempssm/reference/tempssm.md),
[`AIC.tempssm`](https://akihirao.github.io/tempssm/reference/AIC.tempssm.md),
[`get_aic`](https://akihirao.github.io/tempssm/reference/get_aic.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(niigata_sst)

models <- list(
  seasonal = tempssm(
    niigata_sst,
    use_season = TRUE,
    marginal = TRUE
  ),
  no_season = tempssm(
    niigata_sst,
    use_season = FALSE,
    marginal = TRUE
  )
)

compare_tempssm_aic(models)
} # }
```
