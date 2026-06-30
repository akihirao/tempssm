# Residual diagnostics for tempssm models

Compute residual diagnostic statistics for a fitted tempssm model. The
output is returned as a tidy tibble suitable for meta-analysis across
many fitted models.

## Usage

``` r
diagnose_residuals(res, JB_test = FALSE)
```

## Arguments

- res:

  An object of class `"tempssm"` returned by
  [`tempssm()`](https://akihirao.github.io/tempssm/reference/tempssm.md).

- JB_test:

  Logical scalar; if TRUE, the Jarque–Bera test is included.

## Value

A `tibble` with one row containing residual diagnostic statistics,
including Ljung–Box test results and kurtosis. If `JB_test = TRUE`,
Jarque–Bera test statistics are also included.

## Examples

``` r
if (FALSE) { # \dontrun{
data(niigata_sst)

# Fit model
res <- tempssm(niigata_sst)

# Residual diagnostics (tibble output)
diag <- diagnose_residuals(res)

diag
} # }
```
