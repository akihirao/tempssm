# Assign variable names to a ts object

`set_ts_name()` assigns variable (column) names to a time series object
of class `ts`. The function supports both univariate and multivariate
time series and ensures that variable names are handled consistently
within the tempssm framework.

For univariate series, the input is converted to a single-column `ts`
object with the specified name. For multivariate series, column names
are assigned directly.

## Usage

``` r
set_ts_name(ts_in, label, quiet = FALSE)
```

## Arguments

- ts_in:

  A time series object of class `ts`. May be univariate or multivariate.

- label:

  A character string or character vector specifying variable names. For
  a univariate series, must be a length-one character string. For a
  multivariate series, its length must equal the number of columns in
  `ts_in`.

- quiet:

  Logical scalar; if `TRUE`, suppresses the informational message.

## Value

A `ts` object with assigned variable names and the same time scale as
`ts_in`, including `start`, `end`, and `frequency`.

## Details

This function does not modify the time attributes of the input series
(e.g., start time, frequency). It only assigns variable names while
preserving the internal structure required by downstream functions such
as
[`tempssm()`](https://akihirao.github.io/tempssm/reference/tempssm.md).

## See also

[`trim_ts_overlap`](https://akihirao.github.io/tempssm/reference/trim_ts_overlap.md),
[`split_multi_ts`](https://akihirao.github.io/tempssm/reference/split_multi_ts.md),
[`tempssm`](https://akihirao.github.io/tempssm/reference/tempssm.md)

## Examples

``` r
## Univariate example
ts_uni <- ts(
  rnorm(12 * 10),
  start = c(2000, 1),
  frequency = 12
)

ts_uni_named <- set_ts_name(ts_uni, label = "temperature")
#> Assigned 1 variable name to time series

## Multivariate example
ts_multi <- ts(
  matrix(rnorm(240), ncol = 3),
  start = c(2000, 1),
  frequency = 12
)

ts_multi_named <- set_ts_name(
  ts_multi,
  label = c("precip", "solar", "wind")
)
#> Assigned 3 variable names to time series
```
