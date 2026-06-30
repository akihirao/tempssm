# Compute seasonal climatology (mean seasonal cycle)

Compute seasonal climatology (mean seasonal cycle)

## Usage

``` r
compute_monthly_climatology(temp_ts)
```

## Arguments

- temp_ts:

  Univariate temperature time series of class `ts`. The time series must
  have an integer frequency greater than 1. For example,
  `frequency = 12` represents monthly data, `frequency = 24` represents
  twice-monthly data, and `frequency = 4` represents four seasonal
  observations per year.

## Value

A tibble with one row per seasonal period containing the climatological
mean temperature. For monthly data, the first column is named `Month`
for backward compatibility. For other frequencies, the first column is
named `Period`.

## Examples

``` r
temp_ts <- ts(
  rnorm(12 * 30, mean = 10),
  start = c(1981, 1),
  frequency = 12
)

monthly_climatology <- compute_monthly_climatology(temp_ts)
```
