# Compute seasonal temperature anomalies

Temperature anomalies are calculated by subtracting a seasonal
climatology from each observation. The climatology can be computed
either from the full time series or from a user-defined baseline period.

## Usage

``` r
compute_temp_anomaly(temp_ts, baseline = NULL)
```

## Arguments

- temp_ts:

  Univariate temperature time series of class `ts`. The time series must
  have an integer frequency greater than 1. For example,
  `frequency = 12` represents monthly data, `frequency = 24` represents
  twice-monthly data, and `frequency = 4` represents four seasonal
  observations per year.

- baseline:

  Optional numeric vector of length 2 specifying the reference period
  for climatology in complete seasonal cycles (e.g., `c(1981, 2010)`).
  If `NULL`, the full period is used.

## Value

A `ts` object of seasonal temperature anomalies with the same `start`
and `frequency` as `temp_ts`.

## Details

Seasonal climatology is computed using
[`compute_monthly_climatology()`](https://akihirao.github.io/tempssm/reference/compute_monthly_climatology.md).
If `baseline` is provided, climatology is estimated using only data
within the specified reference period. Missing values are ignored when
calculating climatological means.

## Examples

``` r
temp_ts <- ts(
  rnorm(12 * 30, mean = 10),
  start = c(1981, 1),
  frequency = 12
)

# Full-period climatology
anom_all <- compute_temp_anomaly(temp_ts)

# Baseline climatology (1981-2010)
anom_base <- compute_temp_anomaly(temp_ts, baseline = c(1981, 2010))
```
