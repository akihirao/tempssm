# Convert a data frame of monthly temperature time series to a `ts` object

Convert a data frame of monthly temperature time series to a `ts` object

## Usage

``` r
convert_monthly_df_to_ts(df)
```

## Arguments

- df:

  A data frame containing monthly temperature data with columns `Date`
  (class `Date`) and `Temp` (numeric).

## Value

A univariate monthly `ts` object representing the temperature time
series, with `frequency = 12`.

## Details

This function converts a data frame containing a monthly temperature
time series into an R `ts` object with a frequency of 12.

The input data frame must contain the following columns:

- `Date`:

  A date column indicating the time index. If the exact day of the month
  is not uniquely defined, any arbitrary day (e.g., the first day of the
  month) may be used.

- `Temp`:

  Monthly temperature values. Missing values should be represented as
  `NA`.

The function assumes that the input data represent a regularly spaced
monthly time series. If missing months are detected in the `Date`
column, a warning is issued and explicit `NA` values are inserted.

The expected format is:


    Date        Temp
    2001-01-01  10.4
    2001-02-01   8.2
    2001-03-01   NA
    2001-04-01  13.6
    ...

The output is a regular monthly `ts` object. Months are represented as
equally spaced observations with `frequency = 12`; month lengths are not
converted to a fixed number of days.

## Examples

``` r
## Create a data frame of monthly temperature data
df <- data.frame(
  Date = as.Date(c(
    "2001-01-01",
    "2001-02-01",
    "2001-03-01",
    "2001-04-01",
    "2001-05-01"
  )),
  Temp = c(10.4, 8.2, NA, 13.6, 16.1)
)

## Convert to a monthly ts object
temp_ts <- convert_monthly_df_to_ts(df)
#> Converted data frame to monthly time series with 5 observations

## Inspect the result
temp_ts
#>       Jan  Feb  Mar  Apr  May
#> 2001 10.4  8.2   NA 13.6 16.1
```
