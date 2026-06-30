# Convert a daily zoo object to a monthly `ts` object

Convert a daily zoo object to a monthly `ts` object

## Usage

``` r
daily_zoo_to_monthly_ts(zoo_obj, var = "Temp", na.rm = TRUE, na_prop_max = 1)
```

## Arguments

- zoo_obj:

  A `zoo` object with daily observations. The index must be of class
  `Date` or `POSIXt`, and the object must include a numeric column named
  by `var`.

- var:

  A character string specifying the name of the variable to be
  aggregated (default: `"Temp"`).

- na.rm:

  Logical scalar; should missing values be removed before averaging?

- na_prop_max:

  Numeric scalar giving the maximum allowed proportion of NA values
  within a month. If the proportion of missing values exceeds this
  threshold, the monthly mean is set to NA. Default is `1` (no
  additional filtering).

## Value

A univariate monthly `ts` object with `frequency = 12`.

## Details

Daily observations are grouped by calendar month using
[`zoo::as.yearmon()`](https://rdrr.io/pkg/zoo/man/yearmon.html). The
resulting monthly series is represented as a regular base R `ts` object
with `frequency = 12`. Month lengths are not converted to a fixed day
count; the calendar index determines month membership. Months with no
observations between the first and last observed months are represented
explicitly as `NA`.

The missing-value proportion is calculated from records present within
each month. Unobserved calendar days are not counted as missing;
incorporating expected daily coverage may be considered in a future
version.

## Examples

``` r
if (FALSE) { # \dontrun{
sst_zoo <- get_jma_sst_zoo(sea_area_id = 138)
sst_ts_monthly <- daily_zoo_to_monthly_ts(sst_zoo)
} # }
```
