# Read and convert a monthly temperature CSV file to a `ts` object

Read and convert a monthly temperature CSV file to a `ts` object

## Usage

``` r
read_monthly_temp_ts(csv)
```

## Arguments

- csv:

  A length-one character string giving the path to a UTF-8 CSV file with
  columns `Year`, `Month`, and `Temp`.

## Value

A univariate monthly `ts` object representing the temperature time
series, with `frequency = 12`.

## Details

The input CSV file must contain monthly temperature data with three
columns: `Year`, `Month`, and `Temp`. The expected format is:


    Year,Month,Temp
    2001,1,10.4
    2001,2,8.2
    2001,3,NA
    2001,4,13.6
    ...

The output is a regular monthly `ts` object. The `Year` and `Month`
columns define calendar months, which are represented as equally spaced
observations with `frequency = 12`. Month lengths are not converted to a
fixed number of days.

## Examples

``` r
## Create a temporary CSV file with monthly temperature data
tmp_csv <- tempfile(fileext = ".csv")

writeLines(
  c(
    "Year,Month,Temp",
    "2001,1,10.4",
    "2001,2,8.2",
    "2001,3,NA",
    "2001,4,13.6",
    "2001,5,16.1"
  ),
  tmp_csv
)

## Read the CSV file and convert to a monthly ts object
temp_ts <- read_monthly_temp_ts(tmp_csv)
#> Loaded monthly temperature series with 5 observations

## Inspect the result
temp_ts
#>       Jan  Feb  Mar  Apr  May
#> 2001 10.4  8.2   NA 13.6 16.1
```
