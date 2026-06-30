# Retrieve daily mean sea-surface temperature as a monthly ts object from JMA

This function downloads publicly available daily mean sea-surface
temperature (SST) data for Japanese coastal waters provided by the Japan
Meteorological Agency (JMA), and returns the monthly average data as a
`ts` object.

## Usage

``` r
get_jma_sst_ts(sea_area_id, na_prop_max = 1)
```

## Arguments

- sea_area_id:

  Character or numeric scalar giving the JMA sea area ID (numeric values
  are accepted and internally coerced to character). For example, 138
  corresponding to the coastal sea off southern Ibaraki. A list of sea
  area IDs and their corresponding regions is available at:
  <https://www.data.jma.go.jp/kaiyou/data/db/kaikyo/series/engan/eg_areano.html>

- na_prop_max:

  Numeric scalar giving the maximum allowed proportion of NA values
  within a month. If the proportion of missing values exceeds this
  threshold, the monthly mean is set to NA. Default is `1` (no
  additional filtering).

## Value

A univariate monthly `ts` object of mean sea-surface temperature, with
`frequency = 12`.

## Details

The function retrieves a text-format dataset from the JMA website, using
the same retrieval and parsing pathway as
[`get_jma_sst_zoo()`](https://akihirao.github.io/tempssm/reference/get_jma_sst_zoo.md).
The daily data are then aggregated by
[`daily_zoo_to_monthly_ts()`](https://akihirao.github.io/tempssm/reference/daily_zoo_to_monthly_ts.md).
Missing values in the original dataset are represented as `NA`.

To comply with good API usage practices, HTTP requests sent by this
function include a custom *User-Agent* string identifying the tempssm
package. Users may override the default User-Agent by setting the
`"tempssm.user_agent"` option, for example:
`options(tempssm.user_agent = "my-analysis/1.0")`.

## Examples

``` r
if (FALSE) { # \dontrun{
sst_138_ts <- get_jma_sst_ts(sea_area_id = 138)
head(sst_138_ts)
} # }
```
