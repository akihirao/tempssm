# Retrieve daily mean sea-surface temperature as a zoo object from JMA

This function downloads publicly available daily mean sea-surface
temperature (SST) data for Japanese coastal waters provided by the Japan
Meteorological Agency (JMA), and returns the data as a `zoo` object
indexed by date.

## Usage

``` r
get_jma_sst_zoo(sea_area_id)
```

## Arguments

- sea_area_id:

  Character or numeric scalar giving the JMA sea area ID (numeric values
  are accepted and internally coerced to character). For example, 138
  corresponding to the coastal sea off southern Ibaraki. A list of sea
  area IDs and their corresponding regions is available at:
  <https://www.data.jma.go.jp/kaiyou/data/db/kaikyo/series/engan/eg_areano.html>

## Value

A `zoo` object of daily mean sea-surface temperature with a `Date` index
and a single numeric column named `Temp`.

## Details

The function retrieves a text-format dataset from the JMA website,
parses daily observations, and constructs a `zoo` object with calendar
dates as the index. Missing values in the original dataset are
represented as `NA`.

To comply with good API usage practices, HTTP requests sent by this
function include a custom *User-Agent* string identifying the tempssm
package. Users may override the default User-Agent by setting the
`"tempssm.user_agent"` option, for example:
`options(tempssm.user_agent = "my-analysis/1.0")`.

## Examples

``` r
if (FALSE) { # \dontrun{
sst_138_zoo <- get_jma_sst_zoo(sea_area_id = 138)
head(sst_138_zoo)
} # }
```
