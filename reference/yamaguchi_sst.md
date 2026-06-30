# Montly mean temperature time series off Yamaguchi Prefecture

A monthly mean temperature time series observed in the coastal waters
off Yamaguchi Prefecture, Japan. The data are provided as a `zoo` object
and serve as an example dataset for state-space analysis of temperature
time series.

## Usage

``` r
yamaguchi_sst
```

## Format

A `zoo` object with:

- index:

  Date from 1982-01-01 to 2025-12-31

- Temp:

  Monthly mean temperature (°C)

## Source

Japan Meteorological Agency (JMA). Data obtained from the JMA website:
<https://www.jma.go.jp/jma/indexe.html>

## Details

The variable name is unified as `Temp` for consistency with other
example datasets in the package. In this dataset, `Temp` represents
monthly mean sea surface temperature. The observation period spans from
January 1, 1982 to December 31, 2025. Missing values, if any, are
encoded as `NA`.

## Examples

``` r
data(yamaguchi_sst)
plot(yamaguchi_sst,
  ylab = "Temperature (°C)",
  main = "Daily mean sea surface temperature off southern Yamaguchi Pref."
)

```
