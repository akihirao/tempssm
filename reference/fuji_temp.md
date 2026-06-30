# Monthly mean air temperature at the summit of Mt. Fuji

A long-term monthly mean air temperature time series observed at the
summit of Mt. Fuji, Japan. The data are provided as a `ts` object and
can be used as an example dataset for state-space modeling of
temperature time series.

## Usage

``` r
fuji_temp
```

## Format

A `ts` object with:

- frequency:

  12 (monthly data)

- start:

  July 1937

- end:

  December 2025

## Source

Japan Meteorological Agency (JMA). Data obtained from the JMA website:
<https://www.jma.go.jp/jma/indexe.html>

## Details

The time series represents monthly mean air temperature (in degrees
Celsius) at the summit of Mt. Fuji. The observation period spans from
July 1937 to December 2025. Missing values, if any, are encoded as `NA`.

## Examples

``` r
data(fuji_temp)
plot(fuji_temp,
  ylab = "Temperature (°C)",
  main = "Monthly mean air temperature at the summit of Mt. Fuji"
)

```
