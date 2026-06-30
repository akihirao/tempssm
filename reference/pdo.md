# Pacific Decadal Oscillation (PDO) index (JMA)

The PDO index is defined as the projection coefficient of monthly mean
sea surface temperature anomalies onto the leading empirical orthogonal
function (EOF) mode of sea surface temperature variability in the North
Pacific north of 20°N, based on data from 1901 to 2000. The anomalies
are calculated relative to the 1901–2000 mean. Prior to the EOF
analysis, the global mean sea surface temperature anomaly is removed
from the local anomaly at each grid point to reduce the influence of
global warming.

## Usage

``` r
pdo
```

## Format

A `ts` object with:

- frequency:

  12 (monthly data)

- start:

  January 1901

- end:

  December 2025

## Source

Japan Meteorological Agency (JMA). Data obtained from the JMA website:
<https://www.data.jma.go.jp/kaiyou/data/shindan/b_1/pdo/pdo.txt>

## Details

This time series represents the PDO index provided by the Japan
Meteorological Agency (JMA). The observation period spans from January
1901 to December 2025.

## Examples

``` r
data(pdo)
plot(pdo, ylab = "PDO index")

```
