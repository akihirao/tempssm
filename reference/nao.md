# North Atlantic Oscillation: NAO (Hurrell)

The winter (December thru March) station-based index of the NAO is based
on the difference of normalized sea level pressure (SLP) between Lisbon,
Portugal and Stykkisholmur/Reykjavik, Iceland since 1864. Positive
values of the NAO index are typically associated with
stronger-than-average westerlies over the middle latitudes, more intense
weather systems over the North Atlantic and wetter/milder weather over
western Europe. Monthly, seasonal and annual indices using slightly
different data sources for the southern station are also available.

## Usage

``` r
nao
```

## Format

A `ts` object with:

- frequency:

  12 (monthly data)

- start:

  January 1865

- end:

  June 2023

## Source

NSF NCAR. Hurrell North Atlantic Oscillation (NAO) Index
(station-based):
<https://climatedataguide.ucar.edu/climate-data/hurrell-north-atlantic-oscillation-nao-index-station-based>

## Details

This time series represents the NAO index provided by the Climate
Analysis Section, NCAR, Boulder, USA, Hurrell (2003). The observation
period spans from January 1865 to June 2023.

## Examples

``` r
data(nao)
plot(nao)

```
