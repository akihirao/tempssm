# Southern Oscillation Index (SOI)

The Southern Oscillation Index (SOI) is a standardized index based on
the observed sea level pressure (SLP) differences between Tahiti and
Darwin, Australia. The SOI is one measure of the large-scale
fluctuations in air pressure occurring between the western and eastern
tropical Pacific (i.e., the state of the Southern Oscillation) during El
Niño and La Niña episodes.

## Usage

``` r
soi
```

## Format

A `ts` object with:

- frequency:

  12 (monthly data)

- start:

  January 1951

- end:

  December 2025

## Source

National Centers for Environmental Information, National Oceanic and
Atmospheric Administration (NAOO). Data obtained from the below website:
<https://www.ncei.noaa.gov/access/monitoring/enso/soi>

## Details

This time series represents the SOI index provided by National Centers
for Environmental Information, National Oceanic and Atmospheric
Administration (NOAA). The observation period spans from January 1951 to
December 2025.

## Examples

``` r
data(soi)
plot(soi, ylab = "SOI index")

```
