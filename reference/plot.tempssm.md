# Plot method for tempssm objects

Default base R plot method for objects fitted by
[`tempssm()`](https://akihirao.github.io/tempssm/reference/tempssm.md).
This method delegates to
[`autoplot.tempssm()`](https://akihirao.github.io/tempssm/reference/autoplot.tempssm.md)
so that `plot(res)` displays the same component plots as
`autoplot(res)`.

## Usage

``` r
# S3 method for class 'tempssm'
plot(x, ...)
```

## Arguments

- x:

  An object returned by
  [`tempssm()`](https://akihirao.github.io/tempssm/reference/tempssm.md).

- ...:

  Additional arguments passed to
  [`autoplot.tempssm()`](https://akihirao.github.io/tempssm/reference/autoplot.tempssm.md).

## Value

Invisibly returns the plotted `ggplot` object for a single component, or
the `gtable` object combining all component plots.

## Examples

``` r
if (FALSE) { # \dontrun{
data(niigata_sst)
res <- tempssm(niigata_sst)
plot(res)
plot(res, component = "level", ci = FALSE)
} # }
```
