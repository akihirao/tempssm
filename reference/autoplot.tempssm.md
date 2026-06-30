# Autoplot method for tempssm objects

Create diagnostic and component plots for a model fitted by
[`tempssm()`](https://akihirao.github.io/tempssm/reference/tempssm.md).
By default, all available components are plotted.

## Usage

``` r
# S3 method for class 'tempssm'
autoplot(
  object,
  component = NULL,
  ci = TRUE,
  ci_level = 0.95,
  nrow = 2,
  ncol = 2,
  ...
)
```

## Arguments

- object:

  An object returned by
  [`tempssm()`](https://akihirao.github.io/tempssm/reference/tempssm.md).

- component:

  Character string specifying which component to plot. One of `"level"`,
  `"drift"`, `"season"`, or `"ar1"`. If NULL (default), all components
  are plotted.

- ci:

  Logical; if TRUE, pointwise confidence intervals are returned.

- ci_level:

  Numeric confidence level between 0 and 1 (default: 0.95).

- nrow, ncol:

  Positive integers specifying the layout used when all components are
  plotted. The default is a 2 by 2 layout. These arguments are ignored
  when `component` is specified.

- ...:

  Additional arguments passed to the corresponding `autoplot_*()`
  function.

## Value

A `ggplot` object if a single component is requested, or a `gtable`
object combining all component plots.

## Examples

``` r
if (FALSE) { # \dontrun{
data(niigata_sst)
res <- tempssm(niigata_sst)

# plot all components with 95% confidence interval
autoplot(res)

# plot all components without 95% confidence interval
autoplot(res, ci = FALSE)

# plot all components in one column
autoplot(res, nrow = 4, ncol = 1)

# plot each of components
autoplot(res, component = "level", ci = FALSE)
} # }
```
