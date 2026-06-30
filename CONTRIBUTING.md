# Contributing

Thanks for your interest in contributing!

## How to contribute

- Report bugs via GitHub Issues
- Suggest features via Issues
- Submit pull requests

## Development

- Please run
  [`devtools::check()`](https://devtools.r-lib.org/reference/check.html)
- Include tests when possible

## Extended tests

The regular test suite is run with:

``` r

devtools::test()
```

Additional extended tests are available through the same `testthat`
framework. They are skipped by default and can be enabled with:

``` sh
TEMPSSM_EXTENDED_TESTS=true Rscript -e 'devtools::test()'
```

These tests currently run additional parameter-recovery checks across
multiple random seeds. They do not require downloads, large external
data sets, special platforms, or manual inspection of generated
artefacts. On a typical development machine, the current extended test
adds only a few seconds to the test run.

## Code of Conduct

This project follows a Code of Conduct.
