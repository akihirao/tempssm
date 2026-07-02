# tempssm 0.2.4

* Updated package authorship metadata and related documentation to reflect the
  current contributor and co-author list.
* Added rOpenSci statistical software standards documentation and supporting
  tests for parameter recovery, edge cases, noise susceptibility, and extended
  test workflows.
* Clarified missing-value handling for temperature responses and exogenous
  covariates.
* Added a `marginal` option to select diffuse or marginal likelihood
  consistently during KFAS model fitting and likelihood evaluation.
* Added `compare_tempssm_aic()` for validated AIC comparison across fitted
  models using the same response data and likelihood type.

# tempssm 0.2.3.1

* Removed several package dependencies by replacing small helper uses with
  local implementations.
* Replaced patchwork-based `autoplot.tempssm()` composition with an internal
  gtable layout helper.


# tempssm 0.2.2

* Major update of core function tempssm().


# tempssm 0.2.1

* The package was previously named `ThermoSSM` and has been renamed to `tempssm`.
* Generalized the seasonal component to support models with or without seasonality.
