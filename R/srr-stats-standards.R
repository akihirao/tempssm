#' srr_stats
#'
#' All of the following standards initially have `@srrstatsTODO` tags.
#' These may be moved at any time to any other locations in your code.
#' Once addressed, please modify the tag from `@srrstatsTODO` to `@srrstats`,
#' or `@srrstatsNA`, ensuring that references to every one of the following
#' standards remain somewhere within your code.
#' (These comments may be deleted at any time.)
#'
#' @srrstatsVerbose TRUE
#'
#' @srrstats {G1.0} The package lists primary references in README.md and
#' package-level documentation, including Helske (2017) for KFAS and Baba et
#' al. (2024) for the temperature time-series application that motivated the
#' package.
#'
#' @srrstats {G1.1} The "Prior Art and Scope" section of README.md documents
#' that `tempssm` is not the first implementation of linear Gaussian
#' state-space models, Kalman filtering, or Kalman smoothing in R. These
#' computations are delegated to `KFAS`. The package is a domain-focused
#' workflow layer for temperature time-series analysis, integrating model
#' construction, component extraction, uncertainty summaries, diagnostics,
#' visualization, and time-series cross-validation.
#' 
#' @srrstats {G1.2} Package-level documentation includes a "Life Cycle and
#' Development Status" section. It states that `tempssm` is under active
#' development, that the current design intentionally focuses on regular
#' temperature time series represented as base R `ts` objects with integer
#' seasonal frequencies greater than 1, and that possible future direct
#' modelling support for daily or irregular time series such as `zoo` objects
#' is outside the current design goals.
#'  
#' @srrstats {G1.3} Package-level documentation includes a "Statistical
#' Terminology" section defining the main statistical terms used by the package,
#' including linear Gaussian state-space model, latent state, level, seasonal
#' and autoregressive components, exogenous variables, Kalman filtering,
#' Kalman smoothing, standardized recursive residuals, rolling-origin
#' time-series cross-validation, MAE, and MASE.
#' 
#' @srrstats {G1.4} All exported functions are documented with `roxygen2`
#' comments in source files under `R/`. Internal helper functions are either
#' documented with roxygen comments and marked `@noRd`, or kept as local
#' implementation details. The package `NAMESPACE` and `man/*.Rd` files are
#' generated from roxygen comments, with the roxygen2 version recorded in
#' `DESCRIPTION`.
#' 
#' @srrstats {G1.4a} Internal top-level helper functions are documented with
#' standard `roxygen2` comments and marked with `@noRd` to suppress `.Rd`
#' generation. S3 methods and exported functions are documented separately.
#' Local nested functions used only within a parent function are treated as
#' implementation details of that documented parent function.
#'
#' @srrstatsTODO {G2.0} *Implement assertions on lengths of inputs, particularly through asserting that inputs expected to be single- or multi-valued are indeed so.*
#' @srrstatsTODO {G2.0a} Provide explicit secondary documentation of any expectations on lengths of inputs
#'
#' @srrstatsTODO {G2.1} *Implement assertions on types of inputs (see the initial point on nomenclature above).*
#' @srrstatsTODO {G2.1a} *Provide explicit secondary documentation of expectations on data types of all vector inputs.*
#' @srrstatsTODO {G2.2} *Appropriately prohibit or restrict submission of multivariate input to parameters expected to be univariate.*
#' @srrstatsTODO {G2.3} *For univariate character input:*
#' @srrstatsTODO {G2.3a} *Use `match.arg()` or equivalent where applicable to only permit expected values.*
#' @srrstatsTODO {G2.3b} *Either: use `tolower()` or equivalent to ensure input of character parameters is not case dependent; or explicitly document that parameters are strictly case-sensitive.*
#' @srrstatsTODO {G2.4} *Provide appropriate mechanisms to convert between different data types, potentially including:*
#' @srrstatsTODO {G2.4a} *explicit conversion to `integer` via `as.integer()`*
#' @srrstatsTODO {G2.4b} *explicit conversion to continuous via `as.numeric()`*
#' @srrstatsTODO {G2.4c} *explicit conversion to character via `as.character()` (and not `paste` or `paste0`)*
#' @srrstatsTODO {G2.4d} *explicit conversion to factor via `as.factor()`*
#' @srrstatsTODO {G2.4e} *explicit conversion from factor via `as...()` functions*
#' @srrstatsTODO {G2.5} *Where inputs are expected to be of `factor` type, secondary documentation should explicitly state whether these should be `ordered` or not, and those inputs should provide appropriate error or other routines to ensure inputs follow these expectations.* 
#' @srrstatsTODO {G2.6} *Software which accepts one-dimensional input should ensure values are appropriately pre-processed regardless of class structures.* 
#' @srrstatsTODO {G2.7} *Software should accept as input as many of the above standard tabular forms as possible, including extension to domain-specific forms.* 
#' @srrstatsTODO {G2.8} *Software should provide appropriate conversion or dispatch routines as part of initial pre-processing to ensure that all other sub-functions of a package receive inputs of a single defined class or type.*
#' @srrstatsTODO {G2.9} *Software should issue diagnostic messages for type conversion in which information is lost (such as conversion of variables from factor to character; standardisation of variable names; or removal of meta-data such as those associated with [`sf`-format](https://r-spatial.github.io/sf/) data) or added (such as insertion of variable or column names where none were provided).* 
#' @srrstatsTODO {G2.10} *Software should ensure that extraction or filtering of single columns from tabular inputs should not presume any particular default behaviour, and should ensure all column-extraction operations behave consistently regardless of the class of tabular data used as input.* 
#' @srrstatsTODO {G2.11} *Software should ensure that `data.frame`-like tabular objects which have columns which do not themselves have standard class attributes (typically, `vector`) are appropriately processed, and do not error without reason. This behaviour should be tested. Again, columns created by the [`units` package](https://github.com/r-quantities/units/) provide a good test case.*
#' @srrstatsTODO {G2.12} *Software should ensure that `data.frame`-like tabular objects which have list columns should ensure that those columns are appropriately pre-processed either through being removed, converted to equivalent vector columns where appropriate, or some other appropriate treatment such as an informative error. This behaviour should be tested.* 
#' @srrstatsTODO {G2.13} *Statistical Software should implement appropriate checks for missing data as part of initial pre-processing prior to passing data to analytic algorithms.*
#' @srrstatsTODO {G2.14} *Where possible, all functions should provide options for users to specify how to handle missing (`NA`) data, with options minimally including:*
#' @srrstatsTODO {G2.14a} *error on missing data*
#' @srrstatsTODO {G2.14b} *ignore missing data with default warnings or messages issued*
#' @srrstatsTODO {G2.14c} *replace missing data with appropriately imputed values*
#' @srrstatsTODO {G2.15} *Functions should never assume non-missingness, and should never pass data with potential missing values to any base routines with default `na.rm = FALSE`-type parameters (such as [`mean()`](https://stat.ethz.ch/R-manual/R-devel/library/base/html/mean.html), [`sd()`](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/sd.html) or [`cor()`](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/cor.html)).*
#' @srrstatsTODO {G2.16} *All functions should also provide options to handle undefined values (e.g., `NaN`, `Inf` and `-Inf`), including potentially ignoring or removing such values.* 
#' @srrstatsTODO {G3.0} *Statistical software should never compare floating point numbers for equality. All numeric equality comparisons should either ensure that they are made between integers, or use appropriate tolerances for approximate equality.* 
#' @srrstatsTODO {G3.1} *Statistical software which relies on covariance calculations should enable users to choose between different algorithms for calculating covariances, and should not rely solely on covariances from the `stats::cov` function.*
#' @srrstatsTODO {G3.1a} *The ability to use arbitrarily specified covariance methods should be documented (typically in examples or vignettes).* 
#' @srrstatsTODO {G4.0} *Statistical Software which enables outputs to be written to local files should parse parameters specifying file names to ensure appropriate file suffixes are automatically generated where not provided.* 
#' @srrstatsTODO {G5.0} *Where applicable or practicable, tests should use standard data sets with known properties (for example, the [NIST Standard Reference Datasets](https://www.itl.nist.gov/div898/strd/), or data sets provided by other widely-used R packages).*
#' @srrstatsTODO {G5.1} *Data sets created within, and used to test, a package should be exported (or otherwise made generally available) so that users can confirm tests and run examples.* 
#' @srrstatsTODO {G5.2} *Appropriate error and warning behaviour of all functions should be explicitly demonstrated through tests. In particular,*
#' @srrstatsTODO {G5.2a} *Every message produced within R code by `stop()`, `warning()`, `message()`, or equivalent should be unique*
#' @srrstatsTODO {G5.2b} *Explicit tests should demonstrate conditions which trigger every one of those messages, and should compare the result with expected values.*
#' @srrstatsTODO {G5.3} *For functions which are expected to return objects containing no missing (`NA`) or undefined (`NaN`, `Inf`) values, the absence of any such values in return objects should be explicitly tested.* 
#' @srrstatsTODO {G5.4} **Correctness tests** *to test that statistical algorithms produce expected results to some fixed test data sets (potentially through comparisons using binding frameworks such as [RStata](https://github.com/lbraglia/RStata)).*
#' @srrstatsTODO {G5.4a} *For new methods, it can be difficult to separate out correctness of the method from the correctness of the implementation, as there may not be reference for comparison. In this case, testing may be implemented against simple, trivial cases or against multiple implementations such as an initial R implementation compared with results from a C/C++ implementation.*
#' @srrstatsTODO {G5.4b} *For new implementations of existing methods, correctness tests should include tests against previous implementations. Such testing may explicitly call those implementations in testing, preferably from fixed-versions of other software, or use stored outputs from those where that is not possible.*
#' @srrstatsTODO {G5.4c} *Where applicable, stored values may be drawn from published paper outputs when applicable and where code from original implementations is not available*
#' @srrstatsTODO {G5.5} *Correctness tests should be run with a fixed random seed*
#' @srrstatsTODO {G5.6} **Parameter recovery tests** *to test that the implementation produce expected results given data with known properties. For instance, a linear regression algorithm should return expected coefficient values for a simulated data set generated from a linear model.*
#' @srrstatsTODO {G5.6a} *Parameter recovery tests should generally be expected to succeed within a defined tolerance rather than recovering exact values.*
#' @srrstatsTODO {G5.6b} *Parameter recovery tests should be run with multiple random seeds when either data simulation or the algorithm contains a random component. (When long-running, such tests may be part of an extended, rather than regular, test suite; see G5.10-4.12, below).* 
#' @srrstatsTODO {G5.7} **Algorithm performance tests** *to test that implementation performs as expected as properties of data change. For instance, a test may show that parameters approach correct estimates within tolerance as data size increases, or that convergence times decrease for higher convergence thresholds.*
#' @srrstatsTODO {G5.8} **Edge condition tests** *to test that these conditions produce expected behaviour such as clear warnings or errors when confronted with data with extreme properties including but not limited to:*
#' @srrstatsTODO {G5.8a} *Zero-length data*
#' @srrstatsTODO {G5.8b} *Data of unsupported types (e.g., character or complex numbers in for functions designed only for numeric data)*
#' @srrstatsTODO {G5.8c} *Data with all-`NA` fields or columns or all identical fields or columns*
#' @srrstatsTODO {G5.8d} *Data outside the scope of the algorithm (for example, data with more fields (columns) than observations (rows) for some regression algorithms)*
#' @srrstatsTODO {G5.9} **Noise susceptibility tests** *Packages should test for expected stochastic behaviour, such as through the following conditions:*
#' @srrstatsTODO {G5.9a} *Adding trivial noise (for example, at the scale of `.Machine$double.eps`) to data does not meaningfully change results*
#' @srrstatsTODO {G5.9b} *Running under different random seeds or initial conditions does not meaningfully change results* 
#' @srrstatsTODO {G5.10} *Extended tests should included and run under a common framework with other tests but be switched on by flags such as as a `<MYPKG>_EXTENDED_TESTS="true"` environment variable.* - The extended tests can be then run automatically by GitHub Actions for example by adding the following to the `env` section of the workflow: 
#' @srrstatsTODO {G5.11} *Where extended tests require large data sets or other assets, these should be provided for downloading and fetched as part of the testing workflow.*
#' @srrstatsTODO {G5.11a} *When any downloads of additional data necessary for extended tests fail, the tests themselves should not fail, rather be skipped and implicitly succeed with an appropriate diagnostic message.*
#' @srrstatsTODO {G5.12} *Any conditions necessary to run extended tests such as platform requirements, memory, expected runtime, and artefacts produced that may need manual inspection, should be described in developer documentation such as a `CONTRIBUTING.md` or `tests/README.md` file.*
#' 
#' @srrstats {TS1.0} The package uses base R `ts` objects as its explicit
#' time-series class for modelling and related utilities. Core functions such
#' as `tempssm()`, `ts_train_test_split()`, `trim_ts_overlap()`,
#' `split_multi_ts()`, `compute_monthly_climatology()`,
#' `compute_temp_anomaly()`, and `plot_temp_dev()` validate that inputs inherit
#' from class `ts` and reject generic non-time-series inputs. Data-frame and
#' CSV inputs are accepted only by dedicated conversion utilities, which return
#' explicit `ts` objects and are not used as generic substitutes for
#' time-series inputs in modelling functions.
#' 
#' @srrstats {TS1.1} Function documentation explicitly states the expected
#' input types and classes for time-series data and conversion utilities. Core
#' modelling and time-series functions document base R `ts` inputs, including
#' `temp_data`, `exo_data`, and `multi_ts`. Conversion helpers document
#' non-time-series inputs such as data frames and CSV file paths, and document
#' their return values as explicit monthly `ts` objects. Functions for JMA SST
#' data document `zoo` inputs or outputs where daily indexed data are used, and
#' monthly `ts` outputs where data are aggregated.
#' 
#' @srrstats {TS1.2} The package implements explicit validation routines for
#' acceptable time-series classes. The core modelling path uses
#' `.tempssm_check_temp_ts()` to require a univariate base R `ts` object for
#' temperature data and `.tempssm_check_exo_ts()` to require aligned base R
#' `ts` objects for exogenous variables. Time-series utilities such as
#' `trim_ts_overlap()`, `split_multi_ts()`, `compute_monthly_climatology()`,
#' `compute_temp_anomaly()`, and `plot_temp_dev()` also reject non-`ts`
#' inputs. The core modelling path and seasonal climatology helpers accept
#' integer seasonal frequencies greater than 1, not only monthly
#' `frequency = 12` data. Daily SST conversion routines that operate on
#' irregular daily data explicitly require `zoo` inputs before aggregating them
#' to monthly `ts` objects. These validation paths are covered by unit tests
#' for valid and invalid class inputs.
#' 
#' @srrstats {TS1.3} Core model inputs are passed through the single internal
#' pre-processing routine `.tempssm_prepare_model_inputs()`. This routine
#' validates temperature and optional exogenous inputs, standardizes unnamed
#' exogenous variables when allowed by the calling workflow, and returns a
#' uniform list containing validated base R `ts` objects, the common
#' frequency, and the number of observations. The primary model fitting
#' function `tempssm()` and the cross-validation splitter
#' `ts_train_test_split()` both use this routine before passing data to
#' downstream model construction, fitting, or fold-generation code. Conversion
#' helpers for data-frame, CSV, and `zoo` inputs transform those inputs to
#' explicit monthly `ts` objects before they enter the modelling path.
#' 
#' @srrstats {TS1.4} The core pre-processing routine preserves the time-based
#' attributes of accepted `ts` inputs. `.tempssm_prepare_model_inputs()`
#' returns validated `ts` objects without converting them to non-time-series
#' containers, and preserves `start`, `end`, `frequency`, and `time()` values
#' for both temperature and exogenous series, including the branch where
#' default exogenous variable names are assigned. Conversion utilities that
#' start from data-frame, CSV, or `zoo` inputs explicitly construct monthly
#' `ts` outputs with defined start times and `frequency = 12` before those
#' objects enter the modelling path. These behaviours are covered by unit
#' tests for the pre-processing and conversion utilities.
#' 
#' @srrstats {TS1.5} The package checks strict ordering of time indices before
#' model fitting and cross-validation. The core input path calls
#' `.tempssm_check_ts_order()` from `.tempssm_check_temp_ts()` and
#' `.tempssm_check_exo_ts()` to require strictly increasing `time()` values for
#' accepted `ts` inputs. Conversion helpers also validate calendar order before
#' constructing monthly `ts` objects: data-frame input is sorted by `Date` and
#' duplicate months are rejected, while CSV input must have strictly
#' increasing year-month rows with no duplicate year-month combinations.
#' Missing months are warned about separately because they represent
#' regularity/completeness rather than ordering.
#' 
#' @srrstats {TS1.6} Ordering violations are caught during input
#' pre-processing rather than during model fitting. Core `ts` inputs are
#' checked by `.tempssm_check_ts_order()` before they are passed to model
#' construction or cross-validation fold generation. Data-frame conversion
#' catches unordered `Date` values with a warning before sorting, and rejects
#' duplicate monthly indices. CSV conversion rejects duplicate or non-increasing
#' year-month rows before constructing a `ts` object. `zoo` conversion rejects
#' missing or non-increasing indices before monthly aggregation. Unit tests
#' cover these ordering checks in the corresponding pre-processing functions.
#' 
#' @srrstats {TS1.7} The package accepts vector inputs carrying class
#' `units` from the `units` package without adding `units` as a hard runtime
#' dependency. Core model input pre-processing detects `units` attached to
#' `ts` objects, converts values to numeric vectors for downstream state-space
#' modelling, and preserves the original `ts` time attributes. Data-frame and
#' `zoo` conversion helpers also strip `units` from value columns before
#' constructing monthly `ts` objects. The conversion emits an explicit warning
#' so users know that unit metadata are not retained in model inputs. Tests for
#' this behaviour are conditional on the optional `units` package.
#' 
#' @srrstats {TS1.8} Package-level documentation includes a "Time Index and
#' Calendar Conventions" section. Core modelling functions use base R `ts`
#' objects, where `frequency` defines the number of observations per seasonal
#' cycle and is preserved by the core modelling path. Frequencies such as 4,
#' 12, 24, and 36 are treated as regularly spaced observations within a
#' seasonal cycle and are not converted to fixed numbers of days. Daily `zoo`
#' inputs are indexed by `Date` or `POSIXt` and grouped into calendar months
#' with `zoo::as.yearmon()` before aggregation to monthly `ts` objects.
#' Cross-validation window parameters are documented and used as observation
#' counts rather than calendar durations.
#' 
#' @srrstats {TS2.0} Core modelling functions require regular base R `ts`
#' inputs, where missing observations can only be represented explicitly as
#' `NA` values. Conversion utilities that construct monthly `ts` objects from
#' data-frame, CSV, or daily `zoo` inputs detect implicit missing months before
#' constructing the output. Missing months between the first and last observed
#' month are converted to explicit `NA` values with a diagnostic warning, so
#' irregular input is not silently collapsed into a shorter regular `ts`
#' object. Unit tests cover these conversion paths.
#' 
#' @srrstats {TS2.1} The main modelling and cross-validation entry points
#' provide a \code{na_action} argument controlling how explicit missing values
#' in regular `ts` inputs are handled. The argument is passed through the
#' shared input pre-processing routine `.tempssm_prepare_model_inputs()`, so
#' `tempssm()` and `ts_train_test_split()` use the same missing-data policy.
#' Conversion utilities also expose missing-data controls where aggregation is
#' performed, such as `na.rm` and `na_prop_max` in
#' `daily_zoo_to_monthly_ts()`.
#' 
#' @srrstats {TS2.1a} Setting `na_action = "error"` in `tempssm()` or
#' `ts_train_test_split()` stops during input pre-processing if explicit
#' missing values are detected in `temp_data` or `exo_data`.
#' 
#' @srrstats {TS2.1b} The default `na_action = "warn"` issues a diagnostic
#' warning and proceeds with explicit `NA` observations. Setting
#' `na_action = "allow"` proceeds silently. Because these paths preserve the
#' original regular `ts` object and its explicit missing values, results are
#' based on the same time index rather than on an implicitly shortened series.
#' 
#' @srrstatsTODO {TS2.1c} *replace missing data with appropriately imputed values.* 
#' 
#' @srrstats {TS2.2} Stationarity of relevant moments is documented in the
#' package-level help under "Stationarity and Moment Assumptions". The observed
#' temperature series is not required to be stationary in mean because
#' first-order non-stationarity is represented explicitly through latent level,
#' drift, and seasonal components. Stationarity constraints are instead applied
#' to the autoregressive component that represents residual short-term serial
#' dependence after those components have been accounted for. Second-order
#' assumptions are represented by positive, time-invariant Gaussian disturbance
#' variances in the observation and state equations; time-varying volatility
#' and higher-order stationarity checks are outside the current package scope.
#' 
#' @srrstats {TS2.3} Stationarity assumptions and requirements are documented
#' in the package-level help and in the `tempssm()` function help. The observed
#' input series is not required to be stationary in mean; non-stationary mean
#' structure is modelled through level, drift, seasonal, and optional exogenous
#' components. Stationarity is required only for the autoregressive component,
#' and this requirement is enforced by transforming AR coefficients with
#' `KFAS::artransform()` during fitting. The documentation also states the
#' second-order assumption of positive, time-invariant Gaussian disturbance
#' variances and the current absence of automatic stationarity testing for the
#' observed input series.
#' 
#' @srrstats {TS2.4} The relevant stationarity requirement in `tempssm` is
#' the stationarity of the autoregressive component. The package implements an
#' internal AR-root check after transforming unconstrained optimization
#' parameters to AR coefficients. The same transformation-and-check helper is
#' used during model fitting and when extracting or summarising fitted AR
#' coefficients. Unit tests cover stationary and non-stationary AR
#' coefficients and confirm that transformed AR parameters are stationary.
#' 
#' @srrstats {TS2.4b} AR coefficients are transformed with
#' `KFAS::artransform()` and then checked by evaluating whether the AR
#' polynomial roots lie outside the unit circle within numerical tolerance.
#' 
#' @srrstats {TS3.0} Unit tests demonstrate that forecast uncertainty widens
#' with forecast horizon for a fitted temperature state-space model. In
#' `tests/testthat/test-predict_no_exo.R`, the test "prediction intervals widen
#' with forecast horizon" fits `tempssm()` to `niigata_sst`, obtains 24-step
#' ahead prediction intervals from the KFAS backend, and checks that interval
#' widths are non-decreasing and wider at horizon 24 than at horizon 1.
#' 
#' @srrstats {TS3.1} The same test file includes a negative-control test,
#' "prediction interval widening check detects violations", which supplies
#' decreasing and constant synthetic interval-width sequences to the widening
#' check and confirms that they are rejected. This demonstrates that the
#' TS3.0 test condition would detect forecast intervals that fail to widen
#' appropriately with horizon.
#' 
#' @srrstats {TS3.2} Package-level documentation includes a "Forecast
#' Horizons and Error Margins" section. It explains that forecast uncertainty
#' generally increases with horizon because future observations depend on
#' accumulated latent-state and observation uncertainty. The TS3.0 test uses
#' the `niigata_sst` data to demonstrate this widening, and the TS3.1
#' negative-control test demonstrates that decreasing or constant interval
#' widths would be detected as violations.
#' 
#' @srrstats {TS3.3} The package supports trimming forecast values based on
#' prediction interval width through the exported post-processing function
#' `trim_prediction_intervals()`.
#' 
#' @srrstats {TS3.3a} The package-level documentation and the examples for
#' `trim_prediction_intervals()` document how to obtain KFAS prediction
#' intervals from a fitted `tempssm` model and trim them to a maximum allowed
#' interval width.
#' 
#' @srrstats {TS3.3b} The exported function `trim_prediction_intervals()`
#' implements an explicit post-processing mechanism for trimming forecasts to
#' the longest leading horizon whose prediction interval width `upr - lwr`
#' does not exceed a user-specified threshold.
#' 
#' @srrstats {TS4.0} Package-level documentation includes a "Return Value
#' Conventions" section. It explains that the main modelling function returns
#' a classed `tempssm` object, because the result is a fitted statistical model
#' rather than a transformed time series. Accessor functions return `ts`
#' objects for fitted components, and diagnostic or summary functions return
#' documented structured objects appropriate to their purpose.
#' 
#' @srrstats {TS4.0b} The primary modelling return value is a unique S3 class,
#' `tempssm`, with methods for printing, summaries, information criteria, and
#' visualization. The `summary()` method returns a `summary.tempssm` object,
#' and tests verify the expected S3 classes of these primary return values.
#' 
#' @srrstats {TS4.2} Exported functions document the type and class of their
#' return values in roxygen `@return` fields. Structured return values also
#' document their named list elements or columns, including `tempssm`,
#' `summary.tempssm`, `logLik`, `ts`, `zoo`, `ggplot`, `gtable`, `tibble`,
#' `data.frame`, and named list outputs used by accessors, diagnostics, and
#' cross-validation helpers. Package-level documentation summarizes these
#' return-value conventions.
#' 
#' @srrstats {TS4.3} Return values that represent time series explicitly
#' include their time scales. Component accessors and time-series utilities
#' return base R `ts` objects with preserved or explicitly constructed
#' `start`, `end`, `frequency`, and `time()` attributes. Daily SST retrieval
#' returns `zoo` objects indexed by `Date`, and daily-to-monthly conversion
#' returns monthly `ts` objects with `frequency = 12`. Unit attributes from
#' optional `units` inputs are intentionally not propagated, as documented
#' under TS4.1 and in package-level documentation.
#' 
#' @srrstats {TS4.4} Package-level documentation includes a "Forecast Scale
#' and Transformations" section. It states that forecast values are returned
#' on the original numeric temperature scale because the response series is not
#' log-transformed, standardized, differenced, or otherwise transformed before
#' forecasting. It also documents that internal transformations of variance
#' and AR parameters affect model dynamics and uncertainty estimates but do not
#' change the scale of returned forecast values. The post-processing helper
#' `trim_prediction_intervals()` is documented as truncating forecast horizons
#' without modifying forecast means or interval bounds.
#' 
#' @srrstats {TS4.5} Forecast values are returned on the same numeric scale as
#' the input temperature series, so no back-transformation is required to make
#' them commensurate with the original non-stationary input data. The package
#' documentation explicitly states this design and documents forecast
#' limitations.
#' 
#' @srrstats {TS4.5c} Package-level documentation states that forecasts are
#' conditional on the selected state-space model structure, estimated
#' parameters, and any supplied future exogenous variables. It also states that
#' forecasts should not be interpreted as physically constrained forecasts
#' beyond the assumptions of the fitted model.
#' 
#' @srrstats {TS4.6} Forecast uncertainty is available through the fitted
#' `KFAS` model returned by `tempssm()`. Package-level documentation shows the
#' use of `stats::predict(res$model, n.ahead = h,
#' interval = "prediction")`, which returns point forecasts and prediction
#' interval bounds.
#' 
#' @srrstats {TS4.6c} Prediction intervals from the `KFAS` backend provide a
#' general indication of forecast error through lower and upper interval
#' bounds. The exported helper `trim_prediction_intervals()` also operates on
#' these interval bounds and documents their expected `fit`, `lwr`, and `upr`
#' structure.
#' 
#' @srrstats {TS4.7} Forecast values and observed values are clearly
#' distinguished in package outputs. Direct forecasts obtained via
#' `stats::predict()` return forecast values alone, while cross-validation
#' results store observed test values and forecasts in separate named list
#' elements.
#' 
#' @srrstats {TS4.7a} Direct forecasting uses `stats::predict()` on the fitted
#' `KFAS` model returned by `tempssm()`, and the returned object contains
#' forecast values and optional interval bounds rather than observed input
#' values.
#' 
#' @srrstats {TS4.7b} Cross-validation results returned by
#' `ts_cv_run_fold()` and `ts_cv_run()` use distinct named list elements:
#' `y_test` for observed held-out values and `y_pred` for model forecasts.
#' The `@return` documentation for `ts_cv_run_fold()` describes these elements.
#' 
#' @srrstats {TS5.0} The package implements a default S3 `plot.tempssm()`
#' method for fitted `tempssm` objects. The method delegates to
#' `autoplot.tempssm()`, so `plot(res)` displays the same component plots as
#' `autoplot(res)` and accepts the same plotting arguments.
#' 
#' @srrstats {TS5.1} Temporal plots label the horizontal axis as
#' `Time (year)`.
#' Component plots generated by `autoplot_level()`, `autoplot_drift()`,
#' `autoplot_season()`, `autoplot_ar1()`, `autoplot.tempssm()`, and
#' `plot.tempssm()` map the continuous `time()` index of the fitted `ts`
#' series to the x-axis. `plot_temp_dev()` also labels its horizontal axis as
#' `Time (year)`. Unit tests check these axis labels.
#' 
#' @srrstats {TS5.2} Time is placed on the horizontal axis by default.
#' Component plots construct data frames with a `time` column derived from
#' `time(<ts>)` and map that column to the ggplot2 x aesthetic. The default
#' `plot.tempssm()` method delegates to `autoplot.tempssm()`, so the same
#' horizontal time-axis convention is used by both plotting interfaces.
#' `plot_temp_dev()` maps the anomaly time index to the horizontal axis,
#' with that axis labelled as `Time (year)`.
#' 
#' @srrstats {TS5.3} Temporal plots print the time unit on the horizontal
#' axis by default. Component plots and `plot_temp_dev()` use
#' `Time (year)`, reflecting the continuous year scale returned by
#' `time(<ts>)` for the package's seasonal temperature time series.
#' 
#' @srrstats {TS5.5} `plot_temp_dev()` provides a `connect_missing` option
#' controlling how missing values are drawn. The default
#' `connect_missing = FALSE` preserves explicit missing values in the plotting
#' data, so line segments are broken at gaps. Setting
#' `connect_missing = TRUE` removes missing anomaly values before plotting,
#' connecting line segments across gaps. Unit tests cover both behaviours.
#'
#' @noRd
NULL

#' NA_standards
#'
#' Any non-applicable standards can have their tags changed from `@srrstatsTODO`
#' to `@srrstatsNA`, and placed together in this block, along with explanations
#' for why each of these standards have been deemed not applicable.
#' (These comments may also be deleted at any time.)
#'
#' @srrstatsNA {G1.5} The package does not currently make performance claims
#' in associated publications, such as claims of superior speed, accuracy, or
#' predictive performance relative to other software. Although the initial
#' implementation was adapted from supplementary code associated with Baba et
#' al. (2024), `tempssm` does not claim to reproduce all numerical results from
#' that publication. Examples, tests, and vignettes demonstrate package
#' functionality rather than publication-level performance claims.
#'
#' @srrstatsNA {G1.6} The package does not currently make comparative
#' performance claims against alternative implementations in other R packages.
#' `tempssm` uses `KFAS` as its computational backend and is presented as a
#' domain-focused workflow layer for temperature time-series analysis rather
#' than as a faster, more accurate, or otherwise superior alternative to
#' packages such as `KFAS`, `forecast`, `dlm`, `bsts`, or `MARSS`.
#'
#' @srrstatsNA {TS2.4a} The package uses parameter transformations to enforce
#' the relevant stationarity constraint, so routine warnings are not needed for
#' valid user inputs. An internal error is raised only if transformed AR
#' coefficients fail the stationarity check, which should indicate an internal
#' numerical or implementation problem rather than an expected user action.
#'
#' @srrstatsNA {TS2.5} The package does not construct or return
#' auto-covariance matrices whose row and column ordering would need to be
#' linked to the underlying time-series index. Core modelling functions operate
#' on regular `ts` inputs and validate that temperature and exogenous series
#' have matching time indices before fitting, while returned component series
#' preserve the `ts` time attributes of the fitted model.
#'
#' @srrstatsNA {TS2.6} The package does not construct or return
#' auto-covariance matrices, so there are no auto-covariance matrix units to
#' specify. General handling of optional `units` objects in user inputs is
#' documented separately under TS1.7 and tested in the input pre-processing and
#' conversion utilities.
#'
#' @srrstatsNA {TS4.0a} The primary output of `tempssm()` is a fitted model
#' object rather than a transformed version of the input time series, so it is
#' intentionally not returned in the same class as the input `ts` object.
#' Component accessor functions return `ts` objects where the output represents
#' a time series rather than a model object.
#'
#' @srrstatsNA {TS4.1} Unit attributes from optional `units` inputs are
#' intentionally not propagated to return values. As documented under TS1.7
#' and in package-level documentation, such inputs are converted to numeric
#' values with an explicit warning before model fitting or conversion to `ts`.
#' This avoids adding `units` as a hard runtime dependency and keeps the
#' internal state-space calculations numeric. Returned time-series objects
#' preserve time attributes but not unit attributes.
#'
#' @srrstatsNA {TS4.5a} The package does not transform the response series
#' before forecasting, so there is no transformed forecast scale for which a
#' back-transformation routine or option would be required.
#'
#' @srrstatsNA {TS4.5b} The package does not transform the response series
#' before forecasting, so there is no required user workflow for
#' back-transforming forecast values to the original input scale.
#'
#' @srrstatsNA {TS4.6a} The package does not construct a separate forecast
#' distribution object. Forecast uncertainty is exposed through prediction
#' intervals returned by the `KFAS` backend.
#'
#' @srrstatsNA {TS4.6b} The package does not return forecast standard errors
#' directly as a separate second-order moment interface. Forecast uncertainty
#' is instead represented by prediction interval bounds from the `KFAS`
#' backend.
#'
#' @srrstatsNA {TS4.7c} The package does not combine observed and forecast
#' values into a single returned table for forecasting workflows. It instead
#' returns forecasts alone for direct prediction and distinct `y_test` and
#' `y_pred` list elements for cross-validation workflows.
#'
#' @srrstatsNA {TS5.4} The package does not implement frequency-domain
#' visualizations such as spectra, periodograms, or Fourier-domain plots.
#' The `frequency` attribute of `ts` inputs is used to define seasonal
#' periodicity, but no plot uses an angular-frequency abscissa.
#'
#' @srrstatsNA {TS5.6} The package does not currently implement a forecast
#' plotting interface. Direct forecasts are obtained from
#' `stats::predict(res$model, ..., interval = "prediction")`, which returns
#' point forecasts and lower and upper prediction interval bounds. Component
#' plots show fitted latent states and their confidence intervals, but they
#' are not forecast plots.
#'
#' @srrstatsNA {TS5.7} The package does not currently implement a forecast
#' plotting interface. Direct forecasts are returned by `stats::predict()`,
#' and cross-validation helpers return observed held-out values and forecasts
#' as separate `y_test` and `y_pred` elements. Because no default forecast
#' plot is produced, there is no plot in which model input values and forecast
#' output values could be included by default.
#'
#' @srrstatsNA {TS5.8} The package does not currently implement a forecast
#' plotting interface. Direct forecasts are returned by `stats::predict()`,
#' and cross-validation outputs distinguish observed and forecast values
#' structurally through separate `y_test` and `y_pred` elements. Because no
#' default forecast plot is produced, there is no plotted visual distinction
#' to define between model input values and forecast output values.
#'
#' @noRd
NULL
