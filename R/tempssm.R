#' Compute parameter index positions for state-space model
#'
#' Internal helper to generate index positions for parameters in the
#' optimization vector used in \code{tempssm()}. The indices are used
#' to extract autoregressive coefficients and variance parameters.
#'
#' The mapping depends on whether a seasonal component is included.
#'
#' @inheritParams tempssm
#'
#' @return A named list with components:
#' \describe{
#'   \item{ar}{Integer vector of indices for AR coefficients}
#'   \item{var}{Integer index for AR process variance}
#'   \item{H}{Integer index for observation variance}
#' }
#'
#' @details
#' The parameter vector follows different layouts depending on
#' \code{use_season}, and this function centralizes that indexing logic
#' to avoid duplication across model construction and prediction code.
#'
#' @keywords internal
#' @noRd
.get_param_index <- function(ar_order, use_season) {
  if (use_season) {
    list(
      ar = 3:(2 + ar_order),
      var = 3 + ar_order,
      H = 4 + ar_order
    )
  } else {
    list(
      ar = 2:(1 + ar_order),
      var = 2 + ar_order,
      H = 3 + ar_order
    )
  }
}

#' srr standards: model input validation
#'
#' @srrstats {G2.0} Public entry points validate input lengths for scalar,
#' vector, and time-series inputs. Scalar arguments such as `ar_order`,
#' `use_season`, `maxit`, `reltol`, cross-validation window parameters,
#' logical controls, diagnostic options, file paths, JMA area identifiers, and
#' aggregation thresholds are required to have length one. Time-series inputs
#' are validated through shared preprocessing routines, including checks that
#' exogenous series have the same number of observations and time index as the
#' response series. Arguments that may be single- or multi-valued, such as
#' `label` or `exo_name`, are checked against the corresponding number of
#' variables.
#'
#' @srrstats {G2.0a} Function documentation explicitly states length
#' expectations for scalar inputs, time-series inputs, and name or label
#' vectors. Examples include scalar documentation for `ar_order`, `use_season`,
#' cross-validation window parameters, plotting and diagnostic controls, and
#' length-one-or-number-of-variables expectations for labels and exogenous
#' variable names.
#'
#' @srrstats {G2.1} Public entry points assert expected input types, including
#' `ts` objects for temperature and exogenous time series, `zoo` objects for
#' daily data aggregation, `tempssm` objects for model summaries, diagnostics,
#' plotting, and information criteria, numeric scalars or vectors for model
#' controls and parameter inputs, logical scalars for switches, character
#' scalars or vectors for names, file paths, and JMA area identifiers, and
#' lists for cross-validation folds and results.
#'
#' @srrstats {G2.1a} Function documentation states expected input types for
#' public arguments, including class expectations for time-series objects and
#' scalar or vector type expectations for numeric, logical, and character
#' inputs.
#'
#' @srrstats {G2.2} Public functions distinguish univariate and multivariate
#' time-series inputs. Temperature response series supplied to modelling,
#' cross-validation, climatology, anomaly, and plotting functions are required
#' to be univariate `ts` objects. Exogenous series are explicitly allowed to be
#' univariate or multivariate where appropriate, and helper functions that
#' operate on multivariate `ts` objects document and check those expectations.
#'
#' @srrstats {G2.8} Conversion and standardization are performed at package
#' entry points before data are passed to analytical sub-functions. The
#' modelling and cross-validation APIs call `.tempssm_prepare_model_inputs()`,
#' which validates temperature and exogenous inputs, assigns default exogenous
#' names when explicitly allowed, checks frequency and time-index consistency,
#' applies the requested temperature missing-value policy, rejects missing
#' exogenous covariates, and returns a uniform list with validated base R `ts`
#' objects plus common `frequency` and `n_obs` metadata. Tabular, CSV, and
#' `zoo` conversion helpers similarly return regular `ts` objects before
#' those data are used by modelling or utility functions.
#'
#' @srrstats {G2.13} Missing data are checked during initial preprocessing
#' before analytical routines are called. `tempssm()` and
#' `ts_train_test_split()` both call `.tempssm_prepare_model_inputs()`, which
#' checks explicit missing response values through
#' `.tempssm_handle_missing_ts()` and rejects missing exogenous covariates
#' before model fitting, state-space construction, or rolling splits proceed.
#' Monthly tabular and CSV conversion reject missing date/index fields, insert
#' implicit missing months as explicit `NA` values with a warning, and
#' preserve explicit missing temperatures. Daily `zoo` aggregation validates
#' non-missing date indexes, summarizes missing daily values according to
#' `na.rm` and `na_prop_max`, and warns when many monthly values are missing
#' after aggregation.
#'
#' @srrstats {G2.14} User-facing modelling and cross-validation functions
#' provide explicit options for handling missing temperature observations
#' through the `na_action` argument. The shared preprocessing helper
#' `.tempssm_handle_missing_ts()` implements the temperature-series policy
#' before model fitting or rolling-origin splitting proceeds. Missing
#' exogenous covariates are always rejected because KFAS regression terms
#' require complete covariate values. Conversion utilities also expose
#' missing-value controls where relevant, such as `na.rm` and `na_prop_max` in
#' daily-to-monthly `zoo` aggregation.
#'
#' @srrstats {G2.14a} `na_action = "error"` stops when explicit missing
#' values are detected in temperature `ts` inputs. Missing values in
#' exogenous `ts` inputs always stop preprocessing, regardless of
#' `na_action`.
#'
#' @srrstats {G2.14b} The default `na_action = "inform"` reports and proceeds
#' with explicit `NA` temperature observations because linear Gaussian
#' state-space models can treat them as unobserved responses. `na_action =
#' "warn"` raises the notification to a warning, while `na_action = "allow"`
#' proceeds silently.
#'
#' @srrstats {G2.16} Undefined numeric values are handled separately from
#' explicit missing `NA` observations. Model and cross-validation inputs reject
#' `NaN`, `Inf`, and `-Inf` values during shared preprocessing via
#' `.tempssm_check_no_undefined()`. This preserves the package distinction
#' between missing response observations, which linear Gaussian state-space
#' models can handle, and undefined numeric values, which would make
#' likelihood evaluation and regression covariates ill-defined. Tabular and
#' daily `zoo` conversion helpers apply the same check before constructing
#' regular `ts` objects. Unit tests cover undefined values in temperature,
#' exogenous, data-frame, and `zoo` inputs.
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
#' @srrstats {TS2.1} The main modelling and cross-validation entry points
#' provide a \code{na_action} argument controlling how explicit missing values
#' in the regular temperature `ts` input are handled. The argument is passed
#' through the shared input pre-processing routine
#' `.tempssm_prepare_model_inputs()`, so `tempssm()` and
#' `ts_train_test_split()` use the same missing-data policy. Missing
#' exogenous covariates are not controlled by `na_action`; they are rejected
#' because regression covariates must be complete. Conversion utilities also
#' expose missing-data controls where aggregation is performed, such as
#' `na.rm` and `na_prop_max` in `daily_zoo_to_monthly_ts()`.
#'
#' @srrstats {TS2.1a} Setting `na_action = "error"` in `tempssm()` or
#' `ts_train_test_split()` stops during input pre-processing if explicit
#' missing values are detected in `temp_data`. Missing values in `exo_data`
#' always stop input pre-processing.
#'
#' @srrstats {TS2.1b} The default `na_action = "inform"` issues an informational
#' message and proceeds with explicit `NA` temperature observations. Setting
#' `na_action = "warn"` issues a warning instead, and `na_action = "allow"`
#' proceeds silently. These paths preserve the original regular `ts` object
#' and its explicit missing response values, so results use the same time index
#' rather than an implicitly shortened series.
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
#' @noRd
NULL


#' Compute expected length of the initial parameter vector
#'
#' The initial parameter vector contains trend, optional seasonal,
#' autoregressive, AR variance, and observation variance parameters.
#'
#' @inheritParams tempssm
#'
#' @return Integer scalar giving the expected number of initial values.
#'
#' @keywords internal
#' @noRd
.tempssm_inits_length <- function(ar_order, use_season) {
  as.integer(ar_order + if (use_season) 4L else 3L)
}


#' Generate default initial parameter values for \code{tempssm()}
#'
#' The returned vector follows the optimization parameter order:
#' trend variance, optional seasonal variance, AR coefficients, AR variance,
#' and observation variance.
#'
#' @inheritParams tempssm
#'
#' @return Numeric vector of default initial parameter values.
#'
#' @keywords internal
#' @noRd
.default_tempssm_inits <- function(ar_order, use_season) {
  ar_order <- as.integer(ar_order)
  ar_inits <- c(0.5, rep(0, ar_order - 1L))

  if (use_season) {
    c(
      -13, # trend
      -7, # seasonal
      ar_inits,
      -0.3, # AR variance
      -5 # observation variance
    )
  } else {
    c(
      -13, # trend
      ar_inits,
      -0.3, # AR variance
      -5 # observation variance
    )
  }
}


#' Prepare initial parameter values for \code{tempssm()}
#'
#' This helper supplies default initial values when \code{inits = NULL} and
#' validates user-supplied initial values against the model structure implied
#' by \code{ar_order} and \code{use_season}.
#'
#' @inheritParams tempssm
#'
#' @return Numeric vector of initial parameter values.
#'
#' @keywords internal
#' @noRd
.prepare_tempssm_inits <- function(inits, ar_order, use_season) {
  expected_len <- .tempssm_inits_length(ar_order, use_season)

  if (is.null(inits)) {
    .tempssm_cli_debug("Using default initial parameter values")
    inits <- .default_tempssm_inits(ar_order, use_season)
  }

  .tempssm_check_numeric(inits, "inits")

  if (length(inits) != expected_len) {
    cli::cli_abort(
      "{.arg inits} must be a numeric vector of length {expected_len}."
    )
  }

  inits
}


#' Check stationarity of autoregressive coefficients
#'
#' @param ar_coefs Numeric vector of autoregressive coefficients.
#' @param tol Numeric tolerance used for root-modulus comparisons.
#'
#' @return Logical scalar; \code{TRUE} if the AR polynomial roots are outside
#'   the unit circle within numerical tolerance.
#'
#' @keywords internal
#' @noRd
.tempssm_is_stationary_ar <- function(ar_coefs,
                                      tol = sqrt(.Machine$double.eps)) {
  if (length(ar_coefs) == 0) {
    return(TRUE)
  }

  if (any(!is.finite(ar_coefs))) {
    return(FALSE)
  }

  roots <- base::polyroot(c(1, -ar_coefs))
  all(Mod(roots) > 1 - tol)
}


#' Transform AR parameters and check stationarity
#'
#' @param ar_pars Numeric vector of unconstrained autoregressive parameters.
#'
#' @return Numeric vector of stationary autoregressive coefficients.
#'
#' @keywords internal
#' @noRd
.tempssm_transform_ar <- function(ar_pars) {
  ar_coefs <- KFAS::artransform(ar_pars)

  if (!.tempssm_is_stationary_ar(ar_coefs)) {
    cli::cli_abort(
      "Transformed autoregressive coefficients are not stationary."
    )
  }

  ar_coefs
}


#' Transform unconstrained parameters to constrained values
#'
#' Internal helper to apply standard transformations to the parameter vector
#' used in state-space model optimization. This centralizes the logic for
#' exponentiating variance parameters and transforming autoregressive
#' coefficients to ensure stationarity.
#'
#' @param pars Numeric vector of unconstrained parameters
#' @param ar_idx Integer vector of indices for AR coefficients
#' @param var_idx Integer index for AR process variance parameter
#' @param H_idx Integer index for observation variance parameter
#' @inheritParams tempssm
#'
#' @return A named list containing transformed parameters:
#' \describe{
#'   \item{trend_var}{Positive trend variance: \code{exp(pars[[1]])}}
#'   \item{season_var}{Positive seasonal variance if \code{use_season = TRUE},
#'   otherwise \code{NULL}}
#'   \item{ar_coefs}{Stationary AR coefficients via
#'   \code{KFAS::artransform()}}
#'   \item{ar_var}{Positive AR process variance}
#'   \item{H}{Positive observation variance}
#' }
#'
#' @details
#' All variance parameters are exponentiated to ensure positivity.
#' AR coefficients are transformed using \code{KFAS::artransform()} to
#' ensure the autoregressive process satisfies stationarity constraints.
#'
#' This function centralizes the transformation logic to avoid duplication
#' across \code{.define_update_func()} and \code{.build_newdata_ssm()}.
#'
#' @keywords internal
#' @noRd
.transform_parameters <- function(pars, ar_idx, var_idx, H_idx, use_season) {
  list(
    trend_var = exp(pars[1]),
    season_var = if (use_season) exp(pars[2]) else NULL,
    ar_coefs = .tempssm_transform_ar(pars[ar_idx]),
    ar_var = exp(pars[var_idx]),
    H = exp(pars[H_idx])
  )
}


#' Base function for fitting a linear Gaussian state-space model to temperature
#' time series
#'
#' This function estimates a linear Gaussian state-space model (SSM)
#' for regular temperature time series using Kalman filtering and smoothing.
#' The input must be a univariate \code{ts} object with an integer seasonal
#' frequency greater than 1.
#'
#' @details
#' The observed temperature series is not required to be stationary in mean.
#' Long-term mean changes are represented by latent level and drift
#' components, and recurring within-cycle mean variation is represented by the
#' seasonal component when \code{use_season = TRUE}. Stationarity constraints
#' are applied to the autoregressive component, whose coefficients are
#' transformed with \code{KFAS::artransform()} during model fitting.
#'
#' Observation and state disturbances are assumed Gaussian with positive,
#' time-invariant variances. The model does not currently include
#' time-varying volatility or automatic stationarity tests for the observed
#' input series.
#'
#' @param temp_data A temperature time series of class \code{ts}.
#'   The \code{ts} object must be univariate.
#'   The series can have any integer frequency of 2 or higher.
#'   For example, a frequency of 12 represents a monthly \code{ts} object.
#'
#' @param exo_data A data set of exogenous variable(s) of class \code{ts}.
#'   The series may have any integer frequency of 2 or higher,
#'   but it must be the same as that of \code{temp_data}. It must also have
#'   the same number of observations and the same time index as
#'   \code{temp_data}. The default is \code{NULL} when fitting a model without
#'   exogenous variables.
#'
#' @param ar_order Integer scalar specifying the order of the autoregressive
#' (AR)
#' component in the error structure (e.g., 2 for AR(2), 3 for AR(3)).
#' Defaults to 1. Orders from 1 to 4 are intended as the usual working range;
#' larger values are allowed but trigger a warning because they may lead to
#' unstable estimation.
#'
#' @param use_season Logical scalar; should the seasonal component be
#' considered? Defaults to \code{TRUE}.
#'
#' @param inits Optional numeric vector of initial parameter values.
#'  If \code{NULL}, default values are used. When supplied, its length must
#'  match the number of variance and autoregressive parameters implied by
#'  \code{ar_order} and \code{use_season}.
#'
#' @param maxit Optional numeric scalar giving the maximum number of
#' iterations.
#'  If \code{NULL}, default value of 5000 is used.
#'
#' @param reltol Optional numeric scalar giving the relative convergence
#' tolerance.
#'  If \code{NULL}, default value of 1e-16 is used.
#'
#' @param na_action Character scalar specifying how explicit missing
#'   observations in \code{temp_data} should be handled. Use \code{"inform"}
#'   to issue an informational message and proceed, \code{"warn"} to issue a
#'   warning and proceed, \code{"error"} to stop, or \code{"allow"} to proceed
#'   silently. The default is \code{"inform"}. Missing values in
#'   \code{exo_data} are always rejected; exogenous covariates must be completed
#'   before model fitting.
#'
#' @return An object of class \code{"tempssm"}, a named list containing:
#' \describe{
#'   \item{model}{Fitted \code{SSModel} object.}
#'   \item{fit}{Results from \code{fitSSM}.}
#'   \item{kfs}{Kalman filtering and smoothing results from \code{KFS}.}
#'   \item{data_temp}{Temperature time series used for estimation.}
#'   \item{data_exo}{Time series of exogenous variable(s) used for estimation.}
#'   \item{ar_order}{Order of the autoregressive component.}
#'   \item{use_season}{Logical; whether to include a seasonal component.}
#'   \item{call}{Matched function call.}
#' }
#'
#' @import KFAS
#' @export
#'
#' @examples
#' \dontrun{
#' data(niigata_sst)
#' res <- tempssm(niigata_sst)
#' summary(res)
#' }
tempssm <- function(temp_data,
                    exo_data = NULL,
                    ar_order = 1,
                    inits = NULL,
                    maxit = NULL,
                    reltol = NULL,
                    use_season = TRUE,
                    na_action = c("inform", "warn", "error", "allow")) {
  exo_name <- NULL
  state_names <- character(0)

  .tempssm_check_length_one(ar_order, "ar_order")
  .tempssm_check_length_one(use_season, "use_season")
  .tempssm_check_length_one(maxit, "maxit", allow_null = TRUE)
  .tempssm_check_length_one(reltol, "reltol", allow_null = TRUE)
  .tempssm_check_numeric(ar_order, "ar_order")
  .tempssm_check_logical(use_season, "use_season")
  .tempssm_check_numeric(maxit, "maxit", allow_null = TRUE)
  .tempssm_check_numeric(reltol, "reltol", allow_null = TRUE)
  na_action <- match.arg(na_action)

  if (!(is.numeric(ar_order) &&
    !is.na(ar_order) &&
    ar_order >= 1 && .tempssm_is_integerish(ar_order))) {
    cli::cli_abort(
      "The argument {.arg ar_order} must be an integer >= 1."
    )
  }

  if (!is.logical(use_season) || is.na(use_season)) {
    cli::cli_abort(
      "{.arg use_season} must be a logical scalar."
    )
  }

  ## ---- Default initial values --------------------------------------
  inits <- .prepare_tempssm_inits(
    inits = inits,
    ar_order = ar_order,
    use_season = use_season
  )

  ## ---- Default values for optimization parameters --------------------
  if (is.null(maxit)) {
    maxit <- 5000
  }
  if (is.null(reltol)) {
    reltol <- 1e-16
  }

  if (!is.null(maxit) &&
    (!is.numeric(maxit) || is.na(maxit) || !is.finite(maxit))) {
    cli::cli_abort(
      "{.arg maxit} must be a finite numeric scalar."
    )
  }

  if (!is.null(reltol) &&
    (!is.numeric(reltol) || is.na(reltol) || !is.finite(reltol))) {
    cli::cli_abort(
      "{.arg reltol} must be a finite numeric scalar."
    )
  }

  tryCatch(
    {
      ## ---- Start message -----------------------------------------------
      .tempssm_cli_inform(
        paste0(
          "Fitting state-space model (",
          "AR order: {ar_order}, ",
          "seasonal: {use_season})"
        )
      )  

      ## ---- Input checks -------------------------------------------------
      model_inputs <- .tempssm_prepare_model_inputs(
        temp_data = temp_data,
        exo_data = exo_data,
        na_action = na_action
      )

      temp_data <- model_inputs$temp_data
      exo_data <- model_inputs$exo_data
      y <- temp_data
      freq <- model_inputs$frequency

      .tempssm_cli_debug("Time series frequency: {freq}")

      .tempssm_check_length_one(ar_order, "ar_order")
      if (!(is.numeric(ar_order) &&
        !is.na(ar_order) &&
        ar_order >= 1 && .tempssm_is_integerish(ar_order))) {
        cli::cli_abort(
          "The argument {.arg ar_order} must be an integer >= 1."
        )
      }

      .tempssm_check_length_one(use_season, "use_season")
      if (!is.logical(use_season) || is.na(use_season)) {
        cli::cli_abort(
          "{.arg use_season} must be a logical scalar."
        )
      }

      if (ar_order > 4) {
        cli::cli_warn(
          "An {.arg ar_order} greater than 4 may lead to unstable estimation."
        )
      }

      ## ---- Defaults -----------------------------------------------------
      if (is.null(maxit)) {
        maxit <- 5000
      } else {
        .tempssm_check_length_one(maxit, "maxit")
        if (!is.numeric(maxit) || is.na(maxit) || !is.finite(maxit)) {
          cli::cli_abort(
            "{.arg maxit} must be a finite numeric scalar."
          )
        }
      }

      if (is.null(reltol)) {
        reltol <- 1e-16
      } else {
        .tempssm_check_length_one(reltol, "reltol")
        if (!is.numeric(reltol) || is.na(reltol) || !is.finite(reltol)) {
          cli::cli_abort(
            "{.arg reltol} must be a finite numeric scalar."
          )
        }
      }

      .tempssm_cli_debug(
        "Optimization settings: maxit={maxit}, reltol={reltol}"
      )

      ## ---- Parameter indexing ------------------------------------------
      param_idx_list <- .get_param_index(ar_order = ar_order,
                       use_season = use_season)

      ar_idx <- param_idx_list$ar
      var_idx <- param_idx_list$var
      H_idx <- param_idx_list$H

      ## ---- Exogenous handling ------------------------------------------
      if (is.null(exo_data)) {
        .tempssm_cli_debug("No exogenous variables used")

        exo_name <- NULL
        exo_mat <- NULL
      } else {
        .tempssm_cli_inform("Including exogenous variables in the model")

        exo_name <- colnames(exo_data)
        exo_mat <- as.matrix(exo_data)

        .tempssm_cli_debug(
          "Exogenous variables: {paste(exo_name, collapse = ', ')}"
        )
      }

      ## ---- Model definition --------------------------------------------
      build_ssm <- .define_build_model(
        y = y,
        freq = freq,
        use_season = use_season,
        exo_mat = exo_mat,
        ar_order = ar_order
      )

      update_func_common <- .define_update_func(
        y = y,
        freq = freq,
        use_season = use_season,
        exo_mat = exo_mat,
        ar_order = ar_order,
        ar_idx = ar_idx,
        var_idx = var_idx,
        H_idx = H_idx
      )

      ## ---- Optimization -------------------------------------------------
      .tempssm_cli_inform("Optimizing model parameters (stage 1)")

      fit1 <- fitSSM(
        build_ssm,
        inits = inits,
        updatefn = update_func_common,
        method = "Nelder-Mead",
        control = list(maxit = maxit, reltol = reltol)
      )

      .tempssm_cli_inform("Refining optimization (stage 2)")

      fit2 <- fitSSM(
        build_ssm,
        inits = fit1$optim.out$par,
        updatefn = update_func_common,
        method = "BFGS",
        control = list(maxit = maxit, reltol = reltol)
      )

      ## ---- KFS ----------------------------------------------------------
      .tempssm_cli_inform("Running Kalman filtering and smoothing")

      kfs <- KFS(
        fit2$model,
        filtering = c("state", "mean"),
        smoothing = c("state", "mean", "disturbance")
      )

      ## ---- State names --------------------------------------------------
      state_names <- character(ncol(kfs$alphahat))
      idx <- 1

      if (!is.null(exo_name)) {
        state_names[idx:(idx + length(exo_name) - 1)] <- exo_name
        idx <- idx + length(exo_name)
      }

      state_names[idx] <- "level"
      idx <- idx + 1
      state_names[idx] <- "slope"
      idx <- idx + 1

      if (use_season) {
        n_season <- freq - 1
        state_names[idx:(idx + n_season - 1)] <-
          paste0("sea_dummy", seq_len(n_season))
        idx <- idx + n_season
      }

      if (ar_order > 0) {
        state_names[idx:(idx + ar_order - 1)] <-
          paste0("arima", seq_len(ar_order))
      }

      ## ---- Output -------------------------------------------------------
      out <- list(
        model = fit2$model,
        fit = fit2,
        kfs = kfs,
        temp_data = temp_data,
        exogenous_data = exo_data,
        ar_order = ar_order,
        use_season = use_season,
        call = match.call(),
        converged = fit2$optim.out$convergence == 0,
        state_map = list(
          exogenous = exo_name,
          all       = state_names
        )
      )

      colnames(out$kfs$alphahat) <- state_names
      class(out) <- "tempssm"

      ## ---- Completion message ------------------------------------------
      if (out$converged) {
        .tempssm_cli_inform("Model fitting completed successfully")
      } else {
        cli::cli_warn("Model fitting finished but did not converge")
      }

      return(out)
    },
    error = function(e) {
      cli::cli_warn(
        "Model fitting failed: {conditionMessage(e)}"
      )

      out <- list(
        model = NULL,
        fit = NULL,
        kfs = NULL,
        temp_data = temp_data,
        exogenous_data = exo_data,
        ar_order = ar_order,
        use_season = use_season,
        call = match.call(),
        converged = FALSE,
        state_map = list(
          exogenous = exo_name,
          all       = state_names
        )
      )

      class(out) <- "tempssm"
      return(out)
    }
  ) # close tryCatch
}


#' Internal helper to construct a KFAS state-space model
#'
#' This internal function builds and returns a \code{KFAS::SSModel} object
#' corresponding to the specified model structure. It conditionally includes
#' trend, seasonal, autoregressive, and optional exogenous components,
#' depending on the input arguments.
#'
#' The function is primarily used within \code{tempssm()} to generate the
#' baseline state-space model prior to parameter estimation. Model parameters
#' are left unspecified (set to \code{NA} or zero) and are later updated
#' during optimization via the update function.
#'
#' Four model configurations are supported:
#' \itemize{
#'   \item Trend + Seasonal + AR (no exogenous variables)
#'   \item Trend + AR (no seasonal component)
#'   \item Trend + Seasonal + AR + Exogenous variables
#'   \item Trend + AR + Exogenous variables (no seasonal component)
#' }
#'
#' @param y A numeric vector or univariate \code{ts} object representing
#'   the observed time series.
#' @param freq Integer indicating the seasonal frequency of \code{y}.
#' @inheritParams tempssm
#' @param exo_mat Optional numeric matrix of exogenous regressors.
#'   If provided, each column is treated as a separate covariate.
#'
#' @return A \code{KFAS::SSModel} object with unspecified variance parameters.
#'
#' @details
#' The trend component is modeled as a second-order polynomial trend.
#' The seasonal component (if included) uses a dummy-variable formulation
#' with a sum-to-zero constraint. The autoregressive component is implemented
#' via \code{SSMarima()} with coefficients initialized to zero.
#'
#' @keywords internal
#' @noRd
.define_build_model <- function(y = NULL,
                                freq = freq,
                                use_season,
                                exo_mat,
                                ar_order = 1) {
  if (use_season && is.null(exo_mat)) {
    build_ssm <- KFAS::SSModel(
      y ~
        SSMtrend(
          degree = 2,
          Q = c(list(0), list(NA))
        ) +
        SSMseasonal(
          sea.type = "dummy",
          period = freq,
          Q = NA
        ) +
        SSMarima(
          ar = rep(0, ar_order),
          d = 0,
          Q = NA
        ),
      H = NA
    )
  } else if (!use_season && is.null(exo_mat)) {
    build_ssm <- KFAS::SSModel(
      y ~
        SSMtrend(
          degree = 2,
          Q = c(list(0), list(NA))
        ) +
        SSMarima(
          ar = rep(0, ar_order),
          d = 0,
          Q = NA
        ),
      H = NA
    )
  } else if (use_season && !(is.null(exo_mat))) {
    build_ssm <- KFAS::SSModel(
      H = NA,
      y ~ exo_mat +
        SSMtrend(
          degree = 2,
          Q = c(list(0), list(NA))
        ) +
        SSMseasonal(
          sea.type = "dummy",
          period = freq,
          Q = NA
        ) +
        SSMarima(
          ar = rep(0, ar_order),
          d = 0,
          Q = NA
        )
    )
  } else if (!use_season && !(is.null(exo_mat))) {
    build_ssm <- KFAS::SSModel(
      H = NA,
      y ~ exo_mat +
        SSMtrend(
          degree = 2,
          Q = c(list(0), list(NA))
        ) +
        SSMarima(
          ar = rep(0, ar_order),
          d = 0,
          Q = NA
        )
    )
  }
  return(build_ssm)
}


#' Internal helper to generate the parameter update function for KFAS
#'
#' This internal function constructs and returns an update function to be used
#' in \code{KFAS::fitSSM()}, which maps an unconstrained parameter vector
#' to a valid state-space model.
#'
#' The returned function takes model parameters on an unconstrained scale
#' and transforms them to enforce constraints such as:
#' \itemize{
#'   \item Variance parameters are exponentiated to ensure positivity
#'   \item Autoregressive coefficients are transformed using
#'         \code{KFAS::artransform()} to ensure stationarity
#' }
#'
#' The structure of the update function depends on whether seasonal and/or
#' exogenous components are included, ensuring consistency with the model
#' defined in \code{.define_build_model()}.
#'
#' @param y A numeric vector or univariate \code{ts} object representing
#'   the observed time series.
#' @param freq Integer indicating the seasonal frequency.
#' @inheritParams tempssm
#' @param exo_mat Optional matrix of exogenous regressors.
#' @param ar_idx Integer vector indicating positions of AR coefficients
#'   in the parameter vector.
#' @param var_idx Integer indicating the position of the AR process variance.
#' @param H_idx Integer indicating the position of the observation variance.
#'
#' @return A function with signature \code{function(pars, model)} suitable
#'   for use in \code{KFAS::fitSSM()}.
#'
#' @details
#' The returned function rebuilds the full state-space model using the
#' transformed parameters. It ensures that all variance parameters are positive
#' and that the autoregressive process satisfies stationarity constraints.
#'
#' This closure-based design allows efficient reuse of model structure during
#' numerical optimization without reconstructing external inputs.
#'
#' @keywords internal
#' @noRd
.define_update_func <- function(y = NULL,
                                freq = freq,
                                use_season,
                                exo_mat,
                                ar_order = 1,
                                ar_idx,
                                var_idx,
                                H_idx) {
  if (use_season && is.null(exo_mat)) {
    update_func <- function(pars, model) {
      trans <- .transform_parameters(pars, ar_idx, var_idx, H_idx, use_season)
      return(
        KFAS::SSModel(
          y ~
            SSMtrend(
              degree = 2,
              Q = c(list(0), list(trans$trend_var))
            ) +
            SSMseasonal(
              sea.type = "dummy",
              period = freq,
              Q = trans$season_var
            ) +
            SSMarima(
              ar = trans$ar_coefs,
              d = 0,
              Q = trans$ar_var
            ),
          H = trans$H
        )
      )
    }
  } else if (!use_season && is.null(exo_mat)) {
    update_func <- function(pars, model) {
      trans <- .transform_parameters(pars, ar_idx, var_idx, H_idx, use_season)
      return(
        KFAS::SSModel(
          y ~
            SSMtrend(
              degree = 2,
              Q = c(list(0), list(trans$trend_var))
            ) +
            SSMarima(
              ar = trans$ar_coefs,
              d = 0,
              Q = trans$ar_var
            ),
          H = trans$H
        )
      )
    }
  } else if (use_season && !(is.null(exo_mat))) {
    update_func <- function(pars, model) {
      trans <- .transform_parameters(pars, ar_idx, var_idx, H_idx, use_season)
      return(
        KFAS::SSModel(
          H = trans$H,
          y ~ exo_mat +
            SSMtrend(
              degree = 2,
              Q = c(list(0), list(trans$trend_var))
            ) +
            SSMseasonal(
              sea.type = "dummy",
              period = freq,
              Q = trans$season_var
            ) +
            SSMarima(
              ar = trans$ar_coefs,
              d = 0,
              Q = trans$ar_var
            )
        )
      )
    }
  } else if (!use_season && !(is.null(exo_mat))) {
    update_func <- function(pars, model) {
      trans <- .transform_parameters(pars, ar_idx, var_idx, H_idx, use_season)
      return(
        KFAS::SSModel(
          H = trans$H,
          y ~ exo_mat +
            SSMtrend(
              degree = 2,
              Q = c(list(0), list(trans$trend_var))
            ) +
            SSMarima(
              ar = trans$ar_coefs,
              d = 0,
              Q = trans$ar_var
            )
        )
      )
    }
  }
}


#' Check ts object of temperature time series for applying \code{tempssm()}
#'
#' @inheritParams tempssm
#'
#' @return A univariate \code{ts} object.
#'
#' @keywords internal
#' @noRd
.tempssm_check_temp_ts <- function(temp_data) {
  ## ---- type check -----------------------------------------------------
  if (!inherits(temp_data, "ts")) {
    cli::cli_abort(
      "The object {.arg temp_data} must be a {.cls ts} object."
    )
  }

  temp_data <- .strip_units_ts(temp_data, arg_name = "temp_data")
  if (!is.numeric(temp_data)) {
    cli::cli_abort(
      "The object {.arg temp_data} must contain numeric values."
    )
  }
  .tempssm_check_no_undefined(temp_data, "temp_data")

  ## ---- frequency check ------------------------------------------------
  freq <- frequency(temp_data)
  freq_int <- as.integer(round(freq))

  if (freq <= 1 || abs(freq - freq_int) > sqrt(.Machine$double.eps)) {
    cli::cli_abort(
      "The procedure requires a {.cls ts} object with integer frequency > 1."
    )
  }

  ## ---- time order check ----------------------------------------------
  .tempssm_check_ts_order(temp_data, arg_name = "temp_data")

  ## ---- univariate check ----------------------------------------------
  .tempssm_check_univariate_ts(temp_data, "temp_data")

  ## ---- debug message --------------------------------------------------
  .tempssm_cli_debug(
    "Validated temp_data: univariate ts with frequency {freq}"
  )

  return(temp_data)
}


#' Check strict ordering of a \code{ts} time index
#'
#' @param x A \code{ts} object.
#' @param arg_name Name of the argument being checked.
#'
#' @return Invisibly returns \code{x}.
#'
#' @keywords internal
#' @noRd
.tempssm_check_ts_order <- function(x, arg_name) {
  time_index <- time(x)

  if (length(time_index) > 1 && any(diff(time_index) <= 0)) {
    cli::cli_abort(
      "Time index of {.arg {arg_name}} must be strictly increasing."
    )
  }

  invisible(x)
}


#' Handle explicit missing values in a \code{ts} object
#'
#' @param x A \code{ts} object.
#' @param arg_name Name of the argument being checked.
#' @param na_action Missing-value handling policy.
#'
#' @return Invisibly returns \code{x}.
#'
#' @keywords internal
#' @noRd
.tempssm_handle_missing_ts <- function(x, arg_name, na_action) {
  if (!anyNA(x)) {
    return(invisible(x))
  }

  if (all(is.na(x))) {
    cli::cli_abort(
      "{.arg {arg_name}} must contain at least one non-missing value."
    )
  }

  if (identical(na_action, "error")) {
    cli::cli_abort(
      "Missing values detected in {.arg {arg_name}}."
    )
  }

  msg <- c(
    "Missing values detected in {.arg {arg_name}}.",
    "i" = paste0(
      "They will be retained and treated as unobserved responses during ",
      "Kalman filtering and smoothing."
    )
  )

  if (identical(na_action, "inform")) {
    .tempssm_cli_inform(msg)
  }

  if (identical(na_action, "warn")) {
    cli::cli_warn(msg)
  }

  invisible(x)
}


#' Check ts object of exogenous variable(s) for applying \code{tempssm()}
#'
#' @inheritParams tempssm
#'
#' @return A univariate or multivariate \code{ts} object.
#'
#' @keywords internal
#' @noRd
.tempssm_check_exo_ts <- function(temp_data,
                                  exo_data) {
  ## ---- check temp_data -----------------------------------------------
  temp_data_checked <- .tempssm_check_temp_ts(temp_data)
  temp_freq <- frequency(temp_data_checked)

  ## ---- type check -----------------------------------------------------
  if (!inherits(exo_data, "ts")) {
    cli::cli_abort(
      "The object {.arg exo_data} must be a {.cls ts} object."
    )
  }

  exo_data <- .strip_units_ts(exo_data, arg_name = "exo_data")
  if (!is.numeric(exo_data)) {
    cli::cli_abort(
      "The object {.arg exo_data} must contain numeric values."
    )
  }
  .tempssm_check_no_undefined(exo_data, "exo_data")

  ## ---- frequency check ------------------------------------------------
  exo_freq <- frequency(exo_data)
  if (!isTRUE(all.equal(exo_freq, temp_freq))) {
    cli::cli_abort(
      "Frequency of {.arg exo_data} must match that of {.arg temp_data}."
    )
  }

  ## ---- time order check ----------------------------------------------
  .tempssm_check_ts_order(exo_data, arg_name = "exo_data")

  ## ---- length check ---------------------------------------------------
  if (NROW(exo_data) != NROW(temp_data_checked)) {
    cli::cli_abort(
      "Length of {.arg exo_data} must match that of {.arg temp_data}."
    )
  }

  ## ---- time index check ----------------------------------------------
  if (!isTRUE(all.equal(time(exo_data), time(temp_data_checked)))) {
    cli::cli_abort(
      "Time index of {.arg exo_data} must match that of {.arg temp_data}."
    )
  }

  ## ---- dimensionality check ------------------------------------------
  n_col <- NCOL(exo_data)
  n_obs <- NROW(temp_data_checked)
  if (n_col > n_obs) {
    cli::cli_abort(
      paste0(
        "The number of exogenous variables must not exceed ",
        "the number of observations."
      )
    )
  }

  ## ---- column names check --------------------------------------------
  if (is.null(colnames(exo_data))) {
    cli::cli_abort(
      "The object {.arg exo_data} must have column name(s)."
    )
  }

  ## ---- debug message --------------------------------------------------
  uni_multi <- if (n_col > 1) "multivariate" else "univariate"

  .tempssm_cli_debug(
    paste0(
      "Validated exo_data: {uni_multi} ts with {n_col} variable{?s}, ",
      "frequency {exo_freq}"
      )
   )

  return(exo_data)
}


#' Validate and standardize model time-series inputs
#'
#' @inheritParams tempssm
#' @param allow_unnamed_exo Logical; if \code{TRUE}, unnamed exogenous
#'   variables are allowed after assigning default names.
#' @param default_exo_names Logical; if \code{TRUE}, default names are assigned
#'   to unnamed exogenous variables.
#' @param na_action Missing-value handling policy.
#'
#' @return A named list containing validated \code{temp_data}, validated or
#'   \code{NULL} \code{exo_data}, the common frequency, and the number of
#'   observations.
#'
#' @keywords internal
#' @noRd
.tempssm_prepare_model_inputs <- function(temp_data,
                                          exo_data = NULL,
                                          allow_unnamed_exo = FALSE,
                                          default_exo_names = FALSE,
                                          na_action = c("inform",
                                                        "warn",
                                                        "error",
                                                        "allow")) {
  na_action <- match.arg(na_action)

  temp_data <- .tempssm_check_temp_ts(temp_data)
  .tempssm_handle_missing_ts(
    temp_data,
    arg_name = "temp_data",
    na_action = na_action
  )

  if (!is.null(exo_data) && is.null(colnames(exo_data))) {
    if (!allow_unnamed_exo || !default_exo_names) {
      cli::cli_abort(
        "The object {.arg exo_data} must have column name(s)."
      )
    }

    cli::cli_warn(
      "No column names in {.arg exo_data}; assigning default names."
    )
    exo_data <- tempssm::set_ts_name(
      exo_data,
      label = paste0("var", seq_len(NCOL(exo_data))),
      quiet = TRUE
    )
  }

  if (!is.null(exo_data)) {
    exo_data <- .tempssm_check_exo_ts(
      temp_data = temp_data,
      exo_data = exo_data
    )
    if (anyNA(exo_data)) {
      cli::cli_abort(
        c(
          "Missing values detected in {.arg exo_data}.",
          "i" = paste0(
            "Exogenous covariates must be complete for KFAS regression ",
            "terms."
          ),
          "*" = paste0(
            "Please impute, remove, or otherwise complete missing ",
            "exogenous values before calling {.fn tempssm}."
          )
        )
      )
    }
  }

  list(
    temp_data = temp_data,
    exo_data = exo_data,
    frequency = frequency(temp_data),
    n_obs = length(temp_data)
  )
}
