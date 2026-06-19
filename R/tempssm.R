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
#'   but it must be the same as that of \code{temp_data}. The default is
#'   \code{NULL} when fitting a model without exogenous variables.
#'
#' @param ar_order Integer specifying the order of the autoregressive (AR)
#' component in the error structure (e.g., 2 for AR(2), 3 for AR(3)).
#' Defaults to 1.
#'
#' @param use_season Logical; Should seasonal component be considered.
#' Default to TRUE.
#'
#' @param inits Optional numeric vector of initial parameter values.
#'  If \code{NULL}, default values are used.
#'
#' @param maxit Optional numeric of maximum iteration.
#'  If \code{NULL}, default value of 5000 is used.
#'
#' @param reltol Optional numeric of reltol.
#'  If \code{NULL}, default value of 1e-16 is used.
#'
#' @param na_action Character string specifying how explicit missing
#'   observations in \code{temp_data} or \code{exo_data} should be handled.
#'   Use \code{"warn"} to issue a warning and proceed, \code{"error"} to stop,
#'   or \code{"allow"} to proceed silently. The default is \code{"warn"}.
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
                    na_action = c("warn", "error", "allow")) {
  exo_name <- NULL
  state_names <- character(0)

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

      if (!(is.numeric(ar_order) &&
        length(ar_order) == 1 &&
        !is.na(ar_order) &&
        ar_order >= 1 && ar_order == floor(ar_order))) {
        cli::cli_abort(
          "The argument {.arg ar_order} must be an integer >= 1."
        )
      }

      if (ar_order > 4) {
        cli::cli_warn(
          "An {.arg ar_order} greater than 4 may lead to unstable estimation."
        )
      }

      ## ---- Default initial values --------------------------------------
      if (is.null(inits)) {
        .tempssm_cli_debug("Using default initial parameter values")

        ar_rep_length_minus_one <- ar_order - 1
        ar_inits <- c(0.5, rep(0, ar_rep_length_minus_one))

        if (use_season) {
          inits <- c(
            -13, # trend
            -7, # seasonal
            ar_inits,
            -0.3, # AR variance
            -5 # observation variance
          )
        } else {
          inits <- c(
            -13, # trend
            ar_inits,
            -0.3, # AR variance
            -5 # observation variance
          )
        }
      }

      expected_len <- if (use_season) {
        2 + ar_order + 2
      } else {
        1 + ar_order + 2
      }

      if (!is.numeric(inits) || length(inits) != expected_len) {
        cli::cli_abort(
          "{.arg inits} must be a numeric vector of length {expected_len}."
        )
      }

      ## ---- Defaults -----------------------------------------------------
      if (is.null(maxit)) {
        maxit <- 5000
      }

      if (is.null(reltol)) {
        reltol <- 1e-16
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
  if (!is.null(dim(temp_data)) && NCOL(temp_data) != 1) {
    cli::cli_abort(
      "The object {.arg temp_data} must be univariate."
    )
  }

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

  if (identical(na_action, "error")) {
    cli::cli_abort(
      "Missing values detected in {.arg {arg_name}}."
    )
  }

  if (identical(na_action, "warn")) {
    cli::cli_warn(
      c(
        "Missing values detected in {.arg {arg_name}}.",
        "i" = "Proceeding with explicit {.val NA} observations."
      )
    )
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

  ## ---- frequency check ------------------------------------------------
  exo_freq <- frequency(exo_data)
  if (exo_freq != temp_freq) {
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
  if (!all(time(exo_data) == time(temp_data_checked))) {
    cli::cli_abort(
      "Time index of {.arg exo_data} must match that of {.arg temp_data}."
    )
  }

  ## ---- column names check --------------------------------------------
  if (is.null(colnames(exo_data))) {
    cli::cli_abort(
      "The object {.arg exo_data} must have column name(s)."
    )
  }

  ## ---- debug message --------------------------------------------------
  n_col <- NCOL(exo_data)
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
                                          na_action = c("warn",
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
    .tempssm_handle_missing_ts(
      exo_data,
      arg_name = "exo_data",
      na_action = na_action
    )
  }

  list(
    temp_data = temp_data,
    exo_data = exo_data,
    frequency = frequency(temp_data),
    n_obs = length(temp_data)
  )
}
