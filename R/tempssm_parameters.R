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

  if (!all(is.finite(ar_coefs))) {
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


#' Extract transformed parameters from a fitted tempssm object
#'
#' @param res A validated \code{tempssm} object with fitted parameters.
#'
#' @return A named list containing variance parameters and AR coefficients.
#'
#' @keywords internal
#' @noRd
.extract_tempssm_params <- function(res) {
  pars <- res$fit$optim.out$par
  param_idx <- .get_param_index(
    ar_order = res$ar_order,
    use_season = res$use_season
  )
  season_variance <- if (res$use_season) exp(pars[2]) else NA

  list(
    H = exp(pars[param_idx$H]),
    Q_trend = exp(pars[1]),
    Q_season = season_variance,
    Q_ar = exp(pars[param_idx$var]),
    ARs = .tempssm_transform_ar(pars[param_idx$ar])
  )
}
