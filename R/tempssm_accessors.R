#' Check common inputs for tempssm accessor functions
#'
#' @param res An object expected to inherit from \code{"tempssm"}.
#' @param fun Character scalar naming the calling function.
#'
#' @return Invisibly returns \code{NULL}.
#' @noRd
.tempssm_check_accessor_input <- function(res, fun) {
  if (!inherits(res, "tempssm")) {
    stop(
      "`res` must be an object of class 'tempssm' for ",
      fun,
      "().",
      call. = FALSE
    )
  }

  invisible(NULL)
}


#' Check confidence interval arguments for tempssm accessors
#'
#' @param ci Logical; if \code{TRUE}, confidence intervals are requested.
#' @param ci_level Numeric confidence level between 0 and 1.
#' @param fun Character scalar naming the calling function.
#'
#' @return Invisibly returns \code{NULL}.
#' @noRd
.tempssm_check_accessor_ci <- function(ci, ci_level, fun) {
  if (ci) {
    if (!is.numeric(ci_level) || length(ci_level) != 1 ||
        ci_level <= 0 || ci_level >= 1) {
      stop(
        "`ci_level` must be a numeric value between 0 and 1: ",
        fun,
        "().",
        call. = FALSE
      )
    }
  }

  invisible(NULL)
}


#' Check whether a smoothed state is available
#'
#' @param res An object of class \code{"tempssm"}.
#' @param state Character scalar naming a state in \code{res$kfs$alphahat}.
#' @param message Character scalar used when the state is unavailable.
#'
#' @return Invisibly returns \code{NULL}.
#' @noRd
.tempssm_check_state_available <- function(res, state, message) {
  alphahat <- res$kfs$alphahat
  if (is.null(alphahat) || !state %in% colnames(alphahat)) {
    stop(message, call. = FALSE)
  }

  invisible(NULL)
}


#' Convert a smoothed state to a ts object
#'
#' @param res An object of class \code{"tempssm"}.
#' @param state Character scalar naming a state in \code{res$kfs$alphahat}.
#' @param scale Numeric multiplier applied to the extracted state.
#'
#' @return A univariate \code{ts} object.
#' @noRd
.tempssm_state_ts <- function(res, state, scale = 1) {
  ts(
    res$kfs$alphahat[, state] * scale,
    start = start(res$temp_data),
    frequency = frequency(res$temp_data)
  )
}


#' Convert state confidence intervals to a multivariate ts object
#'
#' @param res An object of class \code{"tempssm"}.
#' @param state Character scalar naming a state in \code{stats::confint()}.
#' @param output_name Character scalar used for the point estimate column.
#' @param ci_level Numeric confidence level between 0 and 1.
#' @param scale Numeric multiplier applied to estimates and intervals.
#' @param message Character scalar used when the state interval is unavailable.
#'
#' @return A multivariate \code{ts} object.
#' @noRd
.tempssm_state_ci_ts <- function(res, state, output_name, ci_level,
                                 scale = 1, message) {
  ci_obj <- stats::confint(res$kfs, level = ci_level)

  if (!state %in% names(ci_obj)) {
    stop(message, call. = FALSE)
  }

  point <- .tempssm_state_ts(res, state, scale = scale)
  values <- cbind(
    point,
    lwr = ci_obj[[state]][, "lwr"] * scale,
    upr = ci_obj[[state]][, "upr"] * scale
  )
  colnames(values)[1] <- output_name

  ts(
    values,
    start = start(res$temp_data),
    frequency = frequency(res$temp_data)
  )
}


#' Extract a smoothed state as a time series
#'
#' @param res An object of class \code{"tempssm"}.
#' @param state Character scalar naming the smoothed state.
#' @param output_name Character scalar used for the output column.
#' @param ci Logical; if \code{TRUE}, confidence intervals are returned.
#' @param ci_level Numeric confidence level between 0 and 1.
#' @param fun Character scalar naming the calling function.
#' @param missing_msg Character scalar used when the state is unavailable.
#' @param ci_missing_msg Character scalar used when intervals are unavailable.
#' @param scale_by_frequency Logical; if \code{TRUE}, scale by time frequency.
#'
#' @return A univariate or multivariate \code{ts} object.
#' @noRd
.tempssm_extract_state_ts <- function(res, state, output_name, ci = FALSE,
                                      ci_level = 0.95, fun, missing_msg,
                                      ci_missing_msg,
                                      scale_by_frequency = FALSE) {
  .tempssm_check_accessor_input(res, fun)
  .tempssm_check_accessor_ci(ci, ci_level, fun)
  .tempssm_check_state_available(res, state, missing_msg)

  scale <- 1
  if (scale_by_frequency) {
    scale <- frequency(res$temp_data)
  }

  if (ci) {
    return(
      .tempssm_state_ci_ts(
        res = res,
        state = state,
        output_name = output_name,
        ci_level = ci_level,
        scale = scale,
        message = ci_missing_msg
      )
    )
  }

  .tempssm_state_ts(res, state, scale = scale)
}


#' Extract the level component as a time series
#'
#' @param res An object of class \code{"tempssm"} returned by \code{tempssm()}.
#' @param ci Logical; if TRUE, pointwise confidence intervals are returned.
#' @param ci_level Numeric confidence level between 0 and 1 (default: 0.95).
#'
#' @return
#' A univariate \code{ts} object of the smoothed level component
#' (in degrees Celsius).
#' If \code{ci = TRUE}, a multivariate \code{ts} object with columns
#' \code{level}, \code{lwr}, and \code{upr} is returned.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(niigata_sst)
#' res <- tempssm(niigata_sst)
#' level_ts <- get_level_ts(res)
#' }
get_level_ts <- function(res, ci = FALSE, ci_level = 0.95) {
  .tempssm_extract_state_ts(
    res = res,
    state = "level",
    output_name = "level",
    ci = ci,
    ci_level = ci_level,
    fun = "get_level_ts",
    missing_msg = "Level component not found in the smoothing results.",
    ci_missing_msg = "Level component not found in confidence intervals."
  )
}


#' Extract the smoothed drift (slope) component as a time series
#'
#' @inheritParams get_level_ts
#'
#' @details
#' The drift component is scaled to represent change per year.
#' For example, monthly data (frequency = 12) are multiplied by 12.
#'
#' @return
#' A univariate \code{ts} object of the smoothed drift component
#' (in degrees Celsius per year).
#' If \code{ci = TRUE}, a multivariate \code{ts} object with columns
#' \code{drift}, \code{lwr}, and \code{upr} is returned.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(niigata_sst)
#' res <- tempssm(niigata_sst)
#' drift <- get_drift_ts(res)
#' }
get_drift_ts <- function(res, ci = FALSE, ci_level = 0.95) {
  .tempssm_extract_state_ts(
    res = res,
    state = "slope",
    output_name = "drift",
    ci = ci,
    ci_level = ci_level,
    fun = "get_drift_ts",
    missing_msg = "Drift (slope) component not found in smoothing results.",
    ci_missing_msg = "Slope component not found in confidence intervals.",
    scale_by_frequency = TRUE
  )
}


#' Extract the smoothed seasonal component as a time series
#'
#' @inheritParams get_level_ts
#'
#' @details
#' The seasonal component represents recurrent intra-year variability
#' captured by seasonal dummy state components in the state space model.
#'
#' @return
#' A univariate \code{ts} object of the smoothed seasonal component
#' (in degrees Celsius).
#' If \code{ci = TRUE}, a multivariate \code{ts} object with columns
#' \code{season}, \code{lwr}, and \code{upr} is returned.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(niigata_sst)
#' res <- tempssm(niigata_sst)
#' season_ts <- get_season_ts(res)
#' }
get_season_ts <- function(res, ci = FALSE, ci_level = 0.95) {
  .tempssm_check_accessor_input(res, "get_season_ts")
  if (!res$use_season) {
    stop("Seasonal component is not included in the model.", call. = FALSE)
  }

  .tempssm_extract_state_ts(
    res = res,
    state = "sea_dummy1",
    output_name = "season",
    ci = ci,
    ci_level = ci_level,
    fun = "get_season_ts",
    missing_msg = "Seasonal component not found in the smoothing results.",
    ci_missing_msg = "Seasonal component not found in confidence intervals."
  )
}


#' Extract the smoothed first autoregressive component (AR1) as a time series
#'
#' @inheritParams get_level_ts
#'
#' @details
#' The AR1 component represents short-term autocorrelated deviations
#' from the level and seasonal structure.
#'
#' @return
#' A univariate \code{ts} object of the smoothed AR1 component
#' (in degrees Celsius).
#' If \code{ci = TRUE}, a multivariate \code{ts} object with columns
#' \code{ar1}, \code{lwr}, and \code{upr} is returned.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(niigata_sst)
#' res <- tempssm(niigata_sst)
#' ar1_ts <- get_ar1_ts(res)
#' }
get_ar1_ts <- function(res, ci = FALSE, ci_level = 0.95) {
  .tempssm_extract_state_ts(
    res = res,
    state = "arima1",
    output_name = "ar1",
    ci = ci,
    ci_level = ci_level,
    fun = "get_ar1_ts",
    missing_msg = paste0(
      "First autoregressive component (AR1) not found ",
      "in the smoothing results."
    ),
    ci_missing_msg = paste0(
      "First Autoregressive (AR1) component not found ",
      "in confidence intervals."
    )
  )
}


#' Extract estimated parameters in the fitted models
#'
#' @inheritParams get_level_ts
#'
#' @return A \code{list} object of the estimated parameters.
#'
#' @examples
#' \dontrun{
#' data(niigata_sst)
#' res <- tempssm(niigata_sst)
#' params <- get_tempssm_params(res)
#' }
#' @export
get_tempssm_params <- function(res) {
  if (!inherits(res, "tempssm")) {
    stop("`res` must be an object of class 'tempssm' ",
      "for get_tempssm_params().",
      call. = FALSE
    )
  }

  if (is.null(res$fit) || is.null(res$fit$optim.out) ||
      is.null(res$fit$optim.out$par)) {
    stop("Fitted parameters are not available in the provided tempssm object.",
      call. = FALSE
    )
  }

  model <- res$model
  pars <- res$fit$optim.out$par
  ar_order <- res$ar_order
  use_season <- res$use_season

  if (use_season) {
    ar_idx <- 3:(2 + ar_order)
    var_idx <- 3 + ar_order
    H_idx <- 4 + ar_order

    Q_season_est <- exp(pars[2])
  } else {
    ar_idx <- 2:(1 + ar_order)
    var_idx <- 2 + ar_order
    H_idx <- 3 + ar_order

    Q_season_est <- NA
  }

  params <- list(
    H = exp(pars[H_idx]), # observed error
    Q_trend = exp(pars[1]), # process error for level component
    Q_season = Q_season_est, # process error for seasonal component
    Q_ar = exp(pars[var_idx]), # process error for AR
    ARs = .tempssm_transform_ar(pars[ar_idx]) # the AR(s) coefficient
  )

  return(params)
}


#' Extract coefficients of exogenous variables with confidence intervals
#'
#' Extracts estimated regression coefficients for exogenous variable(s)
#' included in a \code{tempssm} model, together with confidence intervals
#' based on Kalman smoothing results.
#'
#' If the fitted model does not include exogenous variables,
#' the function returns \code{NULL}.
#'
#' @inheritParams get_level_ts
#'
#' @return
#' A \code{data.frame} with the following columns:
#' \describe{
#'   \item{Variable}{Name of the exogenous variable}
#'   \item{Coefficient}{Estimated regression coefficient}
#'   \item{lwr}{Lower bound of the confidence interval}
#'   \item{upr}{Upper bound of the confidence interval}
#' }
#' Returns \code{NULL} if no exogenous variables are included in the model.
#'
#' @examples
#' \dontrun{
#' data(niigata_sst)
#' data(pdo)
#  niigata_sst_pdo <- ts.intersect(niigata_sst,pdo)
#' colnames(niigata_sst_pdo) <- c("Temp", "PDO")
#' niigata_sst_common <- niigata_sst_pdo[, "Temp"]
#' pdo_common <- niigata_sst_pdo[, "PDO"]
#' pdo_common <- set_ts_name(nao_common, label = "PDO")
#' res <- ssm(temp_data = niigata_sst_common, exo_data = pdo_common)
#' get_exo_coef_ci(res)
#' }
#'
#' @importFrom utils head
#' @export
get_exo_coef <- function(res, ci_level = 0.95) {
  if (!inherits(res, "tempssm")) {
    stop("`res` must be an object of class 'tempssm' for get_exo_coef().",
      call. = FALSE
    )
  }

  if (!is.numeric(ci_level) || length(ci_level) != 1 ||
    ci_level <= 0 || ci_level >= 1) {
    stop(
      "`ci_level` must be a numeric value between 0 and 1 for get_exo_coef().",
      call. = FALSE
    )
  }

  ## ---- handle non-convergence or no exogenous variables ----
  if (isFALSE(res$converged)) {
    return(NULL)
  }

  exo_vars <- res$state_map$exogenous
  if (is.null(exo_vars) || length(exo_vars) == 0) {
    return(NULL)
  }

  ## ---- extract states ----
  kfs <- res$kfs
  alpha_hat <- kfs$alphahat
  state_names <- colnames(alpha_hat)

  exo_idx <- match(exo_vars, state_names)
  if (anyNA(exo_idx)) {
    stop(
      "Exogenous states listed in `state_map` were not found in model states.",
      call. = FALSE
    )
  }

  beta_hat <- alpha_hat[, exo_idx, drop = FALSE]

  ## ---- confidence intervals ----
  ci_all <- stats::confint(kfs, level = ci_level)
  ci_exo <- ci_all[exo_idx]

  ci_mat <- do.call(
    rbind,
    lapply(ci_exo, function(x) x[1, c("lwr", "upr")])
  )

  data.frame(
    Variable    = exo_vars,
    Coefficient = as.numeric(beta_hat[1, , drop = TRUE]),
    lwr         = ci_mat[, "lwr"],
    upr         = ci_mat[, "upr"],
    row.names   = NULL
  )
}
