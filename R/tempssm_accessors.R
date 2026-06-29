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


#' Check whether a confidence level is valid
#'
#' @param ci_level Numeric confidence level between 0 and 1.
#'
#' @return A logical scalar.
#' @noRd
.tempssm_is_valid_ci_level <- function(ci_level) {
  if (!is.numeric(ci_level) || length(ci_level) != 1L) {
    return(FALSE)
  }

  all(c(
    !anyNA(ci_level),
    is.finite(ci_level),
    ci_level > 0,
    ci_level < 1
  ))
}


#' Check confidence interval arguments for tempssm accessors
#'
#' @inheritParams .tempssm_is_valid_ci_level
#' @param ci Logical; if \code{TRUE}, confidence intervals are requested.
#' @param fun Character scalar naming the calling function.
#'
#' @return Invisibly returns \code{NULL}.
#' @noRd
.tempssm_check_accessor_ci <- function(ci, ci_level, fun) {
  valid_ci <- all(c(
    is.logical(ci),
    length(ci) == 1L,
    !anyNA(ci)
  ))
  if (!valid_ci) {
    stop(
      "`ci` must be a non-missing logical scalar for ",
      fun,
      "().",
      call. = FALSE
    )
  }

  if (!ci) {
    return(invisible(NULL))
  }

  if (!.tempssm_is_valid_ci_level(ci_level)) {
    stop(
      "`ci_level` must be a numeric value between 0 and 1 for ",
      fun,
      "().",
      call. = FALSE
    )
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

  .extract_tempssm_params(res)
}


#' Resolve exogenous states in smoothing results
#'
#' @param res A validated object of class \code{"tempssm"}.
#'
#' @return \code{NULL} when coefficients are unavailable, otherwise a named
#'   list containing exogenous variable names, smoothing results, smoothed
#'   states, and state positions.
#'
#' @noRd
.tempssm_exo_state_info <- function(res) {
  if (isFALSE(res$converged)) {
    return(NULL)
  }

  exo_vars <- res$state_map$exogenous
  if (length(exo_vars) == 0L) {
    return(NULL)
  }

  kfs <- res$kfs
  alpha_hat <- kfs$alphahat
  if (is.null(kfs) || is.null(alpha_hat) || !is.matrix(alpha_hat)) {
    stop(
      "Smoothing results for exogenous states are not available.",
      call. = FALSE
    )
  }

  exo_idx <- match(exo_vars, colnames(alpha_hat))
  if (anyNA(exo_idx)) {
    stop(
      "Exogenous states listed in `state_map` were not found in model states.",
      call. = FALSE
    )
  }

  list(
    variables = exo_vars,
    kfs = kfs,
    alpha_hat = alpha_hat,
    indices = exo_idx
  )
}


#' Extract the first confidence bounds for one exogenous state
#'
#' @param state_ci Confidence intervals for one smoothed state.
#'
#' @return A numeric vector containing lower and upper bounds.
#'
#' @noRd
.tempssm_first_exo_ci_bounds <- function(state_ci) {
  valid_shape <- all(c(
    is.matrix(state_ci) || is.data.frame(state_ci),
    NROW(state_ci) >= 1L,
    all(c("lwr", "upr") %in% colnames(state_ci))
  ))
  if (!valid_shape) {
    stop(
      "Confidence intervals for all exogenous states are not available.",
      call. = FALSE
    )
  }

  bounds <- unlist(
    state_ci[1L, c("lwr", "upr"), drop = TRUE],
    use.names = FALSE
  )
  if (!is.numeric(bounds) || length(bounds) != 2L) {
    stop(
      "Confidence intervals for all exogenous states must be numeric.",
      call. = FALSE
    )
  }

  as.numeric(bounds)
}


#' Extract confidence intervals for exogenous states by KFAS position
#'
#' KFAS may label exogenous states differently in code{alphahat} and
#' code{confint()} output. The state positions resolved from code{alphahat}
#' are therefore applied to the confidence interval list.
#'
#' @param kfs Kalman filtering and smoothing results.
#' @param indices Integer positions of exogenous states.
#' @param ci_level Numeric confidence level between 0 and 1.
#'
#' @return A numeric matrix with columns code{lwr} and code{upr}.
#'
#' @noRd
.tempssm_exo_ci_matrix <- function(kfs, indices, ci_level) {
  ci_all <- stats::confint(kfs, level = ci_level)
  if (!is.list(ci_all)) {
    stop(
      "Confidence intervals for exogenous states are not available.",
      call. = FALSE
    )
  }

  ci_exo <- ci_all[indices]
  ci_values <- vapply(
    ci_exo,
    .tempssm_first_exo_ci_bounds,
    numeric(2)
  )
  ci_mat <- t(ci_values)
  dimnames(ci_mat) <- list(NULL, c("lwr", "upr"))

  ci_mat
}


#' Construct an exogenous coefficient result
#'
#' @param state_info Exogenous state information returned by
#'   \code{.tempssm_exo_state_info()}.
#' @param ci_mat Numeric confidence interval matrix returned by
#'   \code{.tempssm_exo_ci_matrix()}.
#'
#' @return A data frame containing coefficient estimates and intervals.
#'
#' @noRd
.new_exo_coef_result <- function(state_info, ci_mat) {
  beta_hat <- state_info$alpha_hat[
    1L,
    state_info$indices,
    drop = TRUE
  ]

  data.frame(
    Variable = state_info$variables,
    Coefficient = as.numeric(beta_hat),
    lwr = ci_mat[, "lwr"],
    upr = ci_mat[, "upr"],
    row.names = NULL
  )
}


#' Extract coefficients of exogenous variables with confidence intervals
#'
#' Extracts estimated regression coefficients for exogenous variable(s)
#' included in a \code{tempssm} model, together with confidence intervals
#' based on Kalman smoothing results.
#'
#' If the fitted model did not converge or does not include exogenous
#' variables, the function returns \code{NULL}.
#'
#' @inheritParams get_level_ts
#'
#' @details
#' Exogenous coefficients are represented as static regression states. The
#' returned estimates and confidence limits are taken from the first smoothed
#' time point because these states are constant over time.
#'
#' @return
#' A \code{data.frame} with the following columns:
#' \describe{
#'   \item{Variable}{Name of the exogenous variable}
#'   \item{Coefficient}{Estimated regression coefficient}
#'   \item{lwr}{Lower bound of the confidence interval}
#'   \item{upr}{Upper bound of the confidence interval}
#' }
#' Returns \code{NULL} if the model did not converge or no exogenous variables
#' are included.
#'
#' @examples
#' \dontrun{
#' data(niigata_sst)
#' data(pdo)
#' common_data <- ts.intersect(niigata_sst, pdo)
#' temp_data <- common_data[, "niigata_sst"]
#' exo_data <- common_data[, "pdo", drop = FALSE]
#' res <- tempssm(temp_data = temp_data, exo_data = exo_data)
#' get_exo_coef(res)
#' }
#'
#' @export
get_exo_coef <- function(res, ci_level = 0.95) {
  .tempssm_check_accessor_input(res, "get_exo_coef")
  .tempssm_check_accessor_ci(TRUE, ci_level, "get_exo_coef")

  state_info <- .tempssm_exo_state_info(res)
  if (is.null(state_info)) {
    return(NULL)
  }

  ci_mat <- .tempssm_exo_ci_matrix(
    kfs = state_info$kfs,
    indices = state_info$indices,
    ci_level = ci_level
  )

  .new_exo_coef_result(state_info, ci_mat)
}
