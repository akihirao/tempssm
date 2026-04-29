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

  if (!inherits(res, "tempssm")) {
    stop("`res` must be an object of class 'tempssm'.", call. = FALSE)
  }

  if (ci) {
    if (!is.numeric(ci_level) || length(ci_level) != 1 ||
        ci_level <= 0 || ci_level >= 1) {
      stop("`ci_level` must be a numeric value between 0 and 1.",
           call. = FALSE)
    }
  }

  if (is.null(res$kfs$alphahat) || !"level" %in% colnames(res$kfs$alphahat)) {
    stop("Level component not found in the smoothing results.", call. = FALSE)
  }

  freq <- frequency(res$data_temp)
  
  level <- ts(
    res$kfs$alphahat[, "level"],
    start = start(res$data_temp),
    frequency = freq
  )

  if(ci){
    ci_obj <- stats::confint(res$kfs, level = ci_level)

    if (!"level" %in% names(ci_obj)) {
      stop("Level component not found in confidence intervals.",
           call. = FALSE)
    }

    level <- ts(
      cbind(
        level = level,
        lwr = ci_obj$level[, "lwr"],
        upr = ci_obj$level[, "upr"]
      ),
      start = start(res$data_temp),
      frequency = freq
    )
  }
  level
}



#' Extract the smoothed drift (slope) component as a time series
#'
#' @param res An object of class \code{"tempssm"} returned by \code{tempssm()}.
#' @param ci Logical; if TRUE, pointwise confidence intervals are returned.
#' @param ci_level Numeric confidence level between 0 and 1 (default: 0.95).
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
  
  if (!inherits(res, "tempssm")) {
    stop("`res` must be an object of class 'tempssm'.", call. = FALSE)
  }
  
  if (ci) {
    if (!is.numeric(ci_level) || length(ci_level) != 1 ||
        ci_level <= 0 || ci_level >= 1) {
      stop("`ci_level` must be a numeric value between 0 and 1.",
           call. = FALSE)
    }
  }
  
  if (is.null(res$kfs$alphahat) || !"slope" %in% colnames(res$kfs$alphahat)) {
    stop("Drift (slope) component not found in smoothing results.",
         call. = FALSE)
  }
  
  freq <- frequency(res$data_temp)
  scale <- freq
  
  drift <- ts(
    res$kfs$alphahat[, "slope"] * scale, #scaled per year
    start = start(res$data_temp),
    frequency = freq
  )
  
  if (ci) {
    
    ci_obj <- stats::confint(res$kfs, level = ci_level)
    
    if (!"slope" %in% names(ci_obj)) {
      stop("Slope component not found in confidence intervals.",
           call. = FALSE)
    }
    
    drift <- ts(
      cbind(
        drift = drift,
        lwr = ci_obj$slope[, "lwr"] * scale,#scaled per year 
        upr = ci_obj$slope[, "upr"] * scale #scaled per year
      ),
      start = start(res$data_temp),
      frequency = freq
    )
  }
  
  drift
}



#' Extract the smoothed seasonal component as a time series
#'
#' @param res An object of class \code{"tempssm"} returned by \code{tempssm()}.
#' @param ci Logical; if TRUE, pointwise confidence intervals are returned.
#' @param ci_level Numeric confidence level between 0 and 1 (default: 0.95).
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
  
  if (!inherits(res, "tempssm")) {
    stop("`res` must be an object of class 'tempssm'.", call. = FALSE)
  }
  
  if (ci) {
    if (!is.numeric(ci_level) || length(ci_level) != 1 ||
        ci_level <= 0 || ci_level >= 1) {
      stop("`ci_level` must be a numeric value between 0 and 1.",
           call. = FALSE)
    }
  }
  
  if (is.null(res$kfs$alphahat) || ! "sea_dummy1" %in% colnames(res$kfs$alphahat)) {
    stop("Seasonal component not found in the smoothing results.", call. = FALSE)
  }
  
  freq <- frequency(res$data_temp)

  season <- ts(
    res$kfs$alphahat[, "sea_dummy1"],
    start = start(res$data_temp),
    frequency = freq
  )
  
  if(ci){
    ci_obj <- stats::confint(res$kfs, level = ci_level)
    
    if (!"sea_dummy1" %in% names(ci_obj)) {
      stop("Seasonal component not found in confidence intervals.",
           call. = FALSE)
    }
    
    season <- ts(
      cbind(
        season = season,
        lwr = ci_obj$sea_dummy1[, "lwr"],
        upr = ci_obj$sea_dummy1[, "upr"]
      ),
      start = start(res$data_temp),
      frequency = freq
    )
  }
  season
}



#' Extract the smoothed first autoregressive component (AR1) as a time series
#'
#' @param res An object of class \code{"tempssm"} returned by \code{tempssm()}.
#' @param ci Logical; if TRUE, pointwise confidence intervals are returned.
#' @param ci_level Numeric confidence level between 0 and 1 (default: 0.95).
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
  
  if (!inherits(res, "tempssm")) {
    stop("`res` must be an object of class 'tempssm'.", call. = FALSE)
  }
  
  if (ci) {
    if (!is.numeric(ci_level) || length(ci_level) != 1 ||
        ci_level <= 0 || ci_level >= 1) {
      stop("`ci_level` must be a numeric value between 0 and 1.",
           call. = FALSE)
    }
  }
  
  if (is.null(res$kfs$alphahat) || ! "arima1" %in% colnames(res$kfs$alphahat)) {
    stop("First autoregressive component (AR1) not found in the smoothing results.", call. = FALSE)
  }
  
  freq <- frequency(res$data_temp)
  
  ar1 <- ts(
    res$kfs$alphahat[, "arima1"],
    start = start(res$data_temp),
    frequency = freq
  )
  
  if(ci){
    ci_obj <- stats::confint(res$kfs, level = ci_level)
    
    if (!"arima1" %in% names(ci_obj)) {
      stop("First Autoregressive (AR1) component not found in confidence intervals.",
           call. = FALSE)
    }
    
    ar1 <- ts(
      cbind(
        ar1 = ar1,
        lwr = ci_obj$arima1[, "lwr"],
        upr = ci_obj$arima1[, "upr"]
      ),
      start = start(res$data_temp),
      frequency = freq
    )
  }
  ar1
}



#' Extract estimated parameters in the fitted models
#'
#' @param res An object of class \code{"tempssm"} returned by \code{ssm()}.
#'
#' @return A \code{list} object of the estimated parameters.
#'
#' @examples
#' \dontrun{
#' data(niigata_sst)
#' res <- tempssm(niigata_sst)
#' params <- get_params(res)
#' }
#' @export
get_params <- function(res){

  if (!inherits(res, "tempssm")) {
    stop("`res` must be an object of class 'tempssm.'", call. = FALSE)
  }
  
  model <- res$model
  pars <- res$fit$optim.out$par
  ar_order <- res$ar_order
  use_season <- res$use_season

  if(use_season){
    ar_idx  <- 3:(2 + ar_order)
    var_idx <- 3 + ar_order
    H_idx   <- 4 + ar_order

    Q_season_est <- exp(pars[2])

  }else{
    ar_idx  <- 2:(1 + ar_order)
    var_idx <- 2 + ar_order
    H_idx   <- 3 + ar_order

    Q_season_est <- NA
  }

  params = list(
    H = exp(pars[H_idx]), # observed error
    Q_trend  = exp(pars[1]),# process error for level component
    Q_season = Q_season_est, # process error for seasonal component
    Q_ar = exp(pars[var_idx]),# process error for AR
    ARs = KFAS::artransform(pars[ar_idx])# the AR(s) coefficient
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
#' @param res An object of class \code{"tempssm"} returned by \code{ssm()}.
#'
#' @param ci_level Numeric confidence level between 0 and 1 (default: 0.95).
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
#' niigata_sst_common <- niigata_sst_pdo[,"Temp]
#' pdo_common <- niigata_sst_pdo[,"PDO]
#' pdo_common <- set_ts_name(nao_common, label = "PDO")
#' res <- ssm(temp_data = niigata_sst_common,exo_data = pdo_common)
#' extract_exo_coef_ci(res)
#' }
#'
#' @importFrom utils head
#' @export
get_exo_coef = function(res, ci_level = 0.95) {
  
  if (!inherits(res, "tempssm")) {
    stop("`res` must be an object of class 'tempssm.'", call. = FALSE)
  }
  
  if (!is.numeric(ci_level) || ci_level <= 0 || ci_level >= 1) {
    stop("'CI level' must be a numeric value between 0 and 1.", call. = FALSE)
  }
  
  exo_vars <- colnames(res$data_exogenous)
  kfs      <- res$kfs
  
  # No exogenous variables
  if (is.null(exo_vars)) {
    return(NULL)
  }
  
  n_exo <- length(exo_vars)
  
  # Extract smoothed state estimates
  alpha_hat <- kfs$alphahat
  
  # Exogenous coefficients are assumed to be the first states
  beta_hat <- alpha_hat[1, seq_len(n_exo), drop = FALSE]
  
  # Confidence intervals
  ci_obj <- stats::confint(kfs, level = ci_level)[seq_len(n_exo)]
  ci_mat <- do.call(rbind, lapply(ci_obj, head, n = 1))
  
  data.frame(
    Variable    = exo_vars,
    Coefficient = as.numeric(beta_hat),
    lwr         = ci_mat[, "lwr"],
    upr         = ci_mat[, "upr"],
    row.names   = NULL
  )
}


