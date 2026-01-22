#' Fit a linear Gaussian state-space model to temperature time series
#'
#' This function estimates a linear Gaussian state-space model (LGSSM)
#' for monthly temperature time series using Kalman filtering and smoothing.
#' Both univariate and multivariate \code{ts} objects are supported.
#' If multivariate, a column named \code{"Temp"} is used.
#'
#' @param ts_data Monthly temperature time series of class \code{ts}.
#'   The time series must have a frequency of 12 (monthly data).
#'   May be univariate or multivariate.
#'
#' @param inits Optional numeric vector of initial parameter values.
#'   If \code{NULL}, default values are used.
#'
#' @return An object of class \code{"ThermoSSM"}, a named list containing:
#' \describe{
#'   \item{model}{Fitted \code{SSModel} object.}
#'   \item{fit}{Result of \code{fitSSM}.}
#'   \item{kfs}{Kalman filtering and smoothing results from \code{KFS}.}
#'   \item{data}{Temperature time series used for estimation.}
#'   \item{call}{Matched function call.}
#' }
#'
#' @import KFAS
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- lgssm(fuji_temp)
#' plot(fit$data)
#' }
lgssm <- function(ts_data, inits = NULL) {

  ## ---- Input checks ---------------------------------------------------
  if (!inherits(ts_data, "ts")) {
    stop("ts_data must be a 'ts' object.")
  }

  if (frequency(ts_data) != 12) {
    stop("ts_data must be a monthly time series (frequency = 12).")
  }

  ## Handle univariate / multivariate ts
  if (is.null(dim(ts_data))) {
    y <- ts_data
  } else {
    if (!"Temp" %in% colnames(ts_data)) {
      stop("For multivariate ts objects, a column named 'Temp' is required.")
    }
    y <- ts_data[, "Temp"]
  }

  #browser()
  #colnames(y) <- "Temp"

  ## ---- Default initial values -----------------------------------------
  if (is.null(inits)) {
    # log-variances and AR parameters (on unconstrained scale)
    inits <- c(
      -13,  # trend variance (log)
      -7,   # seasonal variance (log)
      0.9,  # AR(1)
      -0.1, # AR(2)
      -0.3, # AR noise variance (log)
      -5    # observation variance (log)
    )
  }

  if (!is.numeric(inits) || length(inits) != 6) {
    stop("inits must be a numeric vector of length 6.")
  }

  ## ---- Model definition -----------------------------------------------
  build_ssm <- SSModel(
    y ~
      SSMtrend(
        degree = 2,
        Q = c(list(0), list(NA))
      ) +
      SSMseasonal(
        sea.type = "dummy",
        period = 12,
        Q = NA
      ) +
      SSMarima(
        ar = c(0, 0),
        d = 0,
        Q = 0
      ),
    H = NA
  )

  ## ---- Parameter update function --------------------------------------
  update_func <- function(pars, model) {
    model <- SSModel(
      y ~
        SSMtrend(degree = 2,
                 Q = c(list(0), list(exp(pars[1])))) +
        SSMseasonal(
          sea.type = "dummy",
          period = 12,
          Q = exp(pars[2])
        ) +
        SSMarima(
          ar = artransform(pars[3:4]),
          d = 0,
          Q = exp(pars[5])
        ),
      H = exp(pars[6])
    )
  }


  ## ---- Optimization (two-step) ----------------------------------------
  fit1 <- fitSSM(
    build_ssm,
    inits  = inits,
    updatefn = update_func,
    method = "Nelder-Mead",
    control = list(maxit = 5000, reltol = 1e-16)
  )

  fit2 <- fitSSM(
    build_ssm,
    inits  = fit1$optim.out$par,
    updatefn = update_func,
    method = "BFGS",
    control = list(maxit = 5000, reltol = 1e-16)
  )

  ## ---- Kalman filtering & smoothing -----------------------------------
  kfs <- KFS(
    fit2$model,
    filtering = c("state", "mean"),
    smoothing = c("state", "mean", "disturbance")
  )

  ## ---- Output ---------------------------------------------------------
  out <- list(
    model = fit2$model,
    fit   = fit2,
    kfs   = kfs,
    data  = y,
    call  = match.call()
  )

  class(out) <- "ThermoSSM"
  return(out)
}




#' Extract smoothed level component as a time series
#'
#' @param res An object of class \code{"ThermoSSM"} returned by \code{lgssm()}.
#'
#' @return A \code{ts} object of the smoothed level component.
#'
#' @export
extract_level_ts <- function(res) {

  if (!inherits(res, "ThermoSSM")) {
    stop("Input must be a ThermoSSM object.")
  }

  if (is.null(res$kfs$alphahat) || !"level" %in% colnames(res$kfs$alphahat)) {
    stop("Level component not found in the smoothing results.")
  }

  res$kfs$alphahat[, "level"]
}





#' Extract smoothed drift (slope) component as a time series
#'
#' @param res An object of class \code{"ThermoSSM"} returned by \code{lgssm()}.
#'
#' @return A \code{ts} object of the smoothed drift component.
#'
#' @export
extract_drift_ts <- function(res) {

  if (!inherits(res, "ThermoSSM")) {
    stop("Input must be a ThermoSSM object.")
  }

  if (is.null(res$kfs$alphahat) || !"slope" %in% colnames(res$kfs$alphahat)) {
    stop("Drift (slope) component not found in the smoothing results.")
  }

  res$kfs$alphahat[, "slope"]
}




#' Extract estimated parameters in the fitted models
#'
#' @param res An object of class \code{"ThermoSSM"} returned by \code{lgssm()}.
#'
#' @return A \code{list} object of the estimated parameters.
#'
#' @export

extract_param <- function(res){
  pars <- res$fit$optim.out$par

  params <- c(Q_trend  = exp(pars[1]), # 年トレンドの大きさ
              Q_season = exp(pars[2]), # 季節トレンドの大きさ
              AR1      = KFAS::artransform(pars[3:4])[1], # 1次のARの大きさ
              AR2      = KFAS::artransform(pars[3:4])[2], # 2次のARの大きさ
              Q_ar     = exp(pars[5]), # 短期変動の揺らぎ
              H        = exp(pars[6]) # 観察誤差の大きさ
              )
  return(params)
}

