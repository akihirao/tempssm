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
#' res <- lgssm(fuji_temp)
#' plot(res$data)
#' }
lgssm <- function(ts_data, 
                  inits = NULL,
                  exogenous = FALSE) {

  ## ---- Input checks ---------------------------------------------------
  if (!inherits(ts_data, "ts")) {
    stop("ts_data must be a 'ts' object.")
  }

  if (frequency(ts_data) != 12) {
    stop("ts_data must be a monthly time series (frequency = 12).")
  }

  data_df <- as.data.frame(ts_data)
  
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
  
  
  ## Handle univariate / multivariate ts
  # model without exogenous variables
  if(!(exogenous)){ 
    #if (nrow(data_df)==1) {
    #  y <- data_df
    if (is.null(dim(ts_data))) {
      y <- ts_data
    } else {
      if (!"Temp" %in% colnames(ts_data)) {
        stop("For multivariate ts objects, a column named 'Temp' is required.")
      }
        y <- ts_data[, "Temp"]
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
    
  
  #========================================================   
  #======================================================== 
  # model with an exogenous variable
  } else { # 
    variables_lab <- colnames(ts_data)
    if (!"Temp" %in% variables_lab) {
      stop("For multivariate ts objects, a column named 'Temp' is required.")
    }
    num_exo_variable <- length(variables_lab) - 1
    #print(num_exo_variable)
    
    if (num_exo_variable >= 2){
      stop("Sorry. More than two exogenous variables can not be conducted in this packge")
    }
    y <- ts_data[, "Temp"]
    exo_lab <- variables_lab[variables_lab != "Temp"]
    
    exo1 <- ts_data[, exo_lab]
    
    reg_formula <- as.formula(
      paste("~ ", paste(exo_lab, collapse = " + "))
    )
    
    #browser()
    #print(reg_formula)
    
    ## ---- Model definition -----------------------------------------------
    build_ssm <- SSModel(
      H = NA,
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
        ) +
        SSMregression( #exogenous variable
          ~ exo1,
          Q = 0
        )#,
      #data = data_df
    )
    
    ## ---- Parameter update function --------------------------------------
    update_func <- function(pars, model) {
      return(
        SSModel(
          H = exp(pars[6]),
          y ~
            SSMtrend(degree = 2,
                     Q = c(list(0), list(exp(pars[1])))) +
            SSMseasonal(
              sea.type = "dummy",
              period = 12,
              Q = exp(pars[2])) +
            SSMarima(
              ar = artransform(pars[3:4]),
              d = 0,
              Q = exp(pars[5])) +
            SSMregression( #exogenous variable
              ~ exo1,
              Q = 0)#,
        #data = data_df
        )
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
    
      
    
  }  


  ## ---- Output ---------------------------------------------------------
  out <- list(
    model = fit2$model,
    fit   = fit2,
    kfs   = kfs,
    data  = ts_data,
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




#' Extract CI of coefficient for exogenous variable(s)
#'
#' @param res An object of class \code{"ThermoSSM"} returned by \code{lgssm()}.
#'
#' @return A \code{data.frame} object of CI for exogenous variable(s).
#'
#' @export

extract_regression_ci <- function(res, level = 0.95) {
  
  if (!inherits(res, "ThermoSSM")) {
    stop("Input must be a 'ThermoSSM' object.")
  }
  
  model <- res$model
  kfs   <- res$kfs
  
  # 外生変数名（Temp 以外）
  exo_names <- setdiff(colnames(res$data), "Temp")
  
  if (length(exo_names) == 0) {
    stop("No exogenous variables found in the model.")
  }
  
  p <- length(exo_names)
  Tt <- ncol(kfs$alphahat)
  
  # 回帰係数は状態ベクトルの末尾
  idx <- (nrow(kfs$alphahat) - p + 1):nrow(kfs$alphahat)
  
  beta_hat <- kfs$alphahat[idx, Tt]
  beta_var <- diag(kfs$V[idx, idx, Tt])
  
  alpha <- 1 - level
  z <- qnorm(1 - alpha / 2)
  
  out <- data.frame(
    variable = exo_names,
    estimate = beta_hat,
    std_error = sqrt(beta_var),
    lower = beta_hat - z * sqrt(beta_var),
    upper = beta_hat + z * sqrt(beta_var)
  )
  
  return(out)
}

