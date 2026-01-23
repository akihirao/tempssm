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
#' @param exo logical(1). If TRUE, fit a model that incorporates 
#'  exogenous variables, The default is FALSE.
#'  
#'  
#' @param exo_name Optional characteristic vector of names of  
#'   exogenous variables. If \code{NULL}, column names of ts_data of 
#'   class \code{ts}, excluding a column named \code{"Temp"}, are used. 
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
                  exo = FALSE,
                  exo_name = NULL) {

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
       -7,  # seasonal variance (log)
      0.9,  # AR(1)
     -0.1,  # AR(2)
     -0.3,  # AR noise variance (log)
       -5   # observation variance (log)
    )
  }
  
  if (!is.numeric(inits) || length(inits) != 6) {
    stop("inits must be a numeric vector of length 6.")
  }
  
  
  #========================================================
  # Handle univariate / multivariate ts
  
  # ------------------------------------------------------------------------
  # ------------------------------------------------------------------------
  # model without exogenous variables
  if(!(exo)){ 
    if (is.null(dim(ts_data))) {
      y <- ts_data
    } else {
      if (!"Temp" %in% colnames(ts_data)) {
        stop("For multivariate ts objects, a column named 'Temp' is required.")
      }
        y <- ts_data[, "Temp"]
    }
    
    exogenous_lab <- exo_name
    
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
    
   
  # ------------------------------------------------------------------------
  # ------------------------------------------------------------------------
  # model with an exogenous variable
  } else { # 
    variables_lab <- colnames(ts_data)
    if (!"Temp" %in% variables_lab) {
      stop("For multivariate ts objects, a column named 'Temp' is required.")
    }
    
    y <- ts_data[, "Temp"]
    
    if(is.null(exo_name)){ # case_when exo_name = NULL
      num_exo_variable <- length(variables_lab) - 1
      
      if (num_exo_variable == 0){
        stop("At least one exogenous variable(s) should be included in ts object")
      }
      
      exogenous_lab <- variables_lab[variables_lab != "Temp"]
      
      exogenous_ts <- ts_data[, exogenous_lab]
      exogenous_mat <- as.matrix(exogenous_ts)
    
    }else{ # case_when exo_name(s) is/are given
      
      if (!exo_name %in% variables_lab) {
        stop(paste0("Multivariate ts object should be included ",exo_name," as exogenous variables."))
      }
      
      exogenous_lab <- exo_name
      exogenous_ts <- ts_data[, exogenous_lab]
      exogenous_mat <- as.matrix(exogenous_ts)
    }
    
    
    ## ---- Model definition -----------------------------------------------
    build_ssm <- SSModel(
      H = NA,
      y ~ exogenous_mat +
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
        )
    )
    
    ## ---- Parameter update function --------------------------------------
    update_func <- function(pars, model) {
      return(
        SSModel(
          H = exp(pars[6]),
          y ~ exogenous_mat + 
            SSMtrend(degree = 2,
                     Q = c(list(0), list(exp(pars[1])))) +
            SSMseasonal(
              sea.type = "dummy",
              period = 12,
              Q = exp(pars[2])) +
            SSMarima(
              ar = artransform(pars[3:4]),
              d = 0,
              Q = exp(pars[5]))
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

  } # close of model with exogenous variable(s)  


  ## ---- Output ---------------------------------------------------------
  out <- list(
    model = fit2$model,
    fit   = fit2,
    kfs   = kfs,
    data  = ts_data,
    exogenous = exogenous_lab,
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




#' Extract coefficient of exogenous variables with confidence interval
#'
#' @param res An object of class \code{"ThermoSSM"} returned by \code{lgssm()}.
#'
#' @param level Percentage of confidence interval. The default value of 0.95.
#'
#' @return A \code{data.frame} object of CI for exogenous variable(s).
#'
#' @export

extract_exo_coef_ci <- function(res, 
                                level = 0.95) {
  
  if (!inherits(res, "ThermoSSM")) {
    stop("Input must be a 'ThermoSSM' object.")
  }
  
  model <- res$model
  kfs   <- res$kfs
  exo_variable <- res$exogenous
  
  if(is.null(exo_variable)){ # case when model without exogenous variable(s)
  
      beta_coef_ci <- NULL
      
  }else{ # case when model with exogenous variable(s)
  
    num_exo_variable <- length(exo_variable)
  
    alpha_hat <- kfs$alphahat
    colnames(alpha_hat) <- c(exo_variable,
                           "level",
                           "slope",
                           "sea_dummy1",
                           "sea_dummy2",
                           "sea_dummy3",
                           "sea_dummy4",
                           "sea_dummy5",
                           "sea_dummy6",
                           "sea_dummy7",
                           "sea_dummy8",
                           "sea_dummy9",
                           "sea_dummy10",
                           "sea_dummy11",
                           "arima1",
                           "arima2"
                           )
    if(num_exo_variable >=2){
      beta_coef <- alpha_hat[,exo_variable][1,1:num_exo_variable]
    }else{
      beta_coef <- alpha_hat[,exo_variable][1]
    }
    
    ci_general <- confint(kfs, level = 0.95)[1:num_exo_variable] %>% 
      lapply(head, n = 1) %>% 
      as.data.frame
  
    beta_coef_ci = data.frame(Variable=exo_variable,
                              Coefficient= beta_coef,
                              lwr = as.numeric(ci_general[seq(1, by=2, length.out = num_exo_variable)]),
                              upr = as.numeric(ci_general[seq(2, by=2, length.out = num_exo_variable)])
                              )
  }
  
  return(beta_coef_ci)
  
}

