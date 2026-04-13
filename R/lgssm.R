#' Fit a linear Gaussian state-space model to temperature time series
#'
#' This function estimates a linear Gaussian state-space model (LGSSM)
#' for monthly temperature time series using Kalman filtering and smoothing.
#' Both univariate and multivariate \code{ts} objects are supported.
#' If multivariate, a column named \code{"Temp"} is used.
#'
#' @param temp_data A temperature time series of class \code{ts}.
#'   The \code{ts} object must be univariant.
#'   The series can have any arbitrary frequency of 2 or higher.
#'   For example, a frequency of 12 represents a monthly \code{ts} object.
#'
#' @param exo_data A data set of exogenous variable(s) of class \code{ts}.
#'   The series may have any arbitrary frequency of 2 or higher,
#'   but it must be the same as that of \code{temp_data}.
#'   The default is \code{NULL} when fitting a model without exogenous variables.
#'
#' @param ar_order Integer specifying the order of the autoregressive (AR) 
#' component in the error structure (e.g., 2 for AR(2), 3 for AR(3)). 
#' Defaults to 2.
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
lgssm <- function(temp_data, 
                  exo_data = NULL,
                  ar_order = 2,
                  inits = NULL,
                  maxit = NULL,
                  reltol = NULL) {

  tryCatch(
    {
    ## ---- Input checks ---------------------------------------------------
    y <- ThermoSSM::check_temp_ts_lgssm(temp_data)
    freq <- frequency(y)

    ## ---- Default initial values -----------------------------------------
    if (is.null(inits)) {
      # log-variances and AR parameters (on unconstrained scale)
      ar_rep_length_minus_one <- ar_order -1
      ar_inits <- c(0.5,rep(0,ar_rep_length_minus_one))
    
      inits <- c(
        -13,  # trend variance (log)
        -7,  # seasonal variance (log)
        ar_inits,  # AR coefficients     
        -0.3,  # AR noise variance (log)
        -5   # observation variance (log)
      )
    }
  
    expected_len <- 2 + ar_order + 2
    if (!is.numeric(inits) != expected_len) {
      stop(paste("inits must be length ", expected_len))
    }
  
    ar_idx <- 3:(2 + ar_order)
    var_idx <- 3 + ar_order
    H_idx   <- 4 + ar_order
  
    ## ---- Default maxit value -----------------------------------------  
    if (is.null(maxit)) {
      maxit <- 5000
    }
  
    ## ---- Default reltol value -----------------------------------------  
    if (is.null(reltol)) {
      reltol <- 1e-16
    }
  
  
    #========================================================
    # Handle univariate / multivariate ts
  
    # ------------------------------------------------------------------------
    # ------------------------------------------------------------------------
    # model without exogenous variables
    if(is.null(exo_data)){ 
    
    exogenous_lab <- NULL

    ## ---- Model definition -----------------------------------------------
    build_ssm <- SSModel(
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
    
    ## ---- Parameter update function --------------------------------------
    update_func <- function(pars, model) {
      model <- SSModel(
        y ~
          SSMtrend(degree = 2,
                   Q = c(list(0), list(exp(pars[1])))) +
          SSMseasonal(
            sea.type = "dummy",
            period = freq,
            Q = exp(pars[2])
          ) +
          SSMarima(
            ar = artransform(pars[ar_idx]),
            d = 0,
            Q = exp(pars[var_idx])
          ),
        H = exp(pars[H_idx])
      )
    }


    # ------------------------------------------------------------------------
    # ------------------------------------------------------------------------
    # model with an exogenous variable
    } else { # 

      exo_data_checked <- check_exo_ts_lgssm(temp_data = temp_data,
                                             exo_data = exo_data)
    
      if (!inherits(exo_data, "ts")) {
        stop("temp_data must be a 'ts' object.")
      }
    
      if (frequency(exo_data) != freq) {
        stop("frequency of exo_data must be same that of temp_data.")
      }
    
      if (is.null(colnames(exo_data))) {
        stop("exo_data must have column name(s).")
      }
    
      exogenous_lab <- colnames(exo_data_checked)
      exogenous_mat <- as.matrix(exo_data_checked)


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
            period = freq,
            Q = NA
          ) +
          SSMarima(
            ar = rep(0, ar_order),
            d = 0,
            Q = NA
          )
      )
  

      ## ---- Parameter update function --------------------------------------
      update_func <- function(pars, model) {
        return(
          SSModel(
            H = exp(pars[H_idx]),
            y ~ exogenous_mat + 
              SSMtrend(degree = 2,
                       Q = c(list(0), 
                             list(exp(pars[1])))) +
              SSMseasonal(
                sea.type = "dummy",
                period = freq,
                Q = exp(pars[2])) +
              SSMarima(
                ar = artransform(pars[ar_idx]),
                d = 0,
                Q = exp(pars[var_idx]))
          )
        )
      }
    }# close case when model with exogenous variable(s)

    
    ## ---- Optimization (two-step) ----------------------------------------
    fit1 <- fitSSM(
      build_ssm,
      inits  = inits,
      updatefn = update_func,
      method = "Nelder-Mead",
      control = list(maxit = maxit, reltol = reltol)
    )
    
    fit2 <- fitSSM(
      build_ssm,
      inits  = fit1$optim.out$par,
      updatefn = update_func,
      method = "BFGS",
      control = list(maxit = maxit, reltol = reltol)
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
      data_temp  = temp_data,
      data_exogenous = exo_data,
      ar_order = ar_order,
      call  = match.call()
    )

    class(out) <- "ThermoSSM"
    return(out)
  },
  error = function(e) {
    message("Warning(s): NA is returned. ", e$message)
    NA
  }
  )# close tryCatch
}



#' Check ts object of temperature time series for applying \code{lgssm()}
#'
#' @param temp_data A temperature time series of class \code{ts}.
#'   The series can have any arbitrary frequency of 2 or higher.
#'   For example, a frequency of 12 represents a monthly \code{ts} object.
#'
#' @return A univariate \code{ts} object.
#'
#' @export
check_temp_ts_lgssm <- function(temp_data) {
  
  if (!inherits(temp_data, "ts")) {
    stop("The object 'temp_data' must be a 'ts' object.",
         call. = FALSE)
  }
  
  freq <- frequency(temp_data)
  if (freq <= 1) {
    stop("The procedure requires a ts object with frequency > 1.",
         call. = FALSE)
  }
  
  if (!is.null(dim(temp_data)) && NCOL(temp_data) != 1) {
    stop("The object 'temp_data' must be univariate.",
         call. = FALSE)
  }
  
  message("The ts object is univariate with frequency ", freq, ".")
  
  return(temp_data)
}


  
#' Check ts object of exogenous variable(s) for applying \code{lgssm()}
#'
#' @param temp_data A temperature time series of class \code{ts}.
#'   The series can have any arbitrary frequency of 2 or higher.
#'   For example, a frequency of 12 represents a monthly \code{ts} object.
#'
#' @param exo_data A data set of exogenous variable(s) of class \code{ts}.
#'   The series may have any arbitrary frequency of 2 or higher,
#'   but it must be the same as that of \code{temp_data}.
#'   The default is \code{NULL} when fitting a model without exogenous variables.
#'
#' @return A univariate or multivairate \code{ts} object.
#'
#' @export
check_exo_ts_lgssm <- function(temp_data, exo_data) {

  temp_data_checked <- ThermoSSM::check_temp_ts_lgssm(temp_data)

  temp_freq <- frequency(temp_data_checked)

  if (!inherits(exo_data, "ts")) {
    stop("The object 'exo_data' must be a 'ts' object.",
         call. = FALSE)
  }

  exo_freq <- frequency(exo_data)
  if (!(exo_freq == temp_freq)) {
    stop("Frequency of 'exo_data' must be same that of 'temp_data'.",
         call. = FALSE)
  }

  if (!(time(exo_freq) == time(temp_freq))) {
    stop("Time series of 'exo_data' must be same that of 'temp_data'.",
         call. = FALSE)
  }

  if (is.null(colnames(exo_data))) {
      stop("The object 'exo_data' must have column name(s).",
         call. = FALSE)
  }

  if ((dim(exo_data)[2]) > 1) {
    uni_multi <- "multivariate"
  }else{
    uni_multi <- "univariate"
  }
  
  message(paste0("The object 'exo_data' is ", uni_multi, " with frequency ", exo_freq, "."))

  return(exo_data)
}




#' Extract the smoothed level component as a time series
#'
#' @param res An object of class \code{"ThermoSSM"} returned by \code{lgssm()}.
#' @param ci Logical; should confidence intervals be returned? (default: FALSE)
#' @param ci_level Confidence level for intervals (default: 0.95).
#'
#' @return
#' A univariate \code{ts} object of the smoothed level component.
#' If \code{ci = TRUE}, a multivariate \code{ts} object with columns
#' \code{level}, \code{lwr}, and \code{upr} is returned.
#'
#' @export
extract_level_ts <- function(res, ci = FALSE, ci_level = 0.95) {

  if (!inherits(res, "ThermoSSM")) {
    stop("Input must be a ThermoSSM object.", call. = FALSE)
  }

  if (is.null(res$kfs$alphahat) || !"level" %in% colnames(res$kfs$alphahat)) {
    stop("Level component not found in the smoothing results.", call. = FALSE)
  }

  level_ts <- res$kfs$alphahat[, "level"]
  
  if(ci){
    ci_obj <- confint(res$kfs, level = ci_level)

    level_ts <- cbind(
      level = level_ts,
      lwr=ci_obj$level[,"lwr"],
      upr=ci_obj$level[,"upr"]
      )
  }
  return(level_ts)
}




#' Extract the smoothed drift (slope) component as a time series
#'
#' @param res An object of class \code{"ThermoSSM"} returned by \code{lgssm()}.
#' @param ci Logical; should confidence intervals be returned? (default: FALSE)
#' @param ci_level Confidence level for intervals (default: 0.95).
#'
#' @return
#' A univariate \code{ts} object of the smoothed drift component.
#' If \code{ci = TRUE}, a multivariate \code{ts} object with columns
#' \code{level}, \code{lwr}, and \code{upr} is returned.
#'
#' @export
extract_drift_ts <- function(res, ci = FALSE, ci_level = 0.95) {

  if (!inherits(res, "ThermoSSM")) {
    stop("Input must be a ThermoSSM object.", call. = FALSE)
  }

  if (is.null(res$kfs$alphahat) || !"slope" %in% colnames(res$kfs$alphahat)) {
    stop("Drift (slope) component not found in the smoothing results.", call. = FALSE)
  }

  drift_ts <- res$kfs$alphahat[, "slope"]

  if(ci){
  ci_obj <- confint(res$kfs, level = ci_level)

  drift_ts <- cbind(
    level = drift_ts,
    lwr=ci_obj$slope[,"lwr"],
    upr=ci_obj$slope[,"upr"]
    )
  }
  return(drift_ts)
}




#' Extract the Akaike Information Criterion (AIC)
#'
#' This function computes the Akaike Information Criterion (AIC)
#' for a fitted \code{ThermoSSM} model.
#'
#' @param res An object of class \code{"ThermoSSM"} returned by \code{lgssm()}.
#'
#' @return A numeric value representing the AIC of the fitted model.
#'
#' @export
extract_AIC <- function(res) {

  if (!inherits(res, "ThermoSSM")) {
    stop("Input must be a ThermoSSM object.", call. = FALSE)
  }

  k <- length(res$opt$par)

  # if the model includes exogenous variable(s)
  if (!is.null(res$data_exogenous)) {
    k <- k + ncol(res$data_exogenous)
  }

  loglik <- as.numeric(logLik(res$model))

  AIC <- -2 * loglik + 2 * k

  return(AIC)
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

  params <- c(Q_trend  = exp(pars[1]), # process error for level component
              Q_season = exp(pars[2]), # process error for seasonal component
              AR1      = KFAS::artransform(pars[3:4])[1], # the AR(1) coefficient
              AR2      = KFAS::artransform(pars[3:4])[2], # the AR(2) coefficient
              Q_ar     = exp(pars[5]), # process error for AR
              H        = exp(pars[6]) # observed error
              )
  return(params)
}


#' Extract coefficients of exogenous variables with confidence intervals
#'
#' Extracts estimated regression coefficients for exogenous variable(s)
#' included in a \code{ThermoSSM} model, together with confidence intervals
#' based on Kalman smoothing results.
#'
#' If the fitted model does not include exogenous variables,
#' the function returns \code{NULL}.
#'
#' @param res An object of class \code{"ThermoSSM"} returned by \code{lgssm()}.
#'
#' @param level Confidence level for the interval estimation.
#'   Must be a numeric value between 0 and 1.
#'   The default is \code{0.95}, corresponding to a 95\% confidence interval.
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
#' res <- lgssm(temp_ts, exogenous = c("kuroshio")
#' extract_exo_coef_ci(res)
#' }
#'
#' @importFrom utils head
#' @export
extract_exo_coef_ci <- function(res, level = 0.95) {
  
  if (!inherits(res, "ThermoSSM")) {
    stop("Input must be a 'ThermoSSM' object.", call. = FALSE)
  }
  
  if (!is.numeric(level) || level <= 0 || level >= 1) {
    stop("'level' must be a numeric value between 0 and 1.", call. = FALSE)
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
  ci_obj <- confint(kfs, level = level)[seq_len(n_exo)]
  ci_mat <- do.call(rbind, lapply(ci_obj, head, n = 1))
  
  data.frame(
    Variable    = exo_vars,
    Coefficient = as.numeric(beta_hat),
    lwr         = ci_mat[, "lwr"],
    upr         = ci_mat[, "upr"],
    row.names   = NULL
  )
}



#' Assign a column label to a univariate \code{ts} object
#'
#' This function assigns a column name to a univariate time series object
#' of class \code{ts}. It is useful when downstream functions require
#' a labeled series.
#'
#' @param ts_in A univariate time series object of class \code{ts}.
#' @param label A character string specifying the column name to assign
#'   (default: \code{"var"}).
#'
#' @return A \code{ts} object with a single column labeled by \code{label}.
#'
#' @examples
#' ts_in <- ts(
#'   rnorm(12 * 30, mean = 10),
#'   start = c(1981, 1),
#'   frequency = 12
#' )
#'
#' ts_labeled <- label_ts_mono(ts_in, label = "var")
#'
#' @export
label_ts_mono <- function(ts_in, label = "var") {

  if (!inherits(ts_in, "ts")) {
    stop("Input must be a ts object.", call. = FALSE)
  }

  # ts がベクトルでも matrix でも安全に処理
  if (!is.null(dim(ts_in)) && ncol(ts_in) != 1) {
    warning("A univariate ts object is expected.", call. = FALSE)
  }

  ts_labeled <- ts(
    matrix(as.numeric(ts_in), ncol = 1,
           dimnames = list(NULL, label)),
    start = start(ts_in),
    frequency = frequency(ts_in)
  )

  return(ts_labeled)
}
