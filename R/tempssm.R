#' Base function for fitting a linear Gaussian state-space model to temperature time series
#'
#' This function estimates a linear Gaussian state-space model (SSM)
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
                    use_season = TRUE
                    ) {
  tryCatch(
    {
    ## ---- Input checks ---------------------------------------------------
    y <- tempssm::.tempssm_check_temp_ts(temp_data)
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
  

    ## ---- Default maxit value -----------------------------------------  
    if (is.null(maxit)) {
      maxit <- 5000
    }
  
    ## ---- Default reltol value -----------------------------------------  
    if (is.null(reltol)) {
      reltol <- 1e-16
    }
  
    
    # index for parameters
    if(use_season){ 
      ar_idx <- 3:(2 + ar_order)
      var_idx <- 3 + ar_order
      H_idx <- 4 + ar_order
    }else{
      ar_idx <- 2:(1 + ar_order)
      var_idx <- 2 + ar_order
      H_idx <- 3 + ar_order
    }

    
    ## ==== Define Model ====================================================
 
    ### ---- Model with seasonal components--------------------
    if(use_season){

      ### ---- Model with no exogenous variables
      if(is.null(exo_data)){ 

        exogenous_lab <- NULL
        
        #### ---- Build Model -------------------------------
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
      
        #### --- Parameter update function ------------------------
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
      
      }else{ # Model with exogenous variables
       
        exo_data_checked <- tempssm::.tempssm_check_exo_ts(temp_data = temp_data,
                                                           exo_data = exo_data)
        
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

      } # close model with seasonal component & exogenous variables
    
    
    }else{ # model without seasonal component
      
      if(is.null(exo_data)){ # Model with no exogenous variables
        
        exogenous_lab <- NULL
        
        #### ---- Model definition -------------------------------      
        build_ssm <- SSModel(
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
      
        #### ---- Parameter update function ----------------------
        update_func <- function(pars, model) {
          model <- SSModel(
            y ~
              SSMtrend(
                degree = 2,
                Q = c(list(0), list(exp(pars[1])))
              ) +
              SSMarima(
                ar = artransform(pars[ar_idx]),
                d = 0,
                Q = exp(pars[var_idx])
              ),
            H = exp(pars[H_idx])
          )
        }
      
      
      }else{# Model with exogenous variables
        
        exo_data_checked <- tempssm::.tempssm_check_exo_ts(temp_data = temp_data,
                                                           exo_data = exo_data)
        
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
                SSMtrend(
                  degree = 2,
                  Q = c(list(0),list(exp(pars[1])))
                ) +
                SSMarima(
                  ar = artransform(pars[ar_idx]),
                  d = 0,
                  Q = exp(pars[var_idx]))
            )
          )
        }
      }# close case when model with exogenous variable(s)


    }
    ## ==== close Define Model ================================
      

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
      use_season = use_season,
      call  = match.call()
    )

    class(out) <- "tempssm"
    return(out)
  },
  error = function(e) {
    message("Warning(s): NA is returned. ", e$message)
    NA
  }
  )# close tryCatch
}



#' Check ts object of temperature time series for applying \code{tempssm()}
#'
#' @param temp_data A temperature time series of class \code{ts}.
#'   The series can have any arbitrary frequency of 2 or higher.
#'   For example, a frequency of 12 represents a monthly \code{ts} object.
#'
#' @param message Logical; Whether Masseage is showed. 
#' Default to TRUE.
#'
#' @return A univariate \code{ts} object.
#'
#' @export
.tempssm_check_temp_ts <- function(temp_data,
                                   message=TRUE) {
  
  if (!inherits(temp_data, "ts")) {
    stop("The object 'temp_data' must be a 'ts' object.",
         call. = FALSE)
  }
  
  freq = frequency(temp_data)
  if (freq <= 1) {
    stop("The procedure requires a ts object with frequency > 1.",
         call. = FALSE)
  }
  
  if (!is.null(dim(temp_data)) && NCOL(temp_data) != 1) {
    stop("The object 'temp_data' must be univariate.",
         call. = FALSE)
  }
  
  if(message){
    message("The ts object 'temp_data' is univariate with frequency ", freq, ".")
  }
  
  return(temp_data)
}



#' Check ts object of exogenous variable(s) for applying \code{tempssm()}
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
.tempssm_check_exo_ts <- function(temp_data,
                                  exo_data) {

  temp_data_checked <- tempssm::.tempssm_check_temp_ts(temp_data,
                                                       message=FALSE)

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
  
  message(paste0("The ts object 'exo_data' is ", uni_multi, " with frequency ", exo_freq, "."))

  return(exo_data)
}



#' Extract the smoothed level component as a time series
#'
#' @param res An object of class \code{"tempssm"} returned by \code{ssm()}.
#' @param ci Logical; should confidence intervals be returned? (default: FALSE)
#' @param ci_level Confidence level for intervals (default: 0.95).
#'
#' @return
#' A univariate \code{ts} object of the smoothed level component.
#' If \code{ci = TRUE}, a multivariate \code{ts} object with columns
#' \code{level}, \code{lwr}, and \code{upr} is returned.
#'
#' @examples
#' \dontrun{
#' data(niigata_sst)
#' res <- tempssm(niigata_sst)
#' level_ts <- extract_level_ts(res)
#' }
#' @export
extract_level_ts = function(res, ci = FALSE, ci_level = 0.95) {

  if (!inherits(res, "tempssm")) {
    stop("Input must be a tempssm object.", call. = FALSE)
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
#' @param res An object of class \code{"tempssm"} returned by \code{ssm()}.
#' @param ci Logical; should confidence intervals be returned? (default: FALSE)
#' @param ci_level Confidence level for intervals (default: 0.95).
#'
#' @return
#' A univariate \code{ts} object of the smoothed drift component (Celsius per year).
#' If \code{ci = TRUE}, a multivariate \code{ts} object with columns
#' \code{level}, \code{lwr}, and \code{upr} is returned.
#'
#' @export
#' @examples
#' \dontrun{
#' data(niigata_sst)
#' res <- tempssm(niigata_sst)
#' drift_ts <- extract_drift_ts(res)
#' }
extract_drift_ts <- function(res, ci = FALSE, ci_level = 0.95) {

  if (!inherits(res, "tempssm")) {
    stop("Input must be a tempssm object.", call. = FALSE)
  }

  if (is.null(res$kfs$alphahat) || !"slope" %in% colnames(res$kfs$alphahat)) {
    stop("Drift (slope) component not found in the smoothing results.", call. = FALSE)
  }

  # seasonal frequency of ts object
  freq <- frequency(res$data_temp)

  # scaling factor (e.g., 12 for monthly data → per year)
  scale <- freq

  drift_ts <- ts(
    res$kfs$alphahat[, "slope"] * scale,
    start <- start(res$data_temp),
    frequency = freq
    )

  if(ci){

    ci_obj <- confint(res$kfs, level = ci_level)
  
    if (!"slope" %in% names(ci_obj)) {
      stop("Slope component not found in confidence intervals.", call. = FALSE)
    }
    
    drift_ts <- ts(
      cbind(
        drift = drift_ts,
        lwr=ci_obj$slope[,"lwr"] * scale,
        upr=ci_obj$slope[,"upr"] * scale
        ),
      start = start(res$data_temp),
      frequency = freq
      )
  }
  return(drift_ts)
}



#' Extract the smoothed seasonal component as a time series
#'
#' @param res An object of class \code{"tempssm"} returned by \code{ssm()}.
#' @param ci Logical; should confidence intervals be returned? (default: FALSE)
#' @param ci_level Confidence level for intervals (default: 0.95).
#'
#' @return
#' A univariate \code{ts} object of the smoothed level component.
#' If \code{ci = TRUE}, a multivariate \code{ts} object with columns
#' \code{level}, \code{lwr}, and \code{upr} is returned.
#'
#' @export
#' @examples
#' \dontrun{
#' data(niigata_sst)
#' res = tempssm(niigata_sst)
#' season_ts = extract_season_ts(res)
#' }
extract_season_ts <- function(res, ci = FALSE, ci_level = 0.95) {

  use_season <- res$use_season

  if (!inherits(res, "tempssm")) {
    stop("Input must be a tempssm object.", call. = FALSE)
  }

  if (!use_season) {
    stop("Seasonal component is not included in this model.", call. = FALSE)
  }

  if (is.null(res$kfs$alphahat) || !"sea_dummy1" %in% colnames(res$kfs$alphahat)) {
    stop("Seasonal component not found in the smoothing results.", call. = FALSE)
  }

  season_ts <- res$kfs$alphahat[, "sea_dummy1"]
  
  if(ci){
    ci_obj <- confint(res$kfs, level = ci_level)

    season_ts <- cbind(
      season = season_ts,
      lwr=ci_obj$ea_dummy1[,"lwr"],
      upr=ci_obj$ea_dummy1[,"upr"]
      )
  }
  return(season_ts)
}




#' Extract the Akaike Information Criterion (AIC)
#'
#' This function computes the Akaike Information Criterion (AIC)
#' for a fitted \code{tempssm} model.
#'
#' @param res An object of class \code{"tempssm"} returned by \code{ssm()}.
#'
#' @return A numeric value representing the AIC of the fitted model.
#'
#' @examples
#' \dontrun{
#' data(niigata_sst)
#' res <- tempssm(niigata_sst)
#' model_AIC <- extract_AIC(res)
#' }
#' @export
extract_AIC <- function(res) {

  if (!inherits(res, "tempssm")) {
    stop("Input must be a tempssm object.", call. = FALSE)
  }

  k <- length(res$fit$optim.out$par)

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
#' @param res An object of class \code{"tempssm"} returned by \code{ssm()}.
#'
#' @return A \code{list} object of the estimated parameters.
#'
#' @export
extract_param <- function(res){

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
#' data(niigata_sst)
#' data(pdo)
#  niigata_sst_pdo <- ts.intersect(niigata_sst,pdo)
#' colnames(niigata_sst_pdo) <- c("Temp", "PDO")
#' niigata_sst_common <- niigata_sst_pdo[,"Temp]
#' pdo_common <- niigata_sst_pdo[,"PDO]
#' pdo_common <- label_ts_mono(nao_common, label = "PDO")
#' res <- ssm(temp_data = niigata_sst_common,exo_data = pdo_common)
#' extract_exo_coef_ci(res)
#' }
#'
#' @importFrom utils head
#' @export
extract_exo_coef_ci = function(res, level = 0.95) {
  
  if (!inherits(res, "tempssm")) {
    stop("Input must be a 'tempssm' object.", call. = FALSE)
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

  # Safely handle both vector and matrix types of ts object
  if (!is.null(dim(ts_in)) && ncol(ts_in) != 1) {
    warning("A univariate ts object is expected.", call. = FALSE)
  }

  x <- as.numeric(ts_in)

  ts_labeled <- ts(
    matrix(
      x,
      ncol = 1,
      dimnames = list(NULL, label)
    ),
    start = start(ts_in),
    frequency = frequency(ts_in)
  )

  return(ts_labeled)
}



#' Trim and align temperature and exogenous time series over their shared period
#'
#' @description
#' `trim_ts_overlap()` aligns a temperature time series and one or more exogenous
#' time series by trimming them to their shared (overlapping) time period.
#' The function returns the trimmed series as a named list of `ts` objects,
#' with consistent variable labels applied for downstream modeling in
#' \code{tempssm()}.
#'
#' Univariate `ts` objects are labeled using \code{label_ts_mono()} to ensure a
#' consistent handling of variable names within the \pkg{tempssm} framework.
#'
#' @param temp_ts
#' A univariate \code{ts} object representing the temperature time series.
#'
#' @param exo_ts
#' A multivariate or univariate \code{ts} object of exogenous variables.
#' Each column represents a distinct exogenous covariate.
#'
#' @param temp_name
#' Character string giving the variable name for the temperature series.
#' This name is applied using \code{label_ts_mono()}.
#'
#' @param exo_name
#' Optional character vector giving variable names for the exogenous variables.
#' If \code{NULL}, default names of the form \code{var1}, \code{var2}, \dots
#' are assigned with a warning.
#'
#' @details
#' The shared time window is determined using \code{ts.intersect()}, and both
#' temperature and exogenous series are trimmed accordingly.
#'
#' For univariate series, variable names are assigned via
#' \code{label_ts_mono()} rather than \code{colnames()} in order to maintain
#' internal consistency required by \code{tempssm()}.
#'
#' @return
#' A named list with the following elements:
#' \describe{
#'   \item{temperature}{A univariate \code{ts} object of the trimmed temperature series.}
#'   \item{exogenous}{A \code{ts} object of the trimmed exogenous variables.}
#' }
#'
#' @seealso
#' [ts.intersect()], [label_ts_mono()], \code{\link{tempssm}}
#'
#' @examples
#' temp_ts <- ts(rnorm(100), start = c(2000, 1), frequency = 12)
#' exo_ts  <- ts(matrix(rnorm(200), ncol = 2),
#'               start = c(2001, 1), frequency = 12)
#'
#' trim_ts_overlap(
#'   temp_ts,
#'   exo_ts,
#'   temp_name = "mean_temp",
#'   exo_name  = c("precip", "solar")
#' )
#'
#' @importFrom stats ts.intersect
#' @export
trim_ts_overlap <- function(
    temp_ts,
    exo_ts,
    temp_name = "temp",
    exo_name = NULL
) {
  
  num_exo_variable <- NCOL(exo_ts)
  
  if (is.null(exo_name)) {
    warning(
      "`exo_name` is not supplied. Using default names: var1, var2, ..."
    )
    exo_name <- paste0("var", seq_len(num_exo_variable))
  } else {
    if (length(exo_name) != num_exo_variable) {
      stop("Length of `exo_name` must equal number of exogenous variables")
    }
  }
  
  ## overlap
  temp_exo_ts_overlap <- ts.intersect(temp_ts, exo_ts)
  
  temp_ts_overlap <- temp_exo_ts_overlap[, 1, drop = FALSE]
  exo_ts_overlap  <- temp_exo_ts_overlap[, -1, drop = FALSE]
  
  ## labels
    temp_ts_overlap <- label_ts_mono(temp_ts_overlap,
                                   label=temp_name)
  
  if(num_exo_variable==1){
    exo_ts_overlap <- label_ts_mono(exo_ts_overlap,
                                    label=exo_name)
  }else{
    colnames(exo_ts_overlap) <- exo_name
  }
  
  
  ## output
  out <- list(
    temperature = temp_ts_overlap,
    exogenous   = exo_ts_overlap
  )
  
  return(out)
}
