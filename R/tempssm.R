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
    y <- .tempssm_check_temp_ts(temp_data)
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

        exo_name <- NULL
        
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
       
        exo_data_checked <- .tempssm_check_exo_ts(temp_data = temp_data,
                                                           exo_data = exo_data)
        
        exo_name <- colnames(exo_data_checked)
        exo_mat <- as.matrix(exo_data_checked)
        
        ## ---- Model definition -----------------------------------------------
        build_ssm <- SSModel(
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
        
        ## ---- Parameter update function --------------------------------------
        update_func <- function(pars, model) {
          return(
            SSModel(
              H = exp(pars[H_idx]),
              y ~ exo_mat + 
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
        
        exo_name <- NULL
        
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
        
        exo_data_checked <- .tempssm_check_exo_ts(temp_data = temp_data,
                                                           exo_data = exo_data)
        
        exo_name <- colnames(exo_data_checked)
        exo_mat <- as.matrix(exo_data_checked)
        
        
        ## ---- Model definition -----------------------------------------------
        build_ssm <- SSModel(
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
        
        
        ## ---- Parameter update function --------------------------------------
        update_func <- function(pars, model) {
          return(
            SSModel(
              H = exp(pars[H_idx]),
              y ~ exo_mat + 
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


    ## -- Adjust names of parameters ---------------
    state_names <- character(ncol(kfs$alphahat))
    
    idx <- 1
    
    # exogenous
    if (!is.null(exo_name)) {
      state_names[idx:(idx + length(exo_name) - 1)] <- exo_name
      idx <- idx + length(exo_name)
    }
    
    # level
    state_names[idx] <- "level"; idx <- idx + 1
    
    # slope
    state_names[idx] <- "slope"; idx <- idx + 1
    
    # seasonal
    if (use_season) {
      n_season <- freq - 1
      state_names[idx:(idx + n_season - 1)] <-
        paste0("sea_dummy", seq_len(n_season))
      idx <- idx + n_season
    }
    
    # AR
    if (ar_order > 0) {
      state_names[idx:(idx + ar_order - 1)] <-
        paste0("arima", seq_len(ar_order))
    }
    ## -- Adjust names of parameters ---------------
    
        
    ## ---- Output ---------------------------------------------------------
    out <- list(
      model = fit2$model,
      fit   = fit2,
      kfs   = kfs,
      data_temp  = temp_data,
      data_exogenous = exo_data,
      ar_order = ar_order,
      use_season = use_season,
      call  = match.call(),
      converged = fit2$optim.out$convergence == 0,
      state_map = list(
        exogenous = exo_name,
        all       = state_names
        )
    )

    colnames(out$kfs$alphahat) <- state_names
    class(out) <- "tempssm"
    
    return(out)
    
  },
  error = function(e) {
    message("Warning(s): tempssm model did not converge", e$message)
    
    out <- list(
      model = NULL,
      fit = NULL,
      kfs = NULL,
      data_temp = temp_data,
      data_exogenous = exo_data,
      ar_order = ar_order,
      use_season = use_season,
      call = match.call(),
      converged = FALSE,
      state_map = list(
        exogenous = exo_name,
        all       = state_names
      )
    )
    
    colnames(out$kfs$alphahat) <- state_names
    class(out) <- "tempssm"
    return(out)
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

  temp_data_checked <- .tempssm_check_temp_ts(temp_data,
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
