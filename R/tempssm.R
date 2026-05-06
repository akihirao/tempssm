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
    if (!is.numeric(inits) && length(inits) != expected_len) {
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

    
    ### ---- Model with no exogenous variables
    if(is.null(exo_data)){ 

      exo_name <- NULL
      exo_mat <- NULL

    }else{ # Model with exogenous variables
       
      exo_data_checked <- .tempssm_check_exo_ts(temp_data = temp_data,
                                                exo_data = exo_data)
      
      exo_name <- colnames(exo_data_checked)
      exo_mat <- as.matrix(exo_data_checked)
    } 

    
    #  ---- Define model -------------------------------------
    build_ssm <- .define_build_model(y = y,
                                     freq = freq,
                                     use_season=use_season,
                                     exo_mat = exo_mat,
                                     ar_order = ar_order)  

    ## ---- Parameter update function -------------------------------------
    update_func_common <- .define_update_func(y = y,
                                    freq = freq,
                                    use_season = use_season,
                                    exo_mat = exo_mat,
                                    ar_order = ar_order,
                                    ar_idx = ar_idx,
                                    var_idx = var_idx,
                                    H_idx = H_idx)
    
    ## ---- Optimization (two-step) ----------------------------------------
    fit1 <- fitSSM(
      build_ssm,
      inits  = inits,
      updatefn = update_func_common,
      method = "Nelder-Mead",
      control = list(maxit = maxit, reltol = reltol)
    )
    
    fit2 <- fitSSM(
      build_ssm,
      inits  = fit1$optim.out$par,
      updatefn = update_func_common,
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
    message("Warning(s): tempssm model did not converge. ", e$message)
    
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





#' Internal helper to construct a KFAS state-space model
#'
#' This internal function builds and returns a \code{KFAS::SSModel} object
#' corresponding to the specified model structure. It conditionally includes
#' trend, seasonal, autoregressive, and optional exogenous components,
#' depending on the input arguments.
#'
#' The function is primarily used within \code{tempssm()} to generate the
#' baseline state-space model prior to parameter estimation. Model parameters
#' are left unspecified (set to \code{NA} or zero) and are later updated
#' during optimization via the update function.
#'
#' Four model configurations are supported:
#' \itemize{
#'   \item Trend + Seasonal + AR (no exogenous variables)
#'   \item Trend + AR (no seasonal component)
#'   \item Trend + Seasonal + AR + Exogenous variables
#'   \item Trend + AR + Exogenous variables (no seasonal component)
#' }
#'
#' @param y A numeric vector or univariate \code{ts} object representing
#'   the observed time series.
#' @param freq Integer indicating the seasonal frequency of \code{y}.
#' @param use_season Logical; whether to include a seasonal component.
#' @param exo_mat Optional numeric matrix of exogenous regressors.
#'   If provided, each column is treated as a separate covariate.
#' @param ar_order Integer specifying the order of the autoregressive component.
#'
#' @return A \code{KFAS::SSModel} object with unspecified variance parameters.
#'
#' @details
#' The trend component is modeled as a second-order polynomial trend.
#' The seasonal component (if included) uses a dummy-variable formulation
#' with a sum-to-zero constraint. The autoregressive component is implemented
#' via \code{SSMarima()} with coefficients initialized to zero.
#'
#' @keywords internal
.define_build_model <- function(y = NULL,
                             freq = freq,
                             use_season, 
                             exo_mat,
                             ar_order = 1) {
  if(use_season && is.null(exo_mat)){
    
    build_ssm <- KFAS::SSModel(
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
    
  }else if(!use_season && is.null(exo_mat)){
    
    build_ssm <- KFAS::SSModel(
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
    
  }else if(use_season && !(is.null(exo_mat))){
    
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
  }else if(!use_season && !(is.null(exo_mat))){
    build_ssm <- KFAS::SSModel(
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
  }
}




#' Internal helper to generate the parameter update function for KFAS
#'
#' This internal function constructs and returns an update function to be used
#' in \code{KFAS::fitSSM()}, which maps an unconstrained parameter vector
#' to a valid state-space model.
#'
#' The returned function takes model parameters on an unconstrained scale
#' and transforms them to enforce constraints such as:
#' \itemize{
#'   \item Variance parameters are exponentiated to ensure positivity
#'   \item Autoregressive coefficients are transformed using
#'         \code{KFAS::artransform()} to ensure stationarity
#' }
#'
#' The structure of the update function depends on whether seasonal and/or
#' exogenous components are included, ensuring consistency with the model
#' defined in \code{.define_build_model()}.
#'
#' @param y A numeric vector or univariate \code{ts} object representing
#'   the observed time series.
#' @param freq Integer indicating the seasonal frequency.
#' @param use_season Logical; whether to include a seasonal component.
#' @param exo_mat Optional matrix of exogenous regressors.
#' @param ar_order Integer specifying the order of the autoregressive component.
#' @param ar_idx Integer vector indicating positions of AR coefficients
#'   in the parameter vector.
#' @param var_idx Integer indicating the position of the AR process variance.
#' @param H_idx Integer indicating the position of the observation variance.
#'
#' @return A function with signature \code{function(pars, model)} suitable
#'   for use in \code{KFAS::fitSSM()}.
#'
#' @details
#' The returned function rebuilds the full state-space model using the
#' transformed parameters. It ensures that all variance parameters are positive
#' and that the autoregressive process satisfies stationarity constraints.
#'
#' This closure-based design allows efficient reuse of model structure during
#' numerical optimization without reconstructing external inputs.
#'
#' @keywords internal
.define_update_func <- function(y = NULL,
                                freq = freq,
                                use_season,
                                exo_mat,
                                ar_order = 1,
                                ar_idx,
                                var_idx,
                                H_idx) {
  
  if(use_season && is.null(exo_mat)){
  
    update_func <- function(pars, model) {
      return(
        KFAS::SSModel(
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
      )
    }

  }else if(!use_season && is.null(exo_mat)){
    
    update_func <- function(pars, model) {
      return(
        KFAS::SSModel(
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
      )
    }

  }else if(use_season && !(is.null(exo_mat))){
    
    update_func <- function(pars, model) {
      return(
        KFAS::SSModel(
          H = exp(pars[H_idx]),
          y ~ exo_mat + 
            SSMtrend(
              degree = 2,
              Q = c(list(0),list(exp(pars[1])))
            ) + 
            SSMseasonal(
              sea.type = "dummy",
              period = freq,
              Q = exp(pars[2])
            ) + 
            SSMarima(
              ar = artransform(pars[ar_idx]),
              d = 0,
              Q = exp(pars[var_idx])
            )
        )
      )
    }
    
  }else if(!use_season && !(is.null(exo_mat))){
    
    update_func <- function(pars, model) {
      return(
        KFAS::SSModel(
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
  }
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
