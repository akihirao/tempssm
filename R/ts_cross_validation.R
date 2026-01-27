#' Generate rolling-origin train/test splits for time series cross-validation
#'
#' This function splits a monthly time series (\code{frequency = 12})
#' into training and test sets using a rolling-origin (also known as
#' walk-forward) cross-validation scheme.
#'
#' @param temp_data Monthly temperature time series of class \code{ts}.
#'   The time series must have a frequency of 12.
#' @param exo_data Monthly temperature time series of class \code{ts}.
#'   The time series must have a frequency of 12.
#' @param initial Initial length of the training set (number of observations;
#'   e.g., 60). Default is 60.
#' @param horizon Forecast horizon for the test set (number of observations;
#'   e.g., 12). Default is 12.
#' @param step Step size between successive folds (number of observations;
#'   e.g., 1 or 12). Default is 12.
#' @param fixed_window Logical; if \code{TRUE}, a fixed-length training window
#'   of size \code{initial} is used. If \code{FALSE}, an expanding window
#'   is used. Default is \code{FALSE}.
#' @param allow_partial Logical; if \code{TRUE}, the final fold is included
#'   even when the remaining test period is shorter than \code{horizon}.
#'   Default is \code{FALSE}.
#'
#' @return
#' A list of folds. Each element is a list containing:
#' \describe{
#'   \item{fold}{Fold index.}
#'   \item{train_ts}{Training time series (\code{ts}).}
#'   \item{test_ts}{Test time series (\code{ts}).}
#'   \item{train_idx}{Index range of the training set (position-based).}
#'   \item{test_idx}{Index range of the test set (position-based).}
#'   \item{train_range}{Time range of the training set
#'     (\code{c(start_time, end_time)}).}
#'   \item{test_range}{Time range of the test set
#'     (\code{c(start_time, end_time)}).}
#' }
#'
#' @examples
#' # Example: Air temperature (monthly data, 1932Jul-2025Dec)
#' data(fuji_temp)
#' folds <- ts_cv_folds(
#'   fuji_temp,
#'   initial = 60,
#'   horizon = 12,
#'   step = 12,
#'   fixed_window = FALSE
#' )
#' length(folds)        # number of folds
#' folds[[1]]$train_ts  # training data of the first fold
#' folds[[1]]$test_ts   # test data of the first fold
#'
#' @importFrom stats start end frequency time window
#' @export
ts_cv_folds <- function(temp_data,
                        exo_data = NULL,
                        initial = 60,
                        horizon = 12,
                        step = 12,
                        fixed_window = FALSE,
                        allow_partial = FALSE) {
  
  if (!inherits(temp_data, "ts")) {
    stop("The object of temp_data must be a 'ts' object.", call. = FALSE)
  }
  
  if (frequency(temp_data) != 12) {
    stop("The object of temp_data must be a monthly time series (frequency = 12).", call. = FALSE)
  }
  
  if (initial < 1 || horizon < 1 || step < 1) {
    stop("initial, horizon, and step must be positive integers.", call. = FALSE)
  }
  
  n <- length(temp_data)
  if (initial >= n) {
    stop("initial must be smaller than the length of the time series.", call. = FALSE)
  }
  
  # Utility to slice a ts object using position-based indices
  ts_slice <- function(x, i_start, i_end) {
    if(!is.null(x)){
      sliced_ts <- window(x,
             start = stats::time(temp_data)[i_start],
             end   = stats::time(temp_data)[i_end]
      )
    }else{
      sliced_ts <- NULL
    }
    return(sliced_ts)  
  }
    
  folds <- list()
  k <- 1
  train_end <- initial
  
  while (TRUE) {
    
    # Training window
    train_start <- if (fixed_window) {
      max(1, train_end - initial + 1)
    } else {
      1
    }
    
    # Test window
    test_start <- train_end + 1
    test_end   <- train_end + horizon
    
    if (test_start > n) break
    
    if (test_end > n) {
      if (!allow_partial) break
      test_end <- n
    }
    
    if(is.null(exo_data)){
      exo_train_ts <- NULL
      exo_test_ts <- NULL
    }else{
      exo_train_ts <- ts_slice(exo_data, train_start, train_end)
      exo_test_ts <- ts_slice(exo_data, test_start, test_end)
    }
    
    folds[[k]] <- list(
      fold = k,
      train_ts = ts_slice(temp_data, train_start, train_end),
      test_ts  = ts_slice(temp_data, test_start, test_end),
      exo_train_ts = exo_train_ts,
      exo_test_ts = exo_test_ts,
      train_idx = c(train_start, train_end),
      test_idx  = c(test_start, test_end),
      train_range = c(stats::time(temp_data)[train_start],
                      stats::time(temp_data)[train_end]),
      test_range  = c(stats::time(temp_data)[test_start],
                      stats::time(temp_data)[test_end])
    )
    
    train_end <- train_end + step
    if (train_end >= n) break
    k <- k + 1
  }
  
  return(folds)
}




#' Prediction based on fitted model for time series cross-validation
#'
#' This function predict with using test data for in ts cross-validation.
#'
#' @param folds A \code{list} of folds object returned by \code{ts_cv_folds()}.
#' @param fold_ids A \code{vector} of sequential IDs of folds for ts cross-validation.
#' 
#' @return A \code{list} object of ts cross-validation output
#' 
#' @export
tsCV <- function(folds,fold_ids = 1){
  
  exo_judge <- unlist(lapply(folds, `[[`, "exo_train_ts"))
  ts_CV_list <- list()
  
  num_fold <- length(folds)

  for(i in fold_ids){
      
    target_fold <- folds[[i]]
    temp_train_ts <- target_fold$train_ts
    temp_train_mts <- ts(
      as.matrix(temp_train_ts),
      start = start(temp_train_ts),
      frequency = 12
    )
    colnames(temp_train_mts) <- "Temp" 
    
    temp_test_ts <- target_fold$test_ts
    temp_test_mts <- ts(
      as.matrix(temp_test_ts),
      start = start(temp_test_ts),
      frequency = 12
    )
    colnames(temp_train_mts) <- "Temp" 

    if(is.null(exo_judge)){ # if data-set excludes exogenous variables
      
      res_train <- lgssm(temp_data = temp_train_mts,
                        exo_data = NULL)

      convergence <- res_train$fit$optim.out$convergence
      train_model <- res_train$model
      train_pars <- res_train$fit$optim.out$par
 
      #browser()
      
      predict_temp_ts <- predict(
        train_model,
        newdata = SSModel(
          H = exp(train_pars[6]),
          rep(NA, nrow(temp_test_mts)) ~ #nrow(ts_object) 
            SSMtrend(degree = 2,
                     Q = c(list(0), list(exp(train_pars[1])))) + 
            SSMseasonal(sea.type = "dummy",
                        period = 12,
                        Q = exp(train_pars[2])) +
            SSMarima(ar = artransform(train_pars[3:4]),
                     d = 0,
                     Q = exp(train_pars[5])
            ), data = temp_test_mts
          )
      )
    }else{ # if data-set includes exogenous variables
      
      exo_train <- target_fold$exo_train_ts
      exo_test <- target_fold$exo_test_ts
      num_exo <- dim(exo_train)[2]

      temp_exo_test <- cbind(temp_test_mts,
                             exo_test)
      colnames(temp_exo_test) <- c("Temp",colnames(exo_test))
      
      res_train<- lgssm(temp_data = temp_train_mts,
                        exo_data = exo_train)
      
      convergence <- res_train$fit$optim.out$convergence
      train_model <- res_train$model
      train_pars <- res_train$fit$optim.out$par
      
      
      exogenous_mat <- as.matrix(exo_test)
      
      #browser()
      
      predict_temp_ts <- predict(
        train_model,
        newdata = SSModel(
          H = exp(train_pars[6]),
          rep(NA, nrow(temp_exo_test)) ~ exogenous_mat + 
            SSMtrend(degree = 2,
                     Q = c(list(0), list(exp(train_pars[1])))) + 
            SSMseasonal(sea.type = "dummy",
                        period = 12,
                        Q = exp(train_pars[2])) +
            SSMarima(ar = artransform(train_pars[3:4]),
                     d = 0,
                     Q = exp(train_pars[5])
            ), data = temp_exo_test
        )
      )
    }
    
    CV_output <- forecast::accuracy(predict_temp_ts, temp_test_ts) 
  
    ts_CV_list[[i]] <- list(folds[[i]],
                            convergence=convergence,
                            CV = CV_output)
  
  }
  
  return(ts_CV_list)
}
