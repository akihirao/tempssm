#' Generate rolling-origin train/test splits for time series cross-validation
#'
#' This function splits a monthly time series (\code{frequency = 12})
#' into training and test sets using a rolling-origin (also known as
#' walk-forward) cross-validation scheme.
#'
#' @param temp_data A temperature time series of class \code{ts}.
#'   The series can have any arbitrary frequency of 2 or higher.
#'   For example, a frequency of 12 represents a monthly \code{ts} object.
#' @param exo_data Exogenous time series variable(s) of class \code{ts}.
#'   The series may have any arbitrary frequency of 2 or higher.
#'   But it must be the same as that of \code{temp_data}.
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
#' frequency(fuji_temp) # 12: monthly ts object
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
  
  freq = frequency(temp_data)
  if (freq == 1) {
    stop("The procedure requires a ts object with frequency > 1.",
         call. = FALSE)
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





#' Calculate scale coefficient Q for MASE
#'
#' This function computes the scale coefficient \eqn{Q} used in the
#' Mean Absolute Scaled Error (MASE). The scale is based on in-sample
#' naive or seasonal naive forecasts.
#'
#' For the naive method, \eqn{Q} is defined as the mean absolute
#' first difference:
#' \deqn{Q = \frac{1}{n-1} \sum_{t=2}^{n} |y_t - y_{t-1}|}
#'
#' For the seasonal method, \eqn{Q} is defined as the mean absolute
#' seasonal difference with period \eqn{m = frequency(train_ts)}:
#' \deqn{Q = \frac{1}{n-m} \sum_{t=m+1}^{n} |y_t - y_{t-m}|}
#'
#' @param train_ts A \code{ts} object containing the training time series.
#' @param method Character string specifying the scaling method.
#'   Either \code{"naive"} or \code{"seasonal"}.
#'
#' @return A numeric scalar representing the scale coefficient \eqn{Q}.
#'
#' @references
#' Hyndman, R. J., & Koehler, A. B. (2006).
#' Another look at measures of forecast accuracy.
#' \emph{International Journal of Forecasting}, 22(4), 679–688.
#'
#' @examples
#' ts_data <- ts(rnorm(120), frequency = 12)
#' scale_Q(ts_data, method = "naive")
#' scale_Q(ts_data, method = "seasonal")
#'
#' @export
scale_Q <- function(train_ts, method = c("naive", "seasonal")) {

  method <- match.arg(method)

  if (!inherits(train_ts, "ts")) {
    stop("train_ts must be a 'ts' object.", call. = FALSE)
  }

  y <- as.numeric(train_ts)
  n <- length(y)

  if (method == "naive") {

    if (n < 2) {
      stop("At least two observations are required for naive scaling.",
           call. = FALSE)
    }

    Q <- mean(abs(diff(y)), na.rm = TRUE)

  } else {  # seasonal

    m <- frequency(train_ts)

    if (m <= 1) {
      stop("Seasonal scaling requires a ts object with frequency > 1.",
           call. = FALSE)
    }

    if (n <= m) {
      stop("Time series is too short for seasonal scaling.",
           call. = FALSE)
    }

    Q <- mean(abs(y[(m + 1):n] - y[1:(n - m)]), na.rm = TRUE)
  }

  return(Q)
}



#' Execute cross-validation of each of rolling-origin train/test data
#'
#' This function predict with using test data for in ts cross-validation.
#' @param folds A \code{list} of folds object returned by \code{ts_cv_folds()}.
#'
#' @param idx A \code{vector} of sequential IDs of folds for ts cross-validation.
#'
#' @param ar_order_given Integer specifying the order of the autoregressive (AR)
#' component in the error structure (e.g., 2 for AR(2), 3 for AR(3))
#' for applying \code{tempssm()}.
#' 
#' @return A \code{list} object of ts cross-validation output
#' 
#' @export
exec_tsCV <- function(
  folds,
  idx=1,
  ar_order_given){
  
  if(is.na(ar_order_given)){
    stop(paste("Number of order of autoregressive component should be given."))
  }

  ar_idx <- 3:(2 + ar_order_given)
  var_idx <- 3 + ar_order_given
  H_idx <- 4 + ar_order_given

  exo_judge <- unlist(lapply(folds, `[[`, "exo_train_ts"))
  
  target_fold <- folds[[idx]]
  temp_train_ts <- target_fold$train_ts
  freq = frequency(temp_train_ts)
  
  temp_train_mts <- ts(
    as.matrix(temp_train_ts),
    start = start(temp_train_ts),
    frequency = freq
  )
  colnames(temp_train_mts) <- "Temp" 
  
  temp_test_ts <- target_fold$test_ts
  temp_test_mts <- ts(
    as.matrix(temp_test_ts),
    start = start(temp_test_ts),
    frequency = freq
  )
  colnames(temp_train_mts) <- "Temp" 
  
  if(is.null(exo_judge)){ # if data-set excludes exogenous variables
    
    res_train <- tempssm(temp_data = temp_train_mts,
                       exo_data = NULL,
                       ar_order = ar_order_given)
    
    train_model <- res_train$model
    train_pars <- res_train$fit$optim.out$par
    
    predict_temp_ts <- stats::predict(
      train_model,
      newdata = SSModel(
        H = exp(train_pars[H_idx]),
        rep(NA, nrow(temp_test_mts)) ~ #nrow(ts_object) 
          SSMtrend(degree = 2,
                   Q = c(list(0), list(exp(train_pars[1])))) + 
          SSMseasonal(sea.type = "dummy",
                      period = freq,
                      Q = exp(train_pars[2])) +
          SSMarima(ar = artransform(train_pars[ar_idx]),
                   d = 0,
                   Q = exp(train_pars[var_idx])
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
    
    res_train<- tempssm(temp_data = temp_train_mts,
                      exo_data = exo_train,
                      ar_order = ar_order_given)
    
    train_model <- res_train$model
    train_pars <- res_train$fit$optim.out$par
    
    exogenous_mat <- as.matrix(exo_test)
    
    predict_temp_ts <- stats::predict(
      train_model,
      newdata = SSModel(
        H = exp(train_pars[H_idx]),
        rep(NA, nrow(temp_exo_test)) ~ exogenous_mat + 
          SSMtrend(degree = 2,
                   Q = c(list(0), list(exp(train_pars[1])))) + 
          SSMseasonal(sea.type = "dummy",
                      period = freq,
                      Q = exp(train_pars[2])) +
          SSMarima(ar = artransform(train_pars[ar_idx]),
                   d = 0,
                   Q = exp(train_pars[var_idx])
          ), data = temp_exo_test
      )
    )
  }
  
  convergence_status <- res_train$fit$optim.out$convergence
  if(convergence_status==0){
    convergence <- TRUE
  }else{
    convergence <- FALSE
  }
  
  CV_output <- forecast::accuracy(predict_temp_ts, temp_test_ts) 
  
  Q_naive <- scale_Q(temp_train_ts, method="naive")
  Q_seasonal <- scale_Q(temp_train_ts, method="seasonal")
  
  MASE_naive <- CV_output[,"MAE"]/Q_naive
  
  MASE_snaive <- CV_output[,"MAE"]/Q_seasonal
  
  MASE_matrix <- matrix(c(MASE_naive,MASE_snaive), nrow=1)
  colnames(MASE_matrix) <- c("MASE_naive","MASE_snaive")
  
  CV_output <- cbind(CV_output,MASE_matrix)
  
  CV_output_info <- list(folds[[idx]],
                         convergence=convergence,
                         CV = CV_output)
  
  return(CV_output_info)
}



#' Prediction based on fitted model for time series cross-validation
#'
#' This function performs prediction using test data in time series cross-validation.
#'
#' @importFrom tempssm exec_tsCV
#'
#' @param folds A \code{list} of fold objects returned by \code{ts_cv_folds()}.
#'
#' @param fold_ids An integer vector specifying which folds to run.
#'   Default is all folds.
#'
#' @param ar_order Integer specifying the order of the autoregressive (AR)
#' component in the error structure (e.g., 2 for AR(2), 3 for AR(3))
#' for applying \code{tempssm()}. Defaults to 2.
#'
#' @return A \code{list} object of ts cross-validation output.
#' @export
rolling_origin_tsCV <- function(folds, fold_ids = NULL,ar_order){
  
  num_fold <- length(folds)
  
  if (is.null(fold_ids)) {
    fold_ids <- seq_len(num_fold)
  }
  
  exec_fun <- tempssm::exec_tsCV  # Importance!
  
  # ---- safety check ----
  if (!is.numeric(fold_ids) || any(fold_ids < 1) || any(fold_ids > num_fold)) {
    stop("fold_ids must be valid fold indices.")
  }
  
  ts_CV_list <- future.apply::future_lapply(
    stats::setNames(fold_ids, paste0("fold", fold_ids)),
    function(i){
      exec_fun(
        folds,
        idx = i,
        ar_order_given = ar_order
        )
    },
    future.seed = TRUE,
    future.packages = c("tempssm","KFAS"),
    future.globals = c("folds","exec_fun") # Note globals
  )
  
  return(ts_CV_list)
}




