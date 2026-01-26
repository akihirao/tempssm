#' Generate rolling-origin train/test splits for time series cross-validation
#'
#' This function splits a monthly time series (\code{frequency = 12})
#' into training and test sets using a rolling-origin (also known as
#' walk-forward) cross-validation scheme.
#'
#' @param x Monthly temperature time series of class \code{ts}.
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
ts_cv_folds <- function(x,
                        initial = 60,
                        horizon = 12,
                        step = 12,
                        fixed_window = FALSE,
                        allow_partial = FALSE) {
  
  if (!inherits(x, "ts")) {
    stop("x must be a 'ts' object.", call. = FALSE)
  }
  
  if (frequency(x) != 12) {
    stop("x must be a monthly time series (frequency = 12).", call. = FALSE)
  }
  
  if (initial < 1 || horizon < 1 || step < 1) {
    stop("initial, horizon, and step must be positive integers.", call. = FALSE)
  }
  
  n <- length(x)
  if (initial >= n) {
    stop("initial must be smaller than the length of the time series.", call. = FALSE)
  }
  
  # Utility to slice a ts object using position-based indices
  ts_slice <- function(x, i_start, i_end) {
    window(x,
           start = stats::time(x)[i_start],
           end   = stats::time(x)[i_end])
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
    
    folds[[k]] <- list(
      fold = k,
      train_ts = ts_slice(x, train_start, train_end),
      test_ts  = ts_slice(x, test_start, test_end),
      train_idx = c(train_start, train_end),
      test_idx  = c(test_start, test_end),
      train_range = c(stats::time(x)[train_start],
                      stats::time(x)[train_end]),
      test_range  = c(stats::time(x)[test_start],
                      stats::time(x)[test_end])
    )
    
    train_end <- train_end + step
    if (train_end >= n) break
    k <- k + 1
  }
  
  return(folds)
}

