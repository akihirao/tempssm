#' Generate rolling train/test splits for time series
#'
#' @description
#' This function generates training and test splits for time series
#' cross-validation using a rolling-origin (walk-forward) scheme.
#' It is intended for model evaluation, including comparison between
#' simple linear models and state space models.
#'
#' @param temp_data
#' A univariate time series of class \code{ts} representing temperature data.
#'
#' @param exo_data
#' Optional exogenous time series variable(s) of class \code{ts}.
#' Must have the same length and frequency as \code{temp_data}.
#'
#' @param initial
#' Initial length of the training set (number of observations).
#' Default is 60.
#'
#' @param horizon
#' Forecast horizon for the test set (number of observations).
#' Default is 12.
#'
#' @param step
#' Step size between successive folds (number of observations).
#' Default is 12.
#'
#' @param fixed_window
#' Logical; if \code{TRUE}, a fixed-length training window of size
#' \code{initial} is used. If \code{FALSE}, an expanding training window
#' is used. Default is \code{FALSE}.
#'
#' @param allow_partial
#' Logical; if \code{TRUE}, include the final fold even if the remaining
#' test period is shorter than \code{horizon}. Default is \code{FALSE}.
#'
#' @return
#' A list of folds. Each element is a list containing:
#' \describe{
#'   \item{fold}{Fold index.}
#'   \item{train_ts}{Training time series.}
#'   \item{test_ts}{Test time series.}
#'   \item{exo_train_ts}{Training exogenous time series (or \code{NULL}).}
#'   \item{exo_test_ts}{Test exogenous time series (or \code{NULL}).}
#'   \item{train_idx}{Index range of the training set (position-based).}
#'   \item{test_idx}{Index range of the test set (position-based).}
#'   \item{train_range}{Time range of the training set.}
#'   \item{test_range}{Time range of the test set.}
#' }
#'
#' @examples
#' \dontrun{
#' data(yamaguchi_sst)  # monthly SST near Yamaguchi Prefecture
#' data(pdo)            # Pacific Decadal Oscillation
#'
#' # synchronize series
#' common_ts <- ts.intersect(yamaguchi_sst, pdo)
#' temp_ts <- common_ts[, "yamaguchi_sst"]
#' pdo_ts  <- common_ts[, "pdo"]
#' pdo_ts <- set_ts_name(pdo_ts, label="PDO") # label
#'
#' folds <- ts_train_test_split(
#'   temp_data = temp_ts,
#'   exo_data  = pdo_ts,
#'   initial   = 60,
#'   horizon   = 12,
#'   step      = 12
#' )
#'
#' # inspect first fold
#' folds[[1]]$train_ts
#' folds[[1]]$test_ts
#' }
#'
#' @importFrom stats window time frequency
#' @export
ts_train_test_split <- function(temp_data,
                                exo_data = NULL,
                                initial = 60,
                                horizon = 12,
                                step = 12,
                                fixed_window = FALSE,
                                allow_partial = FALSE) {

  ## ---- Input checks ---------------------------------------------------
  if (!inherits(temp_data, "ts")) {
    stop("`temp_data` must be a 'ts' object.", call. = FALSE)
  }

  if (!is.null(dim(temp_data)) && ncol(temp_data) != 1) {
    stop("`temp_data` must be a univariate 'ts' object.", call. = FALSE)
  }

  freq <- frequency(temp_data)
  if (freq <= 1) {
    stop("`temp_data` must have frequency > 1.", call. = FALSE)
  }

  if (any(c(initial, horizon, step) < 1)) {
    stop("`initial`, `horizon`, and `step` must be positive integers.",
         call. = FALSE)
  }

  n <- length(temp_data)
  if (initial >= n) {
    stop("`initial` must be smaller than the length of the time series.",
         call. = FALSE)
  }

  if (!is.null(exo_data)) {

    n_exo <- NCOL(exo_data)

    if (!inherits(exo_data, "ts")) {
      stop("`exo_data` must be a 'ts' object.", call. = FALSE)
    }
    if (length(exo_data) != n) {
      stop("`exo_data` must have the same length as `temp_data`.",
           call. = FALSE)
    }
    if (frequency(exo_data) != freq) {
      stop("`exo_data` must have the same frequency as `temp_data`.",
           call. = FALSE)
    }
    if(is.null(colnames(exo_data))){
      exo_data <- tempssm::set_ts_name(exo_data,level = paste0("var",seq_len(n_exo)))
    }
  }

  ## ---- Utility: slice ts by index ------------------------------------
  ts_slice <- function(x, i_start, i_end) {
    if (is.null(x)) return(NULL)
    window(
      x,
      start = time(x)[i_start],
      end   = time(x)[i_end]
    )
  }

  ## ---- Rolling split --------------------------------------------------
  folds <- list()
  k <- 1
  train_end <- initial

  while (TRUE) {

    train_start <- if (fixed_window) {
      max(1, train_end - initial + 1)
    } else {
      1
    }

    test_start <- train_end + 1
    test_end   <- train_end + horizon

    if (test_start > n) break

    if (test_end > n) {
      if (!allow_partial) break
      test_end <- n
    }

    folds[[k]] <- list(
      fold = k,
      train_ts = ts_slice(temp_data, train_start, train_end),
      test_ts  = ts_slice(temp_data, test_start, test_end),
      exo_train_ts = ts_slice(exo_data, train_start, train_end),
      exo_test_ts  = ts_slice(exo_data, test_start, test_end),
      train_idx = c(train_start, train_end),
      test_idx  = c(test_start, test_end),
      train_range = time(temp_data)[c(train_start, train_end)],
      test_range  = time(temp_data)[c(test_start, test_end)]
    )

    train_end <- train_end + step
    if (train_end >= n) break
    k <- k + 1
  }

  folds
}



#' Run time series cross-validation for a single fold (temporary implementation)
#'
#' @description
#' Fit a state space model using \code{tempssm()} on a single
#' training/test split and generate forecasts for the test period.
#' This temporary implementation uses explicit \code{SSModel()}
#' construction for prediction to ensure functionality.
#'
#' @param fold A single fold object returned by \code{ts_train_test_split()}.
#' @param ar_order Integer specifying the AR order. Default is 1.
#' @param use_season Logical; whether to include seasonal component.
#'
#' @return A list containing fold id, convergence status, training/test data,
#' predicted values, and fitted model.
#'
#' @export
ts_cv_run_fold <- function(fold,
                           ar_order = 1,
                           use_season = TRUE) {

  ## ---- basic checks ---------------------------------------------------
  if (!is.list(fold) || is.null(fold$train_ts) || is.null(fold$test_ts)) {
    stop("`fold` must be a valid fold object from ts_train_test_split().",
         call. = FALSE)
  }

  y_train <- fold$train_ts
  y_test  <- fold$test_ts
  exo_train <- fold$exo_train_ts
  exo_test  <- fold$exo_test_ts
  freq <- frequency(y_train)

  ## ---- prepare training ts (matrix) ----------------------------------
  y_train_mts <- ts(
    as.matrix(y_train),
    start = start(y_train),
    frequency = freq
  )
  colnames(y_train_mts) <- "Temp"

  y_test_mts <- ts(
    as.matrix(y_test),
    start = start(y_test),
    frequency = freq
  )
  colnames(y_test_mts) <- "Temp"

  ## ---- fit model ------------------------------------------------------
  res <- tryCatch(
    tempssm(
      temp_data = y_train_mts,
      exo_data  = exo_train,
      ar_order  = ar_order,
      use_season = use_season
    ),
    error = function(e) NULL
  )

  train_model <- res$model
  
  if (is.null(res) || !isTRUE(res$converged)) {
    return(list(
      fold      = fold$fold,
      converged = FALSE,
      y_train   = y_train,
      y_test    = y_test,
      y_pred    = NULL,
      model     = res$model
    ))
  }


  ## ---- extract fitted parameters -------------------------------------
  pars <- res$fit$optim.out$par

  ar_idx  <- 3:(2 + ar_order)
  var_idx <- 3 + ar_order
  H_idx   <- 4 + ar_order

  ## ---- prediction -----------------------------------------------------
  #
  #y_pred <- tryCatch({

    if (is.null(exo_train)) {

      y_pred <- stats::predict(
        train_model,
        n.ahead = length(y_test_mts)
      )

    } else {

      
      num_exo <- dim(exo_train)[2]
      
      temp_exo_test <- cbind(y_test_mts,
                             exo_test)
      colnames(temp_exo_test) <- c("Temp",colnames(exo_test))
      
      res_train<- tempssm(temp_data = y_train_mts,
                          exo_data = exo_train,
                          ar_order = ar_order)
      
      train_model <- res_train$model
      train_pars <- res_train$fit$optim.out$par
      
      exo_mat <- as.matrix(exo_test)
      
      y_pred <- stats::predict(
        train_model,
        newdata = SSModel(
          H = exp(train_pars[H_idx]),
          rep(NA, nrow(temp_exo_test)) ~ exo_mat + 
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
      
      #exo_mat <- as.matrix(exo_test)
      #y_test_exo <- cbind(y_test_mts, exo_test)

      #y_pred <- stats::predict(
      #  train_model,
      #  newdata = KFAS::SSModel(
      #    H = exp(pars[H_idx]),
      #    rep(NA, nrow(y_test_exo)) ~ exo_mat +
      #      KFAS::SSMtrend(
      #        degree = 2,
      #        Q = list(0, exp(pars[1]))) +
      #      KFAS::SSMseasonal(
      #          sea.type = "dummy",
      #          period = freq,
      #          Q = exp(pars[2])) +
      #      KFAS::SSMarima(
      #        ar = KFAS::artransform(pars[ar_idx]),
      #        d = 0,
      #        Q = exp(pars[var_idx])
      #      ),
      #    data = y_test_exo
      #  )
      #)
    }

  #}, error = function(e) NULL)

  ## ---- return ---------------------------------------------------------
  list(
    fold      = fold$fold,
    converged = TRUE,
    y_train   = y_train,
    y_test    = y_test,
    y_pred    = y_pred,
    model     = res
  )
}




#' Run time series cross-validation over multiple folds
#'
#' @description
#' Apply \code{ts_cv_run_fold()} to multiple train/test splits and
#' collect the results. This function orchestrates time series
#' cross-validation but does not compute evaluation metrics.
#'
#' @param folds
#' A list of fold objects returned by \code{ts_train_test_split()}.
#'
#' @param fold_ids
#' Integer vector specifying which folds to run.
#' Default is \code{NULL}, which runs all folds.
#'
#' @param ar_order
#' Integer specifying the order of the autoregressive (AR) component
#' used in \code{tempssm()}. Default is 1.
#'
#' @param use_season
#' Logical; if \code{TRUE}, include a seasonal component in the model.
#' Default is \code{TRUE}.
#'
#' @return
#' A list of cross-validation results. Each element corresponds to one
#' fold and contains the output of \code{ts_cv_run_fold()}.
#'
#' @export
ts_cv_run <- function(folds,
                      fold_ids = NULL,
                      ar_order = 1,
                      use_season = TRUE) {

  ## ---- basic checks ---------------------------------------------------
  if (!is.list(folds) || length(folds) == 0) {
    stop("`folds` must be a non-empty list returned by ts_train_test_split().",
         call. = FALSE)
  }

  if (!is.logical(use_season) || length(use_season) != 1) {
    stop("`use_season` must be a single logical value.",
         call. = FALSE)
  }

  n_folds <- length(folds)

  if (is.null(fold_ids)) {
    fold_ids <- seq_len(n_folds)
  }

  if (!is.numeric(fold_ids) ||
      any(fold_ids < 1) ||
      any(fold_ids > n_folds)) {
    stop("`fold_ids` must contain valid fold indices.",
         call. = FALSE)
  }

  ## ---- run CV ---------------------------------------------------------
  results <- vector("list", length(fold_ids))
  names(results) <- paste0("fold", fold_ids)

  for (i in seq_along(fold_ids)) {
    idx <- fold_ids[i]
    results[[i]] <- ts_cv_run_fold(
      fold = folds[[idx]],
      ar_order = ar_order,
      use_season = use_season
    )
  }

  results
}



## ---- Assessment Index -------------------------------------------


#' Compute mean absolute error (MAE)
#'
#' @param y_pred Predicted values.
#' @param y_true Observed values.
#'
#' @return Numeric scalar giving MAE.
#' @export
compute_mae <- function(y_pred, y_true) {

  if (is.null(y_pred) || is.null(y_true)) {
    return(NA_real_)
  }

  if (length(y_pred) != length(y_true)) {
    stop("`y_pred` and `y_true` must have the same length.",
         call. = FALSE)
  }

  mean(abs(as.numeric(y_true) - as.numeric(y_pred)), na.rm = TRUE)
}


#' Compute Mean Absolute Scaled Error (MASE)
#'
#' @description
#' Compute Mean Absolute Scaled Error (MASE) using training data
#' for scaling and test data for evaluation.
#'
#' @param y_pred Predicted values for the test period.
#' @param y_true Observed values for the test period.
#' @param y_train Training time series used for scaling.
#' @param method Scaling method: \code{"naive"} or \code{"seasonal"}.
#'
#' @return Numeric scalar giving MASE, or \code{NA} if unavailable.
#'
#' @references
#' Hyndman, R. J., & Koehler, A. B. (2006).
#' Another look at measures of forecast accuracy.
#' International Journal of Forecasting, 22(4), 679–688.
#'
#' @export
compute_mase <- function(y_pred,
                         y_true,
                         y_train,
                         method = c("naive", "seasonal")) {

  method <- match.arg(method)

  if (is.null(y_pred) || is.null(y_true)) {
    return(NA_real_)
  }

  if (!inherits(y_train, "ts")) {
    stop("`y_train` must be a 'ts' object.",
         call. = FALSE)
  }

  mae <- compute_mae(y_pred, y_true)
  if (is.na(mae)) {
    return(NA_real_)
  }

  Q <- tryCatch(
    scale_Q(y_train, method = method),
    error = function(e) NA_real_
  )

  if (is.na(Q) || Q == 0) {
    return(NA_real_)
  }

  mae / Q
}



#' Compute forecast accuracy metrics for a CV fold
#'
#' @param cv_result A single result returned by \code{ts_cv_run_fold()}.
#'
#' @return A named list of accuracy metrics.
#' @export
compute_cv_metrics <- function(cv_result) {

  if (!isTRUE(cv_result$converged)) {
    return(list(
      MAE = NA_real_,
      MASE_naive = NA_real_,
      MASE_seasonal = NA_real_
    ))
  }

  list(
    MAE = compute_mae(cv_result$y_pred, cv_result$y_test),
    MASE_naive = compute_mase(
      cv_result$y_pred,
      cv_result$y_test,
      cv_result$y_train,
      method = "naive"
    ),
    MASE_seasonal = compute_mase(
      cv_result$y_pred,
      cv_result$y_test,
      cv_result$y_train,
      method = "seasonal"
    )
  )
}



#' Collect time series cross-validation results into a tibble
#'
#' @description
#' Combine cross-validation results and evaluation metrics into
#' a tidy \code{tibble} with one row per fold.
#'
#' @param cv_results
#' A list of results returned by \code{ts_cv_run()}.
#'
#' @param metrics
#' A list of evaluation metrics returned by
#' \code{compute_cv_metrics()}, corresponding to \code{cv_results}.
#'
#' @return
#' A \code{tibble} where each row represents a cross-validation fold.
#'
#' @importFrom tibble tibble
#' @export
ts_cv_collect <- function(cv_results, metrics) {

  if (!is.list(cv_results) || !is.list(metrics)) {
    stop("`cv_results` and `metrics` must both be lists.", call. = FALSE)
  }

  if (length(cv_results) != length(metrics)) {
    stop("`cv_results` and `metrics` must have the same length.",
         call. = FALSE)
  }

  tibble::tibble(
    fold = as.numeric(lapply(cv_results, function(x) x$fold)),
    converged = as.numeric(lapply(cv_results, function(x) x$converged)),
    MAE = as.numeric(lapply(metrics, function(x) x$MAE)),
    MASE_naive = as.numeric(lapply(metrics, function(x) x$MASE_naive)),
    MASE_seasonal = as.numeric(lapply(metrics, function(x) x$MASE_seasonal))
  )
}




#################################################################
#################################################################

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






