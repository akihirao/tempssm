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
#' data(yamaguchi_sst) # monthly SST near Yamaguchi Prefecture
#' data(pdo) # Pacific Decadal Oscillation
#'
#' # synchronize series
#' common_ts <- ts.intersect(yamaguchi_sst, pdo)
#' temp_ts <- common_ts[, "yamaguchi_sst"]
#' pdo_ts <- common_ts[, "pdo"]
#' pdo_ts <- set_ts_name(pdo_ts, label = "PDO") # label
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

  ## ---- Start message --------------------------------------------------
  .tempssm_cli_inform(
    "Creating rolling train/test data (initial={initial}, horizon={horizon}, step={step})"
  )

  ## ---- Input checks ---------------------------------------------------
  if (!inherits(temp_data, "ts")) {
    cli::cli_abort("`temp_data` must be a {.cls ts} object.")
  }

  if (!is.null(dim(temp_data)) && ncol(temp_data) != 1) {
    cli::cli_abort("`temp_data` must be a univariate {.cls ts} object.")
  }

  freq <- frequency(temp_data)

  if (freq <= 1) {
    cli::cli_abort("`temp_data` must have frequency > 1.")
  }

  if (any(c(initial, horizon, step) < 1)) {
    cli::cli_abort("`initial`, `horizon`, and `step` must be positive integers.")
  }

  n <- length(temp_data)

  if (initial >= n) {
    cli::cli_abort("`initial` must be smaller than the length of the time series.")
  }

  .tempssm_cli_debug("Series length: {n}, frequency: {freq}")

  ## ---- Exogenous checks ----------------------------------------------
  if (!is.null(exo_data)) {

    if (!inherits(exo_data, "ts")) {
      cli::cli_abort("`exo_data` must be a {.cls ts} object.")
    }

    if (length(exo_data) != n) {
      cli::cli_abort("`exo_data` must have the same length as `temp_data`.")
    }

    if (frequency(exo_data) != freq) {
      cli::cli_abort("`exo_data` must have the same frequency as `temp_data`.")
    }

    if (is.null(colnames(exo_data))) {
      cli::cli_warn("No column names in `exo_data`; assigning default names.")
      exo_data <- tempssm::set_ts_name(
        exo_data,
        label = paste0("var", seq_len(NCOL(exo_data)))
      )
    }

    .tempssm_cli_debug("Exogenous variables: {NCOL(exo_data)}")

  } else {
    .tempssm_cli_debug("No exogenous variables provided")
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

    .tempssm_cli_debug(
      "Fold {k}: train[{train_start}:{train_end}], test[{test_start}:{test_end}]"
    )

    folds[[k]] <- list(
      fold = k,
      train_ts = .ts_slice(temp_data, train_start, train_end),
      test_ts  = .ts_slice(temp_data, test_start, test_end),

      exo_train_ts = .ts_slice(exo_data, train_start, train_end),
      exo_test_ts  = .ts_slice(exo_data, test_start, test_end),

      train_idx = c(train_start, train_end),
      test_idx  = c(test_start, test_end),

      train_range = time(temp_data)[c(train_start, train_end)],
      test_range  = time(temp_data)[c(test_start, test_end)]
    )

    train_end <- train_end + step
    if (train_end >= n) break

    k <- k + 1
  }

  ## ---- Completion -----------------------------------------------------
  .tempssm_cli_inform(
    "Created {length(folds)} rolling fold{?s}"
  )

  return(folds)
}


# internal: slice ts object by index range
.ts_slice <- function(x, i_start, i_end) {

  if (is.null(x)) {
    return(NULL)
  }

  if (!inherits(x, "ts")) {
    cli::cli_abort("Input must be a {.cls ts} object.")
  }

  n <- length(x)

  if (i_start < 1 || i_end > n || i_start > i_end) {
    cli::cli_abort("Invalid slice indices for time series.")
  }

  stats::window(
    x,
    start = stats::time(x)[i_start],
    end   = stats::time(x)[i_end]
  )
}


#' Run time series cross-validation for a single fold
#'
#' @description
#' Fits a state space model using \code{tempssm()} on a single
#' training/test split generated by \code{ts_train_test_split()},
#' and produces forecasts for the corresponding test period.
#'
#' This function is primarily used internally by
#' \code{ts_cv_run()} as part of the cross-validation workflow.
#' However, it is also exported for advanced users who wish to
#' inspect or debug model behavior on individual folds.
#'
#' For most use cases, users should call \code{ts_cv_run()}
#' instead of using this function directly.
#'
#' For models without exogenous variables, forecasts are generated using
#' \code{predict(..., n.ahead = h)}. When exogenous variables are
#' present, a new \code{SSModel} object is constructed for the test
#' period using the estimated parameters.
#'
#' \strong{Note:} The current implementation for models with
#' exogenous variables relies on explicit reconstruction of
#' \code{SSModel()}, and may be sensitive to model specification.
#' This part is subject to future improvement.
#'
#' @param fold
#' A single fold object returned by \code{ts_train_test_split()}.
#'
#' @param ar_order
#' Integer specifying the AR order. Default is 1.
#'
#' @param use_season
#' Logical; whether to include a seasonal component. Default is \code{TRUE}.
#'
#' @return
#' A list with the following components:
#' \describe{
#'   \item{fold}{Fold index.}
#'   \item{converged}{Logical indicating whether model fitting succeeded.}
#'   \item{y_train}{Training time series (\code{ts}).}
#'   \item{y_test}{Test time series (\code{ts}).}
#'   \item{y_pred}{Predicted values for the test period (object returned by \code{predict}).}
#'   \item{model}{Fitted model object returned by \code{tempssm()} (or \code{NULL} if failed).}
#' }
#'
#' @details
#' When no exogenous variables are supplied, forecasts are obtained
#' directly via \code{predict()} with \code{n.ahead}. In this case,
#' future exogenous effects are implicitly treated as zero.
#'
#' When exogenous variables are present, forecasts depend on the
#' provided test-period covariates. Users should ensure that
#' \code{exo_test_ts} represents realistic or observed future values.
#'
#' @seealso
#' \code{\link{ts_cv_run}} for running cross-validation over all folds.
#'
#' @examples
#' \dontrun{
#' data(yamaguchi_sst)
#'
#' # create folds (no exogenous variables)
#' folds <- ts_train_test_split(
#'   temp_data = yamaguchi_sst,
#'   initial   = 60,
#'   horizon   = 12,
#'   step      = 12
#' )
#'
# run CV on first fold (useful for debugging / inspection)
#' res <- ts_cv_run_fold(folds[[1]])
#'
#' # inspect predictions
#' res$y_pred
#'
#' # compare observed vs predicted
#' cbind(
#'   observed = as.numeric(res$y_test),
#'   predicted = as.numeric(res$y_pred)
#' )
#' }
#'
#' @export
ts_cv_run_fold <- function(fold,
                           ar_order = 1,
                           use_season = TRUE) {
  ## ---- Start ----------------------------------------------------------
  .tempssm_cli_inform("Running CV for fold {fold$fold}")

  ## ---- basic checks ---------------------------------------------------
  if (!is.list(fold) || is.null(fold$train_ts) || is.null(fold$test_ts)) {
    cli::cli_abort(
      "`fold` must be a valid object from {.fn ts_train_test_split}."
    )
  }

  y_train <- fold$train_ts
  y_test <- fold$test_ts
  exo_train <- fold$exo_train_ts
  exo_test <- fold$exo_test_ts
  freq <- frequency(y_train)

  .tempssm_cli_debug(
    "Fold {fold$fold}: train length={length(y_train)}, test length={length(y_test)}"
  )

  ## ---- prepare training ts -------------------------------------------
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
  .tempssm_cli_inform("Fitting model for fold {fold$fold}")

  res <- tryCatch(
    tempssm(
      temp_data = y_train_mts,
      exo_data = exo_train,
      ar_order = ar_order,
      use_season = use_season
    ),
    error = function(e) {
      cli::cli_warn(
        "Model fitting failed in fold {fold$fold}: {conditionMessage(e)}"
      )
      NULL
    }
  )

  train_model <- if (!is.null(res)) res$model else NULL

  ## ---- check convergence ---------------------------------------------
  if (is.null(res) || !isTRUE(res$converged)) {
    cli::cli_warn(
      "Model did not converge for fold {fold$fold}"
    )

    return(list(
      fold      = fold$fold,
      converged = FALSE,
      y_train   = y_train,
      y_test    = y_test,
      y_pred    = NULL,
      model     = train_model
    ))
  }

  ## ---- prediction -----------------------------------------------------
  .tempssm_cli_inform(
    "Generating predictions for fold {fold$fold}"
  )

  if (is.null(exo_train)) {
    ## ---- no exogenous -------------------------------------------------
    .tempssm_cli_debug("Prediction without exogenous variables")

    y_pred <- stats::predict(
      train_model,
      n.ahead = length(y_test_mts)
    )
  } else {
    ## ---- with exogenous ----------------------------------------------
    .tempssm_cli_debug("Prediction with exogenous variables")

    pars <- res$fit$optim.out$par

    temp_exo_test <- cbind(y_test_mts, exo_test)
    colnames(temp_exo_test) <- c("Temp", colnames(exo_test))

    res_train <- tempssm(
      temp_data = y_train_mts,
      exo_data = exo_train,
      ar_order = ar_order
    )

    train_model <- res_train$model
    train_pars <- res_train$fit$optim.out$par

    exo_mat <- as.matrix(exo_test)

    if (use_season) {
      .tempssm_cli_debug("Using seasonal model for prediction")

      ar_idx <- 3:(2 + ar_order)
      var_idx <- 3 + ar_order
      H_idx <- 4 + ar_order

      y_pred <- stats::predict(
        train_model,
        newdata = SSModel(
          H = exp(train_pars[H_idx]),
          rep(NA, nrow(temp_exo_test)) ~ exo_mat +
            SSMtrend(
              degree = 2,
              Q = c(list(0), list(exp(train_pars[1])))
            ) +
            SSMseasonal(
              sea.type = "dummy",
              period = freq,
              Q = exp(train_pars[2])
            ) +
            SSMarima(
              ar = artransform(train_pars[ar_idx]),
              d = 0,
              Q = exp(train_pars[var_idx])
            ),
          data = temp_exo_test
        )
      )
    } else {
      .tempssm_cli_debug("Using non-seasonal model for prediction")

      ar_idx <- 2:(1 + ar_order)
      var_idx <- 2 + ar_order
      H_idx <- 3 + ar_order

      y_pred <- stats::predict(
        train_model,
        newdata = SSModel(
          H = exp(train_pars[H_idx]),
          rep(NA, nrow(temp_exo_test)) ~ exo_mat +
            SSMtrend(
              degree = 2,
              Q = c(list(0), list(exp(train_pars[1])))
            ) +
            SSMarima(
              ar = artransform(train_pars[ar_idx]),
              d = 0,
              Q = exp(train_pars[var_idx])
            ),
          data = temp_exo_test
        )
      )
    }
  }

  ## ---- completion -----------------------------------------------------
  .tempssm_cli_inform(
    "Completed fold {fold$fold} successfully"
  )

  ## ---- return ---------------------------------------------------------
  list(
    fold      = fold$fold,
    converged = TRUE,
    y_train   = y_train,
    y_test    = y_test,
    y_pred    = y_pred,
    model     = train_model
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
#' @param ar_order
#' Integer specifying the order of the autoregressive (AR) component
#' used in \code{tempssm()}. Default is 1.
#'
#' @param use_season
#' Logical; if \code{TRUE}, include a seasonal component in the model.
#' Default is \code{TRUE}.
#'
#' @param parallel
#' Logical; if \code{TRUE}, folds are evaluated in parallel using the
#' \pkg{future.apply} framework. If \code{FALSE}, folds are processed
#' sequentially. Default is \code{TRUE}.
#'
#' @param workers
#' Integer specifying the number of parallel workers to use when
#' \code{parallel = TRUE}. The default uses all available cores as
#' returned by \code{\link[future]{availableCores}}.
#'
#' @param progress
#' Logical; if \code{TRUE}, a progress bar is displayed during execution
#' using the \pkg{progressr} package. If \code{FALSE}, no progress bar
#' is shown. Default is \code{TRUE}.
#'
#' @return
#' A list of cross-validation results. Each element corresponds to one
#' fold and contains the output of \code{ts_cv_run_fold()}.
#'
#' @examples
#' \dontrun{
#' data(yamaguchi_sst)
#'
#' # 1. create train/test splits
#' folds <- ts_train_test_split(
#'   temp_data = yamaguchi_sst,
#'   initial   = 60,
#'   horizon   = 12,
#'   step      = 12
#' )
#'
#' # 2. run cross-validation across all folds
#' cv_results <- ts_cv_run(
#'   folds,
#'   ar_order = 1,
#'   use_season = TRUE,
#'   parallel = TRUE
#' )
#'
#' # inspect first fold result
#' cv_results[[1]]
#'
#' # predicted vs observed for first fold
#' cbind(
#'   observed = as.numeric(cv_results[[1]]$y_test),
#'   predicted = as.numeric(cv_results[[1]]$y_pred)
#' )
#' }
#'
#' @export
ts_cv_run <- function(
  folds,
  ar_order = 1,
  use_season = TRUE,
  parallel = TRUE,
  workers = future::availableCores(),
  progress = FALSE
) {
  ## ---- input check ----------------------------------------------------
  if (!is.list(folds)) {
    cli::cli_abort("`folds` must be a list of fold objects.")
  }

  n_folds <- length(folds)

  .tempssm_cli_inform(
    "Running cross-validation on {n_folds} fold{?s}"
  )

  ## ---- plan -----------------------------------------------------------
  if (parallel) {
    future::plan(future::multisession, workers = workers)
    .tempssm_cli_debug("Parallel execution with {workers} workers")
  } else {
    future::plan(future::sequential)
    .tempssm_cli_debug("Sequential execution")
  }

  ## ---- runner function -----------------------------------------------
  run_one <- function(f) {
    ts_cv_run_fold(
      fold = f,
      ar_order = ar_order,
      use_season = use_season
    )
  }

  ## ---- execution ------------------------------------------------------

  if (progress && requireNamespace("progressr", quietly = TRUE)) {

    res <- progressr::with_progress({
      p <- progressr::progressor(along = folds)

      future.apply::future_lapply(
        folds,
        function(f) {
          out <- run_one(f)
          p()
          out
        },
        future.seed = TRUE
      )
    })

  } else {

    if (progress) {
      cli::cli_warn(
        "progressr not installed; running without progress bar"
      )
    }

    res <- future.apply::future_lapply(
      folds,
      run_one,
      future.seed = TRUE
    )
  }

  ## ---- completion -----------------------------------------------------
  .tempssm_cli_inform(
    "Completed cross-validation over {n_folds} fold{?s}"
  )

  return(res)
}


## ---- Assessment Index -------------------------------------------


#' Compute mean absolute error (MAE)
#'
#' @param y_pred Predicted values.
#' @param y_true Observed values.
#'
#' @return Numeric scalar giving MAE.
#'
#' @keywords internal
#' @noRd
.compute_mae <- function(y_pred, y_true) {
  ## ---- fast exit ------------------------------------------------------
  if (is.null(y_pred) || is.null(y_true)) {
    return(NA_real_)
  }

  ## ---- input check ----------------------------------------------------
  if (length(y_pred) != length(y_true)) {
    cli::cli_abort(
      "`y_pred` and `y_true` must have the same length."
    )
  }

  .tempssm_cli_debug(
    "Computing MAE for {length(y_pred)} observation{?s}"
  )

  ## ---- compute --------------------------------------------------------
  mae <- mean(
    abs(as.numeric(y_true) - as.numeric(y_pred)),
    na.rm = TRUE
  )

  .tempssm_cli_debug(
    "MAE computed: {round(mae, 4)}"
  )

  return(mae)
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
#' @keywords internal
#' @noRd
.compute_mase <- function(y_pred,
                         y_true,
                         y_train,
                         method = c("naive", "seasonal")) {
  method <- match.arg(method)

  ## ---- fast exit ------------------------------------------------------
  if (is.null(y_pred) || is.null(y_true)) {
    return(NA_real_)
  }

  ## ---- input check ----------------------------------------------------
  if (!inherits(y_train, "ts")) {
    cli::cli_abort(
      "`y_train` must be a {.cls ts} object."
    )
  }

  .tempssm_cli_debug(
    "Computing MASE (method = {method})"
  )

  ## ---- MAE ------------------------------------------------------------
  mae <- .compute_mae(y_pred, y_true)

  if (is.na(mae)) {
    return(NA_real_)
  }

  ## ---- scaling factor -------------------------------------------------
  Q <- tryCatch(
    .scale_Q(y_train, method = method),
    error = function(e) {
      .tempssm_cli_debug(
        "Failed to compute scaling factor Q"
      )
      NA_real_
    }
  )

  if (is.na(Q) || Q == 0) {
    .tempssm_cli_debug("Invalid scaling factor (Q is NA or zero)")
    return(NA_real_)
  }

  ## ---- result ---------------------------------------------------------
  mase <- mae / Q

  .tempssm_cli_debug("MASE computed: {round(mase, 4)}")

  return(mase)
}


#' Compute forecast accuracy metrics for a CV fold
#'
#' @param cv_result A single result returned by \code{ts_cv_run_fold()}.
#'
#' @examples
#' \dontrun{
#' data(yamaguchi_sst)
#'
#' # 1. create train/test splits
#' folds <- ts_train_test_split(
#'   temp_data = yamaguchi_sst,
#'   initial   = 60,
#'   horizon   = 12,
#'   step      = 12
#' )
#'
#' # 2. run cross-validation on a single fold
#' res <- ts_cv_run_fold(folds[[1]])
#'
#' # 3. compute evaluation metrics
#' metrics <- compute_cv_metrics(res)
#'
#' metrics
#' }
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
    MAE = .compute_mae(cv_result$y_pred, cv_result$y_test),
    MASE_naive = .compute_mase(
      cv_result$y_pred,
      cv_result$y_test,
      cv_result$y_train,
      method = "naive"
    ),
    MASE_seasonal = .compute_mase(
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
#'
#' @examples
#' \dontrun{
#' data(yamaguchi_sst)
#'
#' # 1. create train/test splits
#' folds <- ts_train_test_split(
#'   temp_data = yamaguchi_sst,
#'   initial   = 60,
#'   horizon   = 12,
#'   step      = 12
#' )
#'
#' # 2. run cross-validation over all folds
#' cv_results <- ts_cv_run(folds)
#'
#' # 3. compute metrics for each fold
#' metrics <- lapply(cv_results, compute_cv_metrics)
#'
#' # 4. collect results into a tidy tibble
#' cv_summary <- ts_cv_collect(cv_results, metrics)
#'
#' cv_summary
#' }
#'
#' @export
ts_cv_collect <- function(cv_results, metrics) {
  ## ---- input check ----------------------------------------------------
  if (!is.list(cv_results) || !is.list(metrics)) {
    cli::cli_abort(
      "`cv_results` and `.arg metrics` must both be lists."
    )
  }

  if (length(cv_results) != length(metrics)) {
    cli::cli_abort(
      "`cv_results` and `.arg metrics` must have the same length."
    )
  }

  n <- length(cv_results)

  .tempssm_cli_debug(
    "Collecting cross-validation results for {n} fold{?s}"
  )

  ## ---- build tibble ---------------------------------------------------
  out <- tibble::tibble(
    fold = as.numeric(lapply(cv_results, function(x) x$fold)),
    converged = as.numeric(lapply(cv_results, function(x) x$converged)),
    MAE = as.numeric(lapply(metrics, function(x) x$MAE)),
    MASE_naive = as.numeric(lapply(metrics, function(x) x$MASE_naive)),
    MASE_seasonal = as.numeric(lapply(metrics, function(x) x$MASE_seasonal))
  )

  return(out)
}



#' Compute scale coefficient Q for MASE
#'
#' @description
#' Compute the scaling coefficient \eqn{Q} used in the Mean Absolute
#' Scaled Error (MASE). The scaling is based on in-sample naive or
#' seasonal naive differences of the training time series.
#'
#' @details
#' The scale coefficient \eqn{Q} is defined as follows:
#'
#' For the naive method:
#' \deqn{
#' Q = \frac{1}{n-1} \sum_{t=2}^{n} |y_t - y_{t-1}|
#' }
#'
#' For the seasonal method (with period \eqn{m = frequency(train_ts)}):
#' \deqn{
#' Q = \frac{1}{n-m} \sum_{t=m+1}^{n} |y_t - y_{t-m}|
#' }
#'
#' This function is primarily intended for internal use within the
#' \pkg{tempssm} package to support the computation of MASE.
#'
#' @param train_ts A univariate time series object of class \code{ts}.
#' @param method Character string specifying the scaling method.
#'   Must be either \code{"naive"} or \code{"seasonal"}.
#'
#' @return A numeric scalar representing the scale coefficient \eqn{Q}.
#'
#' @references
#' Hyndman, R. J., & Koehler, A. B. (2006).
#' Another look at measures of forecast accuracy.
#' \emph{International Journal of Forecasting}, 22(4), 679–688.
#'
#' @keywords internal
#' @noRd
.scale_Q <- function(train_ts, method = c("naive", "seasonal")) {
  method <- match.arg(method)

  ## ---- input check ----------------------------------------------------
  if (!inherits(train_ts, "ts")) {
    cli::cli_abort(
      "`train_ts` must be a {.cls ts} object."
    )
  }

  y <- as.numeric(train_ts)
  n <- length(y)

  .tempssm_cli_debug(
    "Computing scale Q (method = {method}, n = {n})"
  )

  ## ---- naive scaling --------------------------------------------------
  if (method == "naive") {
    if (n < 2) {
      cli::cli_abort(
        "At least two observations are required for naive scaling."
      )
    }

    Q <- mean(abs(diff(y)), na.rm = TRUE)
  } else { # seasonal

    m <- stats::frequency(train_ts)

    if (m <= 1) {
      cli::cli_abort(
        "Seasonal scaling requires a {.cls ts} object with frequency > 1."
      )
    }

    if (n <= m) {
      cli::cli_abort(
        "Time series is too short for seasonal scaling."
      )
    }

    Q <- mean(abs(y[(m + 1):n] - y[1:(n - m)]), na.rm = TRUE)
  }

  .tempssm_cli_debug("Scale Q computed: {round(Q, 4)}")

  return(Q)
}
