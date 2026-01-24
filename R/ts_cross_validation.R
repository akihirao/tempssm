#' Generate train/test data of ts object for time series cross-validation
#'
#' 月次 ts（frequency = 12）を rolling-origin で学習・テストに分割します。
#'
#' @param x Monthly temperature time series of class \code{ts}.
#'   The time series must have a frequency of 12 (monthly data).
#' @param initial 学習データの初期長（観測数、例: 60）. Default value = 60
#' @param horizon テストの先読み長（観測数、例: 12）Default value = 12
#' @param step フォールド間のスライド幅（観測数、例: 1 または 12）Default value = 12
#' @param fixed_window TRUE なら固定窓（学習期間は常に initial）、FALSE なら拡張窓. Default = FALSE
#' @param allow_partial TRUE なら最後のテストが horizon に満たなくても採用. Default = FALSE
#'
#' @return 各フォールドの list。各要素は
#' \describe{
#'   \item{fold}{フォールド番号}
#'   \item{train_ts}{学習 ts}
#'   \item{test_ts}{テスト ts}
#'   \item{train_idx}{学習のインデックス範囲（位置ベース）}
#'   \item{test_idx}{テストのインデックス範囲（位置ベース）}
#'   \item{train_range}{c(start_time, end_time)}
#'   \item{test_range}{c(start_time, end_time)}
#' }
#'
#' @examples
#' # 例: AirPassengers（月次、1949–1960）
#' data(AirPassengers)
#' folds <- ts_cv_folds(
#'   AirPassengers, initial = 60, horizon = 12, step = 12, fixed_window = FALSE
#' )
#' length(folds)   # フォールド数
#' folds[[1]]$train_ts
#' folds[[1]]$test_ts
#'
#' @importFrom stats start end frequency time window
#' @export
ts_cv_folds <- function(x,
                        initial = 60,
                        horizon = 12,
                        step = 12,
                        fixed_window = FALSE,
                        allow_partial = FALSE) {
  stopifnot(inherits(x, "ts"))
  n <- length(x)
  if (initial < 1 || horizon < 1 || step < 1) {
    stop("initial, horizon, step は 1 以上で指定してください。")
  }
  if (initial >= n) stop("initial が系列長以上です。")
  if (frequency(x) != 12) {
    stop("ts object must be a monthly time series (frequency = 12).")
  }
  
  # 位置ベースで ts を切り出すユーティリティ
  ts_slice <- function(x, i_start, i_end) {
    window(x, start = stats::time(x)[i_start], end = stats::time(x)[i_end])
  }
  
  folds <- list()
  k <- 1
  train_end <- initial
  
  while (TRUE) {
    # 学習範囲
    if (fixed_window) {
      train_start <- max(1, train_end - initial + 1)
    } else {
      train_start <- 1
    }
    
    # テスト範囲
    test_start <- train_end + 1
    test_end   <- train_end + horizon
    
    # 終了条件（テストが始められない）
    if (test_start > n) break
    
    # 最終フォールドの取り扱い
    if (test_end > n) {
      if (!allow_partial) break
      test_end <- n
    }
    
    # 切り出し
    train_ts <- ts_slice(x, train_start, train_end)
    test_ts  <- ts_slice(x, test_start,  test_end)
    
    folds[[k]] <- list(
      fold = k,
      train_ts = train_ts,
      test_ts  = test_ts,
      train_idx = c(train_start, train_end),
      test_idx  = c(test_start, test_end),
      train_range = c(stats::time(x)[train_start], stats::time(x)[train_end]),
      test_range  = c(stats::time(x)[test_start],  stats::time(x)[test_end])
    )
    
    # 次のフォールドへ
    train_end <- train_end + step
    if (train_end >= n) break
    k <- k + 1
  }
  
  return(folds)
}

