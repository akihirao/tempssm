
#' Function to estimate a linear Gaussian state-space model without exogenous variables
#' @import KFAS
#'
#' @param ts_data time series data
#'
#' @encoding UTF-8
#'
#' @export

lgssm <- function(ts_data) {
  # モデルの構造を決める
  build_ssm <- SSModel(
    H = NA,
    Temp ~
      SSMtrend(degree = 2,                  # 平滑化トレンドモデル
               Q = c(list(0), list(NA))) +
      SSMseasonal(
        sea.type = "dummy", # ダミー変数を利用した季節成分
        period = 12,        # 周期は12とする
        Q = NA
      ) +
      SSMarima(
        ar = c(0, 0),       # 2次のAR成分
        d = 0,
        Q = 0
      ),
    data = ts_data
  )

  # optimに渡す前にパラメータをexpしたりartransformしたり、変換する
  # ほぼbuild_ssmと同じだが、パラメータだけ変更されている
  update_func <- function(pars, model) {
    model <- SSModel(
      H = exp(pars[6]),
      Temp ~
        SSMtrend(degree = 2,
                 Q = c(list(0), list(exp(pars[1])))) +
        SSMseasonal(
          sea.type = "dummy",
          period = 12,
          Q = exp(pars[2])
        ) +
        SSMarima(
          ar = artransform(pars[3:4]),
          d = 0,
          Q = exp(pars[5])
        ),
      data = ts_data
    )
  }


  # 最適化その1。まずはNelder-Mead法を用いて暫定的なパラメータを推定
  fit_ssm_bef <- fitSSM(
    build_ssm,
    #inits = c(-17,-30, 0.5, 0, -1, -3), # パラメータの初期値(任意)
    inits = c(-13,-7, 0.9, -0.1, -0.3, -5), # パラメータの初期値(任意)
    update_func,
    method = "Nelder-Mead",
    control = list(maxit = 5000, reltol = 1e-16)
  )

  # 最適化その2。先ほどの結果を初期値に使ってもう一度最適化する
  fit_ssm <- fitSSM(
    build_ssm,
    inits = fit_ssm_bef$optim.out$par,
    update_func,
    method = "BFGS",
    control = list(maxit = 5000, reltol = 1e-16)
  )

  # フィルタリングとスムージング
  result_ssm <- KFS(
    fit_ssm$model,
    filtering = c("state", "mean"),
    smoothing = c("state", "mean", "disturbance")
  )

  # 結果の出力
  return(list(fit_ssm, result_ssm))

}



#' Function to extract process/observation error parameters
#' @import KFAS
#'
#' @param res return of lgssm
#'
#' @encoding UTF-8
#'
#' @export

extract_pars <- function (res){

  pars <- res[[1]]$optim.out$par #モデルの推定パラメーター

  if (length(pars) != 6) {
    message("Invalid model: the number of parameters must be 6.")
    return(NA)
  }

  pars_comp <- c(Q_trend  = exp(pars[1]), # 年トレンドの大きさ
                  Q_season = exp(pars[2]), # 季節トレンドの大きさ
                  AR1      = KFAS::artransform(pars[3:4])[1], # 1次のARの大きさ
                  AR2      = KFAS::artransform(pars[3:4])[2], # 2次のARの大きさ
                  Q_ar     = exp(pars[5]), # 短期変動の揺らぎ
                  H        = exp(pars[6])) # 観察誤差の大きさ

  return(pars_comp)
}




#' Function to extract smoothing estimate of level component as ts object
#'
#' @param res return of lgssm
#'
#' @encoding UTF-8
#'
#' @export

extract_level <- function (res){

  smooth <- res[[2]] # smoothing
  alpha_hat <- smooth$alphahat
  level <- alpha_hat[,"level"]

  return(level)
}




#' Function to extract smoothing estimate of drift component as ts object
#'
#' @param res return of lgssm
#'
#' @encoding UTF-8
#'
#' @export

extract_drift <- function (res){

  smooth <- res[[2]] # smoothing
  alpha_hat <- smooth$alphahat
  drift <- alpha_hat[,"slope"]

  return(drift)
}

