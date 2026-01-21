## Summary and printing methods for ThermoSSM objects

#' @export
summary.ThermoSSM <- function(object, ...) {

  model <- object$fit$model
  kfs   <- object$kfs
  opt   <- object$fit$optim.out
  pars <- res$fit$optim.out$par

  #params <- c(Q_trend  = exp(pars[1]), # 年トレンドの大きさ
  #            Q_season = exp(pars[2]), # 季節トレンドの大きさ
  #            AR1      = KFAS::artransform(pars[3:4])[1], # 1次のARの大きさ
  #            AR2      = KFAS::artransform(pars[3:4])[2], # 2次のARの大きさ
  #            Q_ar     = exp(pars[5]), # 短期変動の揺らぎ
  #)

  res <- list(
    call        = object$call,
    logLik      = logLik(model),
    k           = length(opt$par),
    AIC         = -2 * as.numeric(logLik(model)) + 2 * length(opt$par),
    convergence = opt$convergence == 0,
    variances   = list(
      H = model$H,
      Q_trend = exp(pars[1]),
      Q_season = exp(pars[2]),
      Q_ar = exp(pars[5])
    ),
    coef_ar = list(
      AR1 = KFAS::artransform(pars[3:4])[1],
      AR2 = KFAS::artransform(pars[3:4])[2]
    )

  )

  class(res) <- "summary.ThermoSSM"
  res
}




#' @export
print.summary.ThermoSSM <- function(x, ...) {

  cat("ThermoSSM summary\n")
  cat("-----------------\n")

  cat("Call:\n")
  print(x$call)
  cat("\n")

  cat("Model fit:\n")
  cat("  Log-likelihood :", round(x$logLik, 2), "\n")
  cat("  k              :", x$k, "\n")
  cat("  AIC            :", round(x$AIC, 2), "\n")
  cat("  Converged      :", x$convergence, "\n\n")

  cat("Variance parameters:\n")
  if (!is.null(x$variances$H)) {
    cat("  Observation (H):", x$variances$H, "\n")
  }
  if (!is.null(x$variances$Q_trend)) {
    cat("  State (Q trend):", x$variances$Q_trend, "\n")
  }
  if (!is.null(x$variances$Q_season)) {
    cat("  State (Q season):", x$variances$Q_season, "\n")
  }
  if (!is.null(x$variances$Q_ar)) {
    cat("  State (Q ar):", x$variances$Q_ar, "\n")
  }
  cat("Coefficient of auto-regression parameters:\n")
  if (!is.null(x$coef_ar$AR1)) {
    cat("  AR1:", x$coef_ar$AR1, "\n")
  }
  if (!is.null(x$coef_ar$AR2)) {
    cat("  AR2:", x$coef_ar$AR2, "\n")
  }

  invisible(x)
}

