#' Summary method for ThermoSSM objects
#'
#' Provides a concise summary of a fitted linear Gaussian
#' state-space model estimated by \code{lgssm()}.
#'
#' @param res An object of class \code{ThermoSSM}.
#' @param ... Additional arguments (currently not used).
#'
#' @return A list containing model diagnostics and summaries.
#'
#' @method summary ThermoSSM
#' @export
summary.ThermoSSM <- function(res, ...) {

  model <- res$fit$model
  kfs   <- res$kfs
  opt   <- res$fit$optim.out
  pars <- res$fit$optim.out$par


  res <- list(
    call        = res$call,
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



#' Print method for summary of ThermoSSM objects
#'
#' Prints a human-readable summary of a fitted
#' linear Gaussian state-space model estimated by
#' \code{lgssm()}.
#'
#' This method is automatically called when a
#' \code{summary.ThermoSSM} object is printed.
#'
#' @param x An object of class \code{summary.ThermoSSM}.
#' @param ... Additional arguments (currently not used).
#'
#' @return The input object \code{x}, invisibly.
#'
#' @method print summary.ThermoSSM
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

