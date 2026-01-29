#' Summary method for ThermoSSM objects
#'
#' Provides a concise summary of a fitted linear Gaussian
#' state-space model estimated by \code{lgssm()}.
#'
#' @param object An object of class \code{"ThermoSSM"} returned by \code{lgssm()}.
#' @param ... Additional arguments (currently not used).
#'
#' @return A list containing model diagnostics and summaries.
#'
#' @method summary ThermoSSM
#' @importFrom stats logLik
#' @export
summary.ThermoSSM <- function(object, ...) {

  model <- object$model
  kfs   <- object$kfs
  opt   <- object$fit$optim.out
  pars <- object$fit$optim.out$par

  exo_data <- object$data_exogenous

  if(is.null(exo_data)){
    exogenous_variable <- NULL
  }else{
    exogenous_variable <- colnames(exo_data)
  }

  exogenous_coef_ci <- extract_exo_coef_ci(object)
  k = length(opt$par) + length(exogenous_variable)
  
  res <- list(
    call        = object$call,
    logLik      = logLik(model),
    k           = k,
    AIC         = -2 * as.numeric(logLik(model)) + 2 * k,
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
    ),
    exogenous = exogenous_variable,
    exogenous_coef = exogenous_coef_ci
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
  cat("\n")
  cat("Coefficient of auto-regression parameters:\n")
  if (!is.null(x$coef_ar$AR1)) {
    cat("  AR1:", x$coef_ar$AR1, "\n")
  }
  if (!is.null(x$coef_ar$AR2)) {
    cat("  AR2:", x$coef_ar$AR2, "\n\n")
  }

  if (!is.null(x$exogenous_coef)) {
  cat("Exogenous variable\t",x$exogenous_coef$Variable, "\n")
  cat("Estimated coefficient\t", x$exogenous_coef$Coefficient, "\n")
  cat("Lower CI\t", x$exogenous_coef$lwr, "\n")
  cat("Upper CI\t", x$exogenous_coef$upr, "\n\n")
  }
  
  invisible(x)
}

