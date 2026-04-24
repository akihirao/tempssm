#' Summary method for tempssm objects
#'
#' Provides a concise summary of a fitted linear Gaussian
#' state-space model estimated by \code{lgssm()}.
#'
#' @param object An object of class \code{"tempssm"} returned by \code{lgssm()}.
#' @param ... Additional arguments (currently not used).
#'
#' @return A list containing model diagnostics and summaries.
#'
#' @method summary tempssm
#' @importFrom stats logLik
#' @export
summary.tempssm <- function(object, ...) {

  model <- object$model
  kfs   <- object$kfs
  opt   <- object$fit$optim.out
  pars <- object$fit$optim.out$par
  ar_order <- object$ar_order
  use_seasonal <- object$use_seasonal

  if(use_seasonal){
    ar_idx <- 3:(2 + ar_order)
    var_idx <- 3 + ar_order
    H_idx   <- 4 + ar_order
    
    Q_season_est <- exp(pars[2])
    
  }else{
    ar_idx <- 2:(1 + ar_order)
    var_idx <- 2 + ar_order
    H_idx   <- 3 + ar_order
    
    Q_season_est <- NA
  }

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
      H = exp(pars[H_idx]),
      Q_trend = exp(pars[1]),
      Q_season = Q_season_est,    
      Q_ar = exp(pars[var_idx])
    ),
    coef_ar = list(
      AR_order = ar_order,
      AR_coef = KFAS::artransform(pars[ar_idx])
    ),
    exogenous = exogenous_variable,
    exogenous_coef = exogenous_coef_ci
  )

  class(res) <- "summary.tempssm"
  res
}



#' Print method for summary of tempssm objects
#'
#' Prints a human-readable summary of a fitted
#' linear Gaussian state-space model estimated by
#' \code{lgssm()}.
#'
#' This method is automatically called when a
#' \code{summary.tempssm} object is printed.
#'
#' @param x An object of class \code{summary.tempssm}.
#' @param ... Additional arguments (currently not used).
#'
#' @return The input object \code{x}, invisibly.
#'
#' @method print summary.tempssm
#' @export
print.summary.tempssm <- function(x, ...) {

  cat("tempssm summary\n")
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
  cat("Components of auto-regression:\n")
  cat("  Order of AR:", x$coef_ar$AR_order, "\n")
  for(i in 1:x$coef_ar$AR_order){
    if (!is.null(x$coef_ar$AR_coef[i])) {
      cat(paste0("  Coefficient of AR",i,":"), x$coef_ar$AR_coef[i], "\n")
    }
  }

  if (!is.null(x$exogenous_coef)) {
  cat("Exogenous variable\t",x$exogenous_coef$Variable, "\n")
  cat("Estimated coefficient\t", x$exogenous_coef$Coefficient, "\n")
  cat("Lower CI\t", x$exogenous_coef$lwr, "\n")
  cat("Upper CI\t", x$exogenous_coef$upr, "\n\n")
  }
  
  invisible(x)
}

