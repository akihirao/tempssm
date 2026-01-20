## Summary and printing methods for ThermoSSM objects

#' @export
summary.ThermoSSM <- function(object, ...) {

  model <- object$fit$model
  kfs   <- object$kfs
  opt   <- object$fit$optim.out

  res <- list(
    call        = object$call,
    logLik      = logLik(model),
    k           = length(opt$par),
    AIC         = -2 * as.numeric(logLik(model)) + 2 * length(opt$par),
    convergence = opt$convergence == 0,
    variances   = list(
      H = model$H,
      Q = model$Q
    ),
    ar          = if (!is.null(model$ar)) model$ar else NULL
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
  cat("  AIC            :", round(x$AIC, 2), "\n")
  cat("  Converged      :", x$convergence, "\n\n")

  cat("Variance parameters:\n")
  if (!is.null(x$variances$H)) {
    cat("  Observation (H):", x$variances$H, "\n")
  }
  if (!is.null(x$variances$Q)) {
    cat("  State (Q):\n")
    print(x$variances$Q)
  }

  if (!is.null(x$ar)) {
    cat("\nAR coefficients:\n")
    print(x$ar)
  }

  invisible(x)
}

