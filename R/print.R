#' @export
print.ThermoSSM <- function(x, ...) {

  cat("ThermoSSM model fit\n")
  cat("------------------\n")

  # Observed time series
  ts_data <- x$data
  cat("Data:\n")
  cat("  Length      :", length(ts_data), "\n")
  cat("  Frequency   :", frequency(ts_data), "\n")
  cat("  Start / End :", 
      paste(start(ts_data), collapse = "-"), " / ",
      paste(end(ts_data), collapse = "-"), "\n\n")

  # Convergence information
  opt <- x$fit$optim.out
  cat("Optimization:\n")
  cat("  Converged   :", opt$convergence == 0, "\n")
  cat("  LogLik      :", round(logLik(x$fit$model), 2), "\n\n")

  cat("Use summary() for detailed results.\n")

  invisible(x)
}
