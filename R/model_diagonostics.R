#' Residual diagnostics for ThermoSSM models
#'
#' Wrapper function for residual diagnostics using functions from
#' the forecast, stats, tseries, and moments packages.
#'
#' @param res An object of class \code{"ThermoSSM"} returned by \code{lgssm()}.
#' @param plot_resid_save Logical; if TRUE, the residual diagnostic plot is saved.
#' @param resid_file_name Character; file name for saving the residual diagnostic plot (default: "checkresiduals.png").
#' @param plot_qq_save Logical; if TRUE, the Q-Q plot of residuals is saved.
#' @param qq_file_name Character; file name for saving the Q-Q plot (default: "qqplot.png").
#' 
#' @return A named \code{list} containing results of residual tests.
#'
#' @export
wrapper_checkresiduals <- function(res,
                                   plot_resid_save = FALSE,
                                   resid_file_name = "checkresiduals.png",
                                   plot_qq_save = FALSE,
                                   qq_file_name = "qqplot.png"
                                   ) {
  
  freq <- frequency(res$data_temp)
  n_ts <- length(res$data_temp)
  # Standardized recursive residuals
  std_obs_resid <- stats::rstandard(res$kfs, type = "recursive")
  
  # Remove NA/Inf values
  std_obs_resid <- std_obs_resid[is.finite(std_obs_resid)]
  
  # ---- Residual diagnostics plot ----
  forecast::checkresiduals(std_obs_resid,
                           test=FALSE)
  
  p <- recordPlot()
  
  if (plot_resid_save) {
    png(resid_file_name, width = 600, height = 400)
    replayPlot(p)
    dev.off()
  }

  if (plot_qq_save) {
    png(qq_file_name, width = 600, height = 400)
    qqnorm(std_obs_resid)
    dev.off()
  }
  
  # ---- Statistical tests ----
  
  # Ljung–Box test
  Ljung_Box_test <- Box.test(std_obs_resid,
                             type = "Ljung-Box",
                             lag = min(2*freq, n_ts/5))
  
  # Jarque–Bera test
  Jarque_Bera_test <- tseries::jarque.bera.test(std_obs_resid)
  
  # 
  kurtosi <- moments::kurtosis(std_obs_resid)
  
  residuals_test_list <- list(
    Ljung_Box = Ljung_Box_test,
    Jarque_Bera = Jarque_Bera_test,
    kurtosi = kurtosi
  )
  
  return(residuals_test_list)
}
