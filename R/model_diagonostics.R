#' Extract standardized recursive residuals
#'
#' @param res An object of class "tempssm".
#' @return A numeric vector of standardized recursive residuals.
#'
#' @examples
#' \dontrun{
#' data(niigata_sst)
#' res <- tempssm(niigata_sst)
#'
#' residuals <- get_residuals(res)
#' }
#' @export
get_residuals <- function(res) {

  if (!inherits(res, "tempssm")) {
    stop("`res` must be an object of class 'tempssm'.", call. = FALSE)
  }

  r <- stats::rstandard(res$kfs, type = "recursive")
  r[is.finite(r)]
}



#' Residual diagnostics for tempssm models
#'
#' @description
#' Compute residual diagnostic statistics for a fitted tempssm model.
#' The output is returned as a tidy tibble suitable for meta-analysis
#' across many fitted models.
#'
#' @param res An object of class \code{"tempssm"} returned by \code{tempssm()}.
#' @param JB_test Logical; if TRUE, the Jarque–Bera test is included.
#'
#' @return
#' A \code{tibble} with one row containing residual diagnostic statistics,
#' including Ljung--Box test results and kurtosis. If \code{JB_test = TRUE},
#' Jarque--Bera test statistics are also included.
#'
#' @examples
#' \dontrun{
#' data(niigata_sst)
#'
#' # Fit model
#' res <- tempssm(niigata_sst)
#'
#' # Residual diagnostics (tibble output)
#' diag <- diagnose_residuals(res)
#'
#' diag
#' }
#' @export
diagnose_residuals <- function(res, JB_test = FALSE) {

  if (!inherits(res, "tempssm")) {
    stop("`res` must be an object of class 'tempssm'.", call. = FALSE)
  }

  r <- get_residuals(res)

  freq <- frequency(res$data_temp)
  n_ts <- length(r)

  lb <- stats::Box.test(
    r,
    type = "Ljung-Box",
    lag = min(2 * freq, n_ts / 5)
  )

  kurt <- moments::kurtosis(r)

  if (JB_test) {
    jb <- tseries::jarque.bera.test(r)

    return(
      tibble::tibble(
        lb_stat    = unname(lb$statistic),
        lb_df      = unname(lb$parameter),
        lb_pvalue  = lb$p.value,
        kurtosis   = kurt,
        jb_stat    = unname(jb$statistic),
        jb_pvalue  = jb$p.value
      )
    )
  }

  tibble::tibble(
    lb_stat    = unname(lb$statistic),
    lb_df      = unname(lb$parameter),
    lb_pvalue  = lb$p.value,
    kurtosis   = kurt
  )
}



#' Plot residual diagnostics for tempssm models
#'
#' @param res An object of class "tempssm".
#' @param save Logical; if TRUE, plots are saved.
#' @param prefix Character prefix for file names.
#'
#' @examples
#' \dontrun{
#' data(niigata_sst)
#' res <- tempssm(niigata_sst)
#'
#' model_diagnose <- diagnose_residuals(res)
#' plot_residual_diagnostics(model_diagnose)
#' }
#' @export
plot_residual_diagnostics <- function(res,
                                      save = FALSE,
                                      prefix = "residuals") {

  r <- get_residuals(res)

  # Standard residual plot
  forecast::checkresiduals(r, test = FALSE)

  if (save) {
    grDevices::png(paste0(prefix, "_check.png"), 600, 400)
    forecast::checkresiduals(r, test = FALSE)
    grDevices::dev.off()

    grDevices::png(paste0(prefix, "_qq.png"), 600, 400)
    stats::qqnorm(r)
    stats::qqline(r)
    grDevices::dev.off()
  }

  invisible(NULL)
}
