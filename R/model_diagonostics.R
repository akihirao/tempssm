#' Extract standardized recursive residuals
#'
#' @inheritParams get_level_ts
#' @return A numeric vector of standardized recursive residuals.
#'
#' @examples
#' \dontrun{
#' data(niigata_sst)
#' res <- tempssm(niigata_sst)
#'
#' residuals <- get_tempssm_residuals(res)
#' }
#' @export
get_tempssm_residuals <- function(res) {
  if (!inherits(res, "tempssm")) {
    cli::cli_abort(
      "`res` must be an object of class {.cls tempssm}."
    )
  }

  r <- stats::rstandard(res$kfs, type = "recursive")
  r[is.finite(r)]
}


#' Compute kurtosis
#'
#' Uses non-excess kurtosis, m4 / m2^2.
#'
#' @inheritParams .check_na_ratio
#' @param na.rm Logical scalar; if TRUE, missing values are removed.
#'
#' @return A numeric scalar.
#' @keywords internal
#' @noRd
.kurtosis <- function(x, na.rm = FALSE) {
  if (na.rm) {
    x <- x[!is.na(x)]
  }

  m2 <- mean((x - mean(x))^2)
  m4 <- mean((x - mean(x))^4)

  m4 / m2^2
}


#' Residual diagnostics for tempssm models
#'
#' @description
#' Compute residual diagnostic statistics for a fitted tempssm model.
#' The output is returned as a tidy tibble suitable for meta-analysis
#' across many fitted models.
#'
#' @inheritParams get_level_ts
#' @param JB_test Logical scalar; if TRUE, the Jarque–Bera test is included.
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
    cli::cli_abort(
      "`res` must be an object of class {.cls tempssm}."
    )
  }

  .tempssm_check_length_one(JB_test, "JB_test")
  .tempssm_check_logical(JB_test, "JB_test")
  if (!is.logical(JB_test) || is.na(JB_test)) {
    cli::cli_abort("{.arg JB_test} must be a logical scalar.")
  }

  r <- get_tempssm_residuals(res)

  freq <- frequency(res$data_temp)
  n_ts <- length(r)

  lb <- stats::Box.test(
    r,
    type = "Ljung-Box",
    lag = min(2 * freq, n_ts / 5)
  )

  kurt <- .kurtosis(r)

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
#' @inheritParams get_level_ts
#' @param save Logical scalar; if TRUE, plots are saved.
#' @param prefix Character scalar used as the prefix for file names.
#'
#' @return
#' Invisibly returns NULL. Called for its side effects (plots).
#'
#' @examples
#' \dontrun{
#' data(niigata_sst)
#' res <- tempssm(niigata_sst)
#'
#' model_diagnose <- diagnose_residuals(res)
#' plot_tempssm_residual_diagnostics(model_diagnose)
#' }
#' @export
plot_tempssm_residual_diagnostics <- function(res,
                                              save = FALSE,
  prefix = "residuals") {
  .tempssm_check_length_one(save, "save")
  .tempssm_check_logical(save, "save")
  if (!is.logical(save) || is.na(save)) {
    cli::cli_abort("{.arg save} must be a logical scalar.")
  }

  .tempssm_check_length_one(prefix, "prefix")
  .tempssm_check_character(prefix, "prefix")
  if (!is.character(prefix) || is.na(prefix)) {
    cli::cli_abort("{.arg prefix} must be a character scalar.")
  }

  r <- get_tempssm_residuals(res)

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
