######################################
#' @keywords internal
#' @noRd
.validate_tempssm_for_ic <- function(res) {
  ## ---- class check ----------------------------------------------------
  if (!inherits(res, "tempssm")) {
    cli::cli_abort(
      "`res` must be an object of class {.cls tempssm}."
    )
  }

  ## ---- convergence check ----------------------------------------------
  if (!isTRUE(res$converged)) {
    cli::cli_abort(
      paste(
        "Information criteria are not available",
        "because the model did not converge."
      )
    )
  }

  ## ---- component check ------------------------------------------------
  if (is.null(res$model) || is.null(res$fit)) {
    cli::cli_abort(
      paste(
        "Information criteria are not available",
        "because the fitted model is missing."
      )
    )
  }

  ## ---- debug message --------------------------------------------------
  .tempssm_cli_debug("Validated tempssm object for information criteria")

  return(invisible(res))
}


######################################
#' @keywords internal
#' @noRd
.internal_logLik_tempssm <- function(res) {
  ## ---- validation ----------------------------------------------------
  .validate_tempssm_for_ic(res)

  ## ---- logLik extraction ---------------------------------------------
  ll <- tryCatch(
    as.numeric(stats::logLik(res$model)),
    error = function(e) {
      cli::cli_abort(
        "Failed to extract log-likelihood from the fitted model."
      )
    }
  )

  ## ---- degrees of freedom --------------------------------------------
  df <- length(res$fit$optim.out$par)

  if (!is.null(res$exogenous_data)) {
    df <- df + ncol(res$exogenous_data)
  }

  ## ---- number of observations ----------------------------------------
  nobs <- length(res$temp_data)

  ## ---- debug message --------------------------------------------------
  .tempssm_cli_debug(
    "Computed logLik components: logLik={round(ll, 3)}, df={df}, nobs={nobs}"
  )

  ## ---- return --------------------------------------------------------
  return(list(
    logLik = ll,
    df     = df,
    nobs   = nobs
  ))
}


#' Log-likelihood method for tempssm objects (S3 method)
#'
#' @param object
#' An object of class \code{"tempssm"}, typically returned by
#' \code{tempssm()}.
#'
#' @param ...
#' Additional arguments passed to the generic \code{logLik()} function.
#' These are currently ignored but are included for compatibility
#' with the generic interface.
#'
#' @method logLik tempssm
#'
#' @return
#' A numeric value representing the log-likelihood.
#'
#' @seealso
#' \code{\link{AIC.tempssm}} for computing AIC,
#' \code{\link{get_aic}} as a convenience wrapper.
#'
#' @examples
#' \dontrun{
#' data(niigata_sst)
#'
#' # fit model
#' res <- tempssm(niigata_sst)
#'
#' # extract log-likelihood
#' ll <- logLik(res)
#'
#' ll
#'
#' # access attributes
#' attr(ll, "df") # number of parameters
#' attr(ll, "nobs") # number of observations
#' }
#'
#' @export
logLik.tempssm <- function(object, ...) {
  info <- .internal_logLik_tempssm(object)

  structure(
    info$logLik,
    class = "logLik",
    df    = info$df,
    nobs  = info$nobs
  )
}


#' AIC method for tempssm objects (S3 method)
#'
#' @description
#' Compute the Akaike Information Criterion (AIC) for a model fitted
#' by \code{tempssm()}. This method extends the generic
#' \code{\link[stats]{AIC}} function.
#'
#' @param object
#' An object of class \code{"tempssm"}, typically returned by
#' \code{tempssm()}.
#'
#' @param ...
#' Additional arguments passed to the generic \code{AIC()} function.
#' These are currently ignored but are included for compatibility
#' with the generic interface.
#'
#' @param k
#' Numeric penalty coefficient for the number of parameters.
#' This argument is included for compatibility with
#' \code{\link[stats]{AIC}} but is not used in the \pkg{tempssm} method,
#' where the standard AIC definition (\code{k = 2}) is applied.
#'
#' @return
#' A numeric value giving the AIC of the fitted \code{tempssm} model.
#'
#' @examples
#' \dontrun{
#' data(niigata_sst)
#' res <- tempssm(niigata_sst)
#' aic <- AIC(res)
#' }
#'
#' @method AIC tempssm
#' @export
AIC.tempssm <- function(object, ..., k = 2) {
  ll <- logLik.tempssm(object)
  -2 * as.numeric(ll) + k * attr(ll, "df")
}


#' Extract the Akaike Information Criterion (AIC)
#'
#' @description
#' Compute the Akaike Information Criterion (AIC) for a fitted
#' \code{tempssm} model using the model log-likelihood and the
#' number of estimated parameters.
#'
#' @param res An object of class \code{"tempssm"} returned by \code{tempssm()}.
#'
#' @details
#' The number of parameters is determined from the optimization results
#' stored in the fitted model. If exogenous variables are included,
#' their coefficients are added to the parameter count.
#'
#' @return
#' A single numeric value representing the AIC of the fitted model.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(niigata_sst)
#' res <- tempssm(niigata_sst)
#' aic <- get_aic(res)
#' }
get_aic <- function(res) {
  AIC.tempssm(res)
}
