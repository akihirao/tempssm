# ---- internal helper: validate tempssm object ----
.validate_tempssm_for_ic <- function(res) {

  if (!inherits(res, "tempssm")) {
    stop("`res` must be an object of class 'tempssm'.", call. = FALSE)
  }

  if (!isTRUE(res$converged)) {
    stop(
      "Information criteria are not available because the model did not converge.",
      call. = FALSE
    )
  }

  if (is.null(res$model) || is.null(res$fit)) {
    stop(
      "Information criteria are not available because the fitted model is missing.",
      call. = FALSE
    )
  }

  invisible(res)
}



# ---- internal helper: logLik core definition ----
.internal_logLik_tempssm <- function(res) {

  .validate_tempssm_for_ic(res)

  ll <- tryCatch(
    as.numeric(stats::logLik(res$model)),
    error = function(e) {
      stop(
        "Failed to extract log-likelihood from the fitted model.",
        call. = FALSE
      )
    }
  )

  df <- length(res$fit$optim.out$par)
  if (!is.null(res$data_exogenous)) {
    df <- df + ncol(res$data_exogenous)
  }

  nobs <- length(res$data_temp)

  list(
    logLik = ll,
    df     = df,
    nobs   = nobs
  )
}








######################################
# internal helper for log-likeliihod
.get_loglik_tempssm <- function(res) {

  if (!inherits(res, "tempssm")) {
    stop("`res` must be an object of class 'tempssm'.", call. = FALSE)
  }

  loglik <- tryCatch(
    as.numeric(logLik(res$model)),
    error = function(e) {
      stop("Failed to extract log-likelihood from the fitted model.",
           call. = FALSE)
    }
  )

  loglik
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