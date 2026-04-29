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
#' @importFrom stats logLik
logLik.tempssm <- function(object, ...) {

  # extract log-likelihood
  ll <- .get_loglik_tempssm(object)

  # degrees of freedom (number of estimated parameters)
  if (is.null(object$fit$optim.out$par)) {
    stop("Optimization results not found in the fitted model.",
         call. = FALSE)
  }

  df <- length(object$fit$optim.out$par)

  if (!is.null(object$data_exogenous)) {
    df <- df + ncol(object$data_exogenous)
  }

  # number of observations
  nobs <- length(object$data_temp)

  structure(
    ll,
    class = "logLik",
    df = df,
    nobs = nobs
  )
}



# internal helper for AIC
.compute_aic_tempssm <- function(res) {

  if (!inherits(res, "tempssm")) {
    stop("`res` must be an object of class 'tempssm'.", call. = FALSE)
  }

  if (is.null(res$fit$optim.out$par)) {
    stop("Optimization results not found in the fitted model.",
         call. = FALSE)
  }

  k <- length(res$fit$optim.out$par)

  if (!is.null(res$data_exogenous)) {
    k <- k + ncol(res$data_exogenous)
  }

  loglik <- tryCatch(
    as.numeric(logLik(res$model)),
    error = function(e) {
      stop("Failed to extract log-likelihood from the fitted model.",
           call. = FALSE)
    }
  )

  -2 * loglik + 2 * k
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
  .compute_aic_tempssm(res)
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
#' @method AIC tempssm
#' @export
#' @importFrom stats AIC
AIC.tempssm <- function(object, ..., k = 2) {
  .compute_aic_tempssm(object)
}



