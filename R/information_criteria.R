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


#' Resolve the likelihood setting for a fitted model
#'
#' @param res A fitted \code{tempssm} object.
#' @inheritParams logLik.tempssm
#'
#' @return A logical scalar.
#'
#' @keywords internal
#' @noRd
.resolve_tempssm_marginal <- function(res, marginal = NULL) {
  if (is.null(marginal)) {
    # Objects created before marginal was stored used KFAS's FALSE default.
    return(isTRUE(res$marginal))
  }

  .validate_marginal(marginal)
  marginal
}


######################################
#' Extract log-likelihood components from a fitted model
#'
#' @param res A fitted \code{tempssm} object.
#' @inheritParams logLik.tempssm
#'
#' @return A named list containing the log-likelihood, degrees of freedom,
#'   number of observations, and resolved likelihood setting.
#'
#' @keywords internal
#' @noRd
.internal_logLik_tempssm <- function(res, marginal = NULL) {
  ## ---- validation ----------------------------------------------------
  .validate_tempssm_for_ic(res)
  marginal <- .resolve_tempssm_marginal(res, marginal)

  ## ---- logLik extraction ---------------------------------------------
  ll <- tryCatch(
    # stats::logLik() dispatches SSModel objects to KFAS's registered
    # logLik.SSModel method; KFAS does not export the method directly.
    as.numeric(stats::logLik(res$model, marginal = marginal)),
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
    nobs   = nobs,
    marginal = marginal
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
#' These are currently ignored.
#'
#' @param marginal Logical scalar or \code{NULL}. If \code{NULL}, use the
#' likelihood setting stored when the model was fitted. Set to \code{TRUE} or
#' \code{FALSE} to evaluate the fitted parameters with the marginal or diffuse
#' likelihood, respectively. An explicit value does not refit the model.
#'
#' @method logLik tempssm
#'
#' @details
#' The implementation calls the public \code{stats::logLik()} generic on the
#' fitted \code{SSModel} object. S3 dispatch then invokes KFAS's registered
#' \code{logLik.SSModel} method; that method is not exported for direct use.
#' The \code{marginal} argument is passed to the KFAS method.
#'
#' @return
#' An object of class \code{"logLik"} containing the numeric log-likelihood,
#' with \code{df} and \code{nobs} attributes.
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
#' marginal_ll <- logLik(res, marginal = TRUE)
#'
#' ll
#'
#' # access attributes
#' attr(ll, "df") # number of parameters
#' attr(ll, "nobs") # number of observations
#' }
#'
#' @export
logLik.tempssm <- function(object, ..., marginal = NULL) {
  info <- .internal_logLik_tempssm(object, marginal = marginal)

  structure(
    info$logLik,
    class = "logLik",
    df    = info$df,
    nobs  = info$nobs,
    marginal = info$marginal
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
#' These are currently ignored.
#'
#' @inheritParams logLik.tempssm
#'
#' @param k
#' Numeric penalty coefficient for the number of parameters.
#' The default \code{k = 2} gives the standard AIC definition.
#'
#' @return
#' A numeric scalar giving the AIC of the fitted \code{tempssm} model.
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
AIC.tempssm <- function(object, ..., k = 2, marginal = NULL) {
  ll <- logLik.tempssm(object, marginal = marginal)
  -2 * as.numeric(ll) + k * attr(ll, "df")
}


#' Extract the Akaike Information Criterion (AIC)
#'
#' @description
#' Compute the Akaike Information Criterion (AIC) for a fitted
#' \code{tempssm} model using the model log-likelihood and the
#' number of estimated parameters.
#'
#' @inheritParams get_level_ts
#' @inheritParams logLik.tempssm
#'
#' @details
#' The number of parameters is determined from the optimization results
#' stored in the fitted model. If exogenous variables are included,
#' their coefficients are added to the parameter count.
#'
#' @return
#' A numeric scalar representing the AIC of the fitted model.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(niigata_sst)
#' res <- tempssm(niigata_sst)
#' aic <- get_aic(res)
#' }
get_aic <- function(res, marginal = NULL) {
  AIC.tempssm(res, marginal = marginal)
}
