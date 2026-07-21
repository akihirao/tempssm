######################################
#' @keywords internal
#' @noRd
.validate_tempssm_for_logLik <- function(res) {
  ## ---- class check ----------------------------------------------------
  if (!inherits(res, "tempssm")) {
    cli::cli_abort(
      "`res` must be an object of class {.cls tempssm}."
    )
  }

  ## ---- convergence check ----------------------------------------------
  if (!isTRUE(res$converged)) {
    cli::cli_abort(
      "Log-likelihood is not available because the model did not converge."
    )
  }

  ## ---- component check ------------------------------------------------
  if (is.null(res$model) || is.null(res$fit)) {
    cli::cli_abort(
      "Log-likelihood is not available because the fitted model is missing."
    )
  }

  ## ---- debug message --------------------------------------------------
  .tempssm_cli_debug("Validated tempssm object for log-likelihood")

  return(invisible(res))
}


#' Resolve the likelihood setting for a fitted model
#'
#' @inheritParams get_level_ts
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
#' @inheritParams get_level_ts
#' @inheritParams logLik.tempssm
#'
#' @return A named list containing the log-likelihood, degrees of freedom,
#'   number of observations, and resolved likelihood setting.
#'
#' @keywords internal
#' @noRd
.internal_logLik_tempssm <- function(res, marginal = NULL) {
  ## ---- validation ----------------------------------------------------
  .validate_tempssm_for_logLik(res)
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
#' AIC is intentionally not computed for models fitted by \code{tempssm()}.
#' The method is registered to prevent automatic AIC calculation from the
#' \code{logLik()} method.
#'
#' @param object
#' An object of class \code{"tempssm"}, typically returned by
#' \code{tempssm()}.
#'
#' @param ...
#' Additional arguments passed to the generic \code{AIC()} function.
#' These are currently ignored.
#'
#' @param k
#' Numeric penalty coefficient accepted for compatibility with
#' \code{\link[stats]{AIC}}. This argument is ignored.
#'
#' @param marginal
#' Logical scalar or \code{NULL} accepted for backward compatibility.
#' This argument is ignored.
#'
#' @return
#' This function always raises an error.
#'
#' @method AIC tempssm
#' @export
AIC.tempssm <- function(object, ..., k = 2, marginal = NULL) {
  cli::cli_abort(
    c(
      "AIC is not computed for {.cls tempssm} objects.",
      "i" = "Use {.fn logLik} to extract the log-likelihood.",
      "i" = "Use {.code attr(logLik(x), \"df\")} for the parameter count.",
      "i" = "If needed, compute AIC explicitly under your own assumptions."
    )
  )
}


#' Deprecated AIC helper for tempssm objects
#'
#' @description
#' This function is retained for backward compatibility but is deprecated.
#' AIC is intentionally not computed for \code{tempssm} objects.
#'
#' @inheritParams get_level_ts
#'
#' @param marginal
#' Logical scalar or \code{NULL} accepted for backward compatibility.
#' This argument is ignored.
#'
#' @return
#' This function always raises an error.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(niigata_sst)
#' res <- tempssm(niigata_sst)
#' logLik(res)
#' }
get_aic <- function(res, marginal = NULL) {
  cli::cli_abort(
    c(
      "{.fn get_aic} is deprecated and no longer computes AIC.",
      "i" = "AIC is no longer computed by tempssm.",
      "i" = "Use {.fn logLik} and {.code attr(logLik(x), \"df\")} instead."
    )
  )
}
