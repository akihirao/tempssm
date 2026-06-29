#' Validate input for a tempssm summary
#'
#' @param object An object expected to contain fitted \code{tempssm} results.
#'
#' @return Invisibly returns \code{object}.
#'
#' @keywords internal
#' @noRd
.validate_tempssm_summary_input <- function(object) {
  if (!inherits(object, "tempssm")) {
    cli::cli_abort("`object` must be an object of class {.cls tempssm}.")
  }

  if (is.null(object$fit) || is.null(object$fit$optim.out) ||
      is.null(object$fit$optim.out$par)) {
    cli::cli_abort(
      "Summary is not available because the fitted model results are missing."
    )
  }

  invisible(object)
}


#' Summary method for tempssm objects
#'
#' Provides a concise summary of a fitted linear Gaussian
#' state-space model estimated by \code{tempssm()}.
#'
#' @param object An object of class \code{"tempssm"} returned
#' by \code{tempssm()}.
#' @param ... Additional arguments (currently not used).
#'
#' @return
#' An object of class \code{"summary.tempssm"}, implemented as a named list
#' with components \code{call}, \code{logLik}, \code{k}, \code{AIC},
#' \code{convergence}, \code{variances}, \code{coef_ar},
#' \code{exogenous}, and \code{exogenous_coef}.
#'
#' @method summary tempssm
#' @importFrom stats logLik
#'
#' @examples
#' \dontrun{
#' data(niigata_sst)
#'
#' # fit model
#' res <- tempssm(niigata_sst)
#'
#' # compute summary
#' s <- summary(res)
#'
#' s
#'
#' # access components programmatically
#' s$AIC
#' s$variances
#' }
#'
#' @export
summary.tempssm <- function(object, ...) {
  .validate_tempssm_summary_input(object)

  opt <- object$fit$optim.out
  params <- .extract_tempssm_params(object)
  ll <- logLik.tempssm(object)
  exo_data <- object$exogenous_data
  exogenous_variable <- if (is.null(exo_data)) NULL else colnames(exo_data)

  res <- list(
    call = object$call,
    logLik = as.numeric(ll),
    k = attr(ll, "df"),
    AIC = AIC.tempssm(object),
    convergence = opt$convergence == 0,
    variances = params[c("H", "Q_trend", "Q_season", "Q_ar")],
    coef_ar = list(
      AR_order = object$ar_order,
      AR_coef = params$ARs
    ),
    exogenous = exogenous_variable,
    exogenous_coef = get_exo_coef(object)
  )

  class(res) <- "summary.tempssm"
  res
}


#' Print method for summary of tempssm objects
#'
#' Prints a human-readable summary of a fitted
#' linear Gaussian state-space model estimated by
#' \code{tempssm()}.
#'
#' This method is automatically called when a
#' \code{summary.tempssm} object is printed.
#'
#' @param x An object of class \code{summary.tempssm}.
#' @param ... Additional arguments (currently not used).
#'
#' @return
#' The input object \code{x}, invisibly. The returned object has class
#' \code{"summary.tempssm"}.
#'
#' @method print summary.tempssm
#'
#' @examples
#' \dontrun{
#' data(niigata_sst)
#'
#' # fit model and compute summary
#' res <- tempssm(niigata_sst)
#' s <- summary(res)
#'
#' # print summary (explicit)
#' print(s)
#'
#' # equivalent (implicit print)
#' s
#' }
#'
#' @export
print.summary.tempssm <- function(x, ...) {
  cat("tempssm summary\n")
  cat("-----------------\n")

  cat("Call:\n")
  print(x$call)
  cat("\n")

  cat("Model fit:\n")
  cat("  Log-likelihood :", round(x$logLik, 2), "\n")
  cat("  k              :", x$k, "\n")
  cat("  AIC            :", round(x$AIC, 2), "\n")
  cat("  Converged      :", x$convergence, "\n\n")

  cat("Variance parameters:\n")
  if (!is.null(x$variances$H)) {
    cat("  Observation (H):", x$variances$H, "\n")
  }
  if (!is.null(x$variances$Q_trend)) {
    cat("  State (Q trend):", x$variances$Q_trend, "\n")
  }
  if (!is.null(x$variances$Q_season)) {
    cat("  State (Q season):", x$variances$Q_season, "\n")
  }
  if (!is.null(x$variances$Q_ar)) {
    cat("  State (Q ar):", x$variances$Q_ar, "\n")
  }
  cat("\n")
  cat("Components of auto-regression:\n")
  cat("  Order of AR:", x$coef_ar$AR_order, "\n")
  for (i in 1:x$coef_ar$AR_order) {
    if (!is.null(x$coef_ar$AR_coef[i])) {
      cat(paste0("  Coefficient of AR", i, ":"), x$coef_ar$AR_coef[i], "\n")
    }
  }

  if (!is.null(x$exogenous_coef)) {
    cat("Exogenous variable\t", x$exogenous_coef$Variable, "\n")
    cat("Estimated coefficient\t", x$exogenous_coef$Coefficient, "\n")
    cat("Lower CI\t", x$exogenous_coef$lwr, "\n")
    cat("Upper CI\t", x$exogenous_coef$upr, "\n\n")
  }

  invisible(x)
}
