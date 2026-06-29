#' Trim forecast intervals by maximum interval width
#'
#' Trims forecast or prediction intervals to the longest leading forecast
#' horizon for which the interval width remains within a user-specified
#' threshold.
#'
#' @param pred A matrix-like object, data frame, or multivariate \code{ts}
#'   object containing prediction interval columns named \code{"lwr"} and
#'   \code{"upr"}. This matches the object returned by
#'   \code{stats::predict(..., interval = "prediction")} for the KFAS model
#'   used by \code{tempssm()}.
#' @param max_width Numeric scalar giving the maximum permitted interval width,
#'   calculated as \code{upr - lwr}.
#'
#' @return An object of the same basic class as \code{pred}, truncated to the
#'   longest leading sequence whose prediction interval width is less than or
#'   equal to \code{max_width}. If \code{pred} is a \code{ts} object and at
#'   least one row is retained, the returned object preserves its time scale.
#'   A zero-row result is returned as a matrix because base R does not support
#'   zero-length \code{ts} objects.
#'
#' @details
#' Forecast uncertainty in state-space models generally increases with forecast
#' horizon because future observations depend on accumulated state and
#' observation uncertainty. This helper provides a simple post-processing rule:
#' retain forecasts only until their prediction interval width exceeds a
#' user-defined error margin. It does not modify forecast means, lower bounds,
#' or upper bounds; it only truncates the returned forecast horizon.
#'
#' @examples
#' \dontrun{
#' data(niigata_sst)
#' res <- tempssm(niigata_sst, na_action = "allow")
#' pred <- stats::predict(
#'   res$model,
#'   n.ahead = 24,
#'   interval = "prediction",
#'   level = 0.95
#' )
#'
#' trim_prediction_intervals(pred, max_width = 4)
#' }
#'
#' @export
trim_prediction_intervals <- function(pred, max_width) {
  .validate_max_interval_width(max_width)
  interval_width <- .validate_prediction_intervals(pred)
  keep_n <- .prediction_interval_keep_n(interval_width, max_width)

  .slice_prediction_intervals(pred, keep_n)
}


#' Validate the maximum permitted prediction interval width
#'
#' @param max_width Numeric scalar giving the maximum interval width.
#'
#' @return Invisibly returns \code{max_width}.
#'
#' @keywords internal
#' @noRd
.validate_max_interval_width <- function(max_width) {
  valid <- is.numeric(max_width) &&
    length(max_width) == 1L &&
    !is.na(max_width) &&
    is.finite(max_width) &&
    max_width >= 0

  if (!valid) {
    cli::cli_abort(
      "{.arg max_width} must be a single non-negative finite number."
    )
  }

  invisible(max_width)
}


#' Validate prediction intervals and calculate their widths
#'
#' @param pred A matrix-like prediction interval object.
#'
#' @return A numeric vector containing interval widths.
#'
#' @keywords internal
#' @noRd
.validate_prediction_intervals <- function(pred) {
  .validate_prediction_interval_structure(pred)

  lower <- pred[, "lwr"]
  upper <- pred[, "upr"]
  .validate_prediction_interval_types(lower, upper)
  .validate_prediction_interval_finiteness(lower, upper)
  .validate_prediction_interval_order(lower, upper)

  as.numeric(upper - lower)
}


#' Validate the structure of a prediction interval object
#'
#' @param pred A matrix-like prediction interval object.
#'
#' @return Invisibly returns \code{pred}.
#'
#' @keywords internal
#' @noRd
.validate_prediction_interval_structure <- function(pred) {
  if (is.null(dim(pred))) {
    cli::cli_abort("{.arg pred} must be a matrix-like object.")
  }

  pred_names <- colnames(pred)
  if (is.null(pred_names) || !all(c("lwr", "upr") %in% pred_names)) {
    cli::cli_abort(
      "{.arg pred} must contain columns named {.val lwr} and {.val upr}."
    )
  }

  invisible(pred)
}


#' Validate prediction interval column types
#'
#' @param lower Numeric vector of lower interval bounds.
#' @param upper Numeric vector of upper interval bounds.
#'
#' @return Invisibly returns \code{NULL}.
#'
#' @keywords internal
#' @noRd
.validate_prediction_interval_types <- function(lower, upper) {
  if (!is.numeric(lower) || !is.numeric(upper)) {
    cli::cli_abort(
      "The {.val lwr} and {.val upr} columns of {.arg pred} must be numeric."
    )
  }

  invisible(NULL)
}


#' Validate that prediction interval bounds are finite
#'
#' @inheritParams .validate_prediction_interval_types
#'
#' @return Invisibly returns \code{NULL}.
#'
#' @keywords internal
#' @noRd
.validate_prediction_interval_finiteness <- function(lower, upper) {
  if (anyNA(lower) || anyNA(upper) ||
      !all(is.finite(lower)) || !all(is.finite(upper))) {
    cli::cli_abort(
      "Columns {.val lwr} and {.val upr} in {.arg pred} must be finite."
    )
  }

  invisible(NULL)
}


#' Validate prediction interval bound ordering
#'
#' @inheritParams .validate_prediction_interval_types
#'
#' @return Invisibly returns \code{NULL}.
#'
#' @keywords internal
#' @noRd
.validate_prediction_interval_order <- function(lower, upper) {
  if (any(lower > upper)) {
    cli::cli_abort(
      "Each {.val lwr} value must be less than or equal to {.val upr}."
    )
  }

  invisible(NULL)
}


#' Find the length of the leading acceptable prediction horizon
#'
#' @param interval_width Numeric vector of prediction interval widths.
#' @inheritParams trim_prediction_intervals
#'
#' @return An integer giving the number of leading rows to retain.
#'
#' @keywords internal
#' @noRd
.prediction_interval_keep_n <- function(interval_width, max_width) {
  first_exceed <- match(TRUE, interval_width > max_width)

  if (is.na(first_exceed)) {
    return(length(interval_width))
  }

  first_exceed - 1L
}


#' Slice a prediction interval object to a leading horizon
#'
#' @inheritParams trim_prediction_intervals
#' @param keep_n Number of leading rows to retain.
#'
#' @return A prediction interval object truncated to \code{keep_n} rows.
#'
#' @keywords internal
#' @noRd
.slice_prediction_intervals <- function(pred, keep_n) {
  if (keep_n == NROW(pred)) {
    return(pred)
  }

  if (inherits(pred, "ts") && keep_n > 0L) {
    return(stats::window(pred, end = stats::time(pred)[keep_n]))
  }

  pred[seq_len(keep_n), , drop = FALSE]
}
