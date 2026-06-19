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
#'   equal to \code{max_width}. If \code{pred} is a \code{ts} object, the
#'   returned object preserves its time scale.
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
  if (!is.numeric(max_width) ||
      length(max_width) != 1 ||
      is.na(max_width) ||
      !is.finite(max_width) ||
      max_width < 0) {
    cli::cli_abort(
      "{.arg max_width} must be a single non-negative finite number."
    )
  }

  if (is.null(dim(pred))) {
    cli::cli_abort("{.arg pred} must be a matrix-like object.")
  }

  pred_names <- colnames(pred)
  if (is.null(pred_names) || !all(c("lwr", "upr") %in% pred_names)) {
    cli::cli_abort(
      "{.arg pred} must contain columns named {.val lwr} and {.val upr}."
    )
  }

  interval_width <- as.numeric(pred[, "upr"] - pred[, "lwr"])
  keep <- interval_width <= max_width
  first_exceed <- match(FALSE, keep)

  if (is.na(first_exceed)) {
    return(pred)
  }

  keep_n <- first_exceed - 1L
  if (inherits(pred, "ts") && keep_n > 0L) {
    return(stats::window(pred, end = stats::time(pred)[keep_n]))
  }

  pred[seq_len(keep_n), , drop = FALSE]
}
