#' Check that an argument has length one
#'
#' @param x Object to check.
#' @param arg_name Name of the argument being checked.
#' @param allow_null Logical scalar indicating whether \code{NULL} is allowed.
#'
#' @return Invisibly returns \code{x}.
#'
#' @keywords internal
#' @noRd
.tempssm_check_length_one <- function(x, arg_name, allow_null = FALSE) {
  if (is.null(x) && isTRUE(allow_null)) {
    return(invisible(x))
  }

  if (length(x) != 1L) {
    cli::cli_abort(
      "{.arg {arg_name}} must have length one."
    )
  }

  invisible(x)
}


#' Check that an argument is numeric
#'
#' @inheritParams .tempssm_check_length_one
#'
#' @return Invisibly returns \code{x}.
#'
#' @keywords internal
#' @noRd
.tempssm_check_numeric <- function(x, arg_name, allow_null = FALSE) {
  if (is.null(x) && isTRUE(allow_null)) {
    return(invisible(x))
  }

  if (!is.numeric(x)) {
    cli::cli_abort(
      "{.arg {arg_name}} must be numeric."
    )
  }

  invisible(x)
}


#' Check that numeric-like values are not undefined
#'
#' @inheritParams .tempssm_check_length_one
#'
#' @return Invisibly returns \code{x}.
#'
#' @keywords internal
#' @noRd
.tempssm_check_no_undefined <- function(x, arg_name) {
  if (any(is.nan(x)) || any(is.infinite(x))) {
    cli::cli_abort(
      c(
        "{.arg {arg_name}} must not contain {.val NaN},",
        "{.val Inf}, or {.val -Inf}."
      )
    )
  }

  invisible(x)
}


#' Check whether numeric values are integer-like
#'
#' @param x Numeric vector.
#'
#' @return A logical vector indicating whether each value is integer-like.
#'
#' @keywords internal
#' @noRd
.tempssm_is_integerish <- function(x) {
  is.finite(x) & abs(x - round(x)) <= sqrt(.Machine$double.eps)
}


#' Validate an autoregressive order
#'
#' @inheritParams tempssm
#'
#' @return Invisibly returns \code{ar_order}.
#'
#' @keywords internal
#' @noRd
.validate_ar_order <- function(ar_order) {
  .tempssm_check_length_one(ar_order, "ar_order")
  .tempssm_check_numeric(ar_order, "ar_order")

  valid <- !is.na(ar_order) &&
    ar_order >= 1 &&
    .tempssm_is_integerish(ar_order)
  if (!valid) {
    cli::cli_abort(
      "The argument {.arg ar_order} must be an integer >= 1."
    )
  }

  invisible(ar_order)
}


#' Validate whether a seasonal component should be used
#'
#' @inheritParams tempssm
#'
#' @return Invisibly returns \code{use_season}.
#'
#' @keywords internal
#' @noRd
.validate_use_season <- function(use_season) {
  .tempssm_check_length_one(use_season, "use_season")
  .tempssm_check_logical(use_season, "use_season")

  if (is.na(use_season)) {
    cli::cli_abort(
      "{.arg use_season} must be a logical scalar."
    )
  }

  invisible(use_season)
}


#' Validate the marginal-likelihood control
#'
#' @inheritParams tempssm
#'
#' @return Invisibly returns \code{marginal}.
#'
#' @keywords internal
#' @noRd
.validate_marginal <- function(marginal) {
  .tempssm_check_length_one(marginal, "marginal")
  .tempssm_check_logical(marginal, "marginal")

  if (is.na(marginal)) {
    cli::cli_abort(
      "{.arg marginal} must be a logical scalar."
    )
  }

  invisible(marginal)
}


#' Check that an argument is logical
#'
#' @inheritParams .tempssm_check_numeric
#'
#' @return Invisibly returns \code{x}.
#'
#' @keywords internal
#' @noRd
.tempssm_check_logical <- function(x, arg_name, allow_null = FALSE) {
  if (is.null(x) && isTRUE(allow_null)) {
    return(invisible(x))
  }

  if (!is.logical(x)) {
    cli::cli_abort(
      "{.arg {arg_name}} must be logical."
    )
  }

  invisible(x)
}


#' Check that an argument is character
#'
#' @inheritParams .tempssm_check_numeric
#'
#' @return Invisibly returns \code{x}.
#'
#' @keywords internal
#' @noRd
.tempssm_check_character <- function(x, arg_name, allow_null = FALSE) {
  if (is.null(x) && isTRUE(allow_null)) {
    return(invisible(x))
  }

  if (!is.character(x)) {
    cli::cli_abort(
      "{.arg {arg_name}} must be character."
    )
  }

  invisible(x)
}


#' Check that an argument inherits from a class
#'
#' @param class Expected class name.
#' @inheritParams .tempssm_check_numeric
#'
#' @return Invisibly returns \code{x}.
#'
#' @keywords internal
#' @noRd
.tempssm_check_class <- function(x, arg_name, class, allow_null = FALSE) {
  if (is.null(x) && isTRUE(allow_null)) {
    return(invisible(x))
  }

  if (!inherits(x, class)) {
    cli::cli_abort(
      "{.arg {arg_name}} must be an object of class {.cls {class[[1L]]}}."
    )
  }

  invisible(x)
}


#' Check that an argument is a univariate time series
#'
#' @inheritParams .tempssm_check_length_one
#'
#' @return Invisibly returns \code{x}.
#'
#' @keywords internal
#' @noRd
.tempssm_check_univariate_ts <- function(x, arg_name) {
  .tempssm_check_class(x, arg_name, "ts")

  if (!is.null(dim(x)) && NCOL(x) != 1L) {
    cli::cli_abort(
      "{.arg {arg_name}} must be univariate."
    )
  }

  invisible(x)
}
