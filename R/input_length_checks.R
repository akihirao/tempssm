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


#' Check that an argument has an expected length
#'
#' @param x Object to check.
#' @param arg_name Name of the argument being checked.
#' @param expected_len Expected length.
#'
#' @return Invisibly returns \code{x}.
#'
#' @keywords internal
#' @noRd
.tempssm_check_length <- function(x, arg_name, expected_len) {
  if (length(x) != expected_len) {
    cli::cli_abort(
      "{.arg {arg_name}} must have length {expected_len}."
    )
  }

  invisible(x)
}


#' Check that an argument is numeric
#'
#' @param x Object to check.
#' @param arg_name Name of the argument being checked.
#' @param allow_null Logical scalar indicating whether \code{NULL} is allowed.
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
#' @param x Object to check.
#' @param arg_name Name of the argument being checked.
#'
#' @return Invisibly returns \code{x}.
#'
#' @keywords internal
#' @noRd
.tempssm_check_no_undefined <- function(x, arg_name) {
  if (any(is.nan(x)) || any(is.infinite(x))) {
    cli::cli_abort(
      "{.arg {arg_name}} must not contain {.val NaN}, {.val Inf}, or {.val -Inf}."
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
#' @param x Object to check.
#' @param arg_name Name of the argument being checked.
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
