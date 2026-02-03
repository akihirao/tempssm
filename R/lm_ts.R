#' Fit a linear model to temperature time series
#'
#' This function estimates a linear model for monthly temperature time series.
#' Univariate \code{ts} object is supported.
#' If multivariate, a column named \code{"Temp"} is used.
#'
#' @param temp_data A temperature time series of class \code{ts}.
#'   The series can have any arbitrary frequency of 2 or higher.
#'   For example, a frequency of 12 represents a monthly \code{ts} object.
#'
#'  
#' @return An object of class \code{list} containing:
#' \describe{
#'   \item{model}{Fitted \code{SSModel} object.}
#'   \item{fit}{Result of \code{fitSSM}.}
#'   \item{kfs}{Kalman filtering and smoothing results from \code{KFS}.}
#'   \item{data}{Temperature time series used for estimation.}
#'   \item{call}{Matched function call.}
#' }
#'
#' @importFrom stats lm
#' @export
lm_ts <- function(temp_data) {

  ## ---- Input checks ---------------------------------------------------
  if (!inherits(temp_data, "ts")) {
    stop("The object of temp_data must be a 'ts' object.")
  }

  freq = frequency(temp_data)
  if (freq == 1) {
    stop("The procedure requires a ts object with frequency > 1.",
         call. = FALSE)
  }
  
  
  if(is.null(dim(temp_data))) {
    y <- as.numeric(temp_data)
    x <- time(temp_data)
  }else if(dim(temp_data)[2]!=1){ # if 'temp_data' has more than two variables
    stop("The object of temp_data must be univariant.")
  }else{
    y <- as.numeric(temp_data)
    x <- time(temp_data)
  }
  
  res <- lm(y~x)
  
  
  return(res)
}
