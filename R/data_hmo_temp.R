#' Monthly mean air temperature at Hohenpeissenberg Meteorological Observatory
#'
#' A long-term monthly mean air temperature time series observed at
#' the Hohenpeissenberg Meteorological Observatory, Germany. The data are
#' are provided as a \code{ts} object and can be used as an example dataset
#' for state-space modeling of temperature time series.
#'
#' @format
#' A \code{ts} object with:
#' \describe{
#'   \item{frequency}{12 (monthly data)}
#'   \item{start}{January 1781}
#'   \item{end}{December 2025}
#' }
#'
#' @details
#' The time series represents monthly mean air temperature (in degrees
#' Celsius) at the Hohenpeissenberg Meteorological Observatory.
#' The observation period spans from January 1781 to December 2025.
#' Missing values, if any, are encoded as \code{NA}.
#'
#' @source
#' Global Historical Climatology Network monthly (GHCNm).
#' Data obtained from the NOAA website:
#' \url{https://www.ncei.noaa.gov/data/global-historical-climatology-network-monthly/v4}
#'
#' @examples
#' data(hmo_temp)
#' plot(hmo_temp)
#'
"hmo_temp"
